# Snakefile: Simulation study genotype + PGS pipeline
# ------------------------------------------------------
# Generates PLINK binary files for three cohorts (gwas, target, panel),
# simulates phenotypes, runs GWAS + LDPred2/PRS-CS, and evaluates polygenic
# scores across panel sizes and ancestries.
#
# External tools required:
#   bcftools >= 1.15, tabix, plink2 >= 2.0, wget
#   python3: msprime, stdpopsim
#   R: bigsnpr, bigsparser, optparse, data.table
#
# Quick test run (two chromosomes, two replicates, matched-admixed panel only):
#   snakemake --config chromosomes=[chr21,chr22] n_replicates=2 \
#             panel_ancestries=[matched_admixed] sim_methods=[msprime] -j4
#
# Full run:
#   snakemake all_pgs -j <n_cores> --resources mem_mb=<total_mb>

configfile: "config/simulation.yaml"

# ── Core dimensions ────────────────────────────────────────────────────────────

CHROMS      = config["chromosomes"]
N_REPS      = config["n_replicates"]
BASE_SEED   = config["base_seed"]
N_GWAS      = config["n_gwas"]
N_TARGET    = config["n_target"]
N_PANEL_MAX = config["n_panel_max"]
MAF         = config["maf_threshold"]
TRAIN_FRAC  = config.get("train_frac", 0.80)
N_PCS       = config.get("n_pcs", 10)

KG3_URL     = config.get("kg3_base_url",
                          "https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3")
KG3_EBI_URL = config.get("kg3_ebi_url",
                          "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502")

SIM_METHODS          = config.get("sim_methods", ["msprime"])
HAPNEST_CONTAINER    = config.get("hapnest_container",
                                   "resources/hapnest/intervene-synthetic-data_latest.sif")
HAPNEST_DATA_DIR     = config.get("hapnest_data_dir", "resources/hapnest/data")
HAPNEST_SUPERPOP     = config.get("hapnest_superpopulation", "AMR")
HAPNEST_USE_DOCKER   = config.get("hapnest_use_docker", False)

COHORTS = ["gwas", "target", "panel"]
REPS    = list(range(1, N_REPS + 1))

# ── Phenotype parameters ───────────────────────────────────────────────────────

H2_LEVELS       = config["h2_levels"]
P_CAUSAL_LEVELS = config["p_causal_levels"]
EFFECT_DISTS    = config["effect_dists"]
TRAIT_CONFIGS   = config["trait_configs"]
TRAIT_LABELS    = list(TRAIT_CONFIGS.keys())

# ── Panel & method parameters ──────────────────────────────────────────────────

PANEL_SIZES      = config.get("panel_sizes",
                               [100, 250, 500, 1000, 2500, 5000, 10000])
PANEL_ANCESTRIES = config.get("panel_ancestries",
                               ["matched_admixed", "EUR_1kg", "AFR_1kg",
                                "AMR_1kg", "oracle"])
PANEL_1KG_MAX    = config.get("panel_1kg_max",
                               {"EUR_1kg": 503, "AFR_1kg": 661, "AMR_1kg": 347})
PGS_METHODS      = config.get("pgs_methods", ["ldpred2", "prscs"])
KG_ANCESTRIES    = [a for a in PANEL_ANCESTRIES if a.endswith("_1kg")]


def _valid_panel_combos():
    """Return list of (panel_ancestry, panel_n) tuples that are feasible."""
    combos = []
    for n in PANEL_SIZES:
        combos.append(("matched_admixed", n))
    for anc, max_n in PANEL_1KG_MAX.items():
        if anc in PANEL_ANCESTRIES:
            for n in PANEL_SIZES:
                if n <= max_n:
                    combos.append((anc, n))
    if "oracle" in PANEL_ANCESTRIES:
        combos.append(("oracle", N_GWAS))
    return combos


VALID_PANEL_COMBOS = _valid_panel_combos()


# ── Wildcard constraints ───────────────────────────────────────────────────────

wildcard_constraints:
    sim_method     = r"msprime|hapgen2|hapnest",
    panel_ancestry = r"matched_admixed|EUR_1kg|AFR_1kg|AMR_1kg|oracle",
    panel_n        = r"\d+",
    method         = r"ldpred2|prscs",
    rep            = r"\d+",
    chrom_n        = r"\d+",


# ── Helper functions ───────────────────────────────────────────────────────────

def trait_type(label):
    return TRAIT_CONFIGS[label]["type"]


def trait_prevalence(label):
    prev = TRAIT_CONFIGS[label]["prevalence"]
    return prev if prev is not None else 0.10


def is_binary(wildcards):
    return trait_type(wildcards.trait) == "binary"


def sim_seed(wildcards):
    return BASE_SEED + int(wildcards.rep)


def chrom_num(chrom):
    """Strip 'chr' prefix: 'chr1' → '1'."""
    return str(chrom).lstrip("chr")


def panel_bed(wildcards):
    """Return the merged.bed path for the given panel_ancestry wildcard."""
    a  = wildcards.panel_ancestry
    sm = wildcards.sim_method
    if a == "matched_admixed":
        return f"results/plink/{sm}/panel/rep{wildcards.rep}/merged.bed"
    elif a == "oracle":
        return f"results/plink/{sm}/gwas/rep{wildcards.rep}/merged.bed"
    else:
        return f"resources/panels/{a}/merged.bed"


def _ref_sumstats(wildcards):
    """First available quantitative linear sumstats for the LD reference."""
    return (
        f"results/gwas/{wildcards.sim_method}/rep{wildcards.rep}/quantitative"
        f"/h2_{H2_LEVELS[0]}/pc_{P_CAUSAL_LEVELS[0]}/gaussian/sumstats.tsv"
    )


def _n_train(wildcards):
    return round(N_GWAS * TRAIN_FRAC)


# ── Target rules ───────────────────────────────────────────────────────────────

rule all:
    """Genotypes + phenotypes (no PGS)."""
    input:
        expand(
            "results/plink/{sim_method}/{cohort}/rep{rep}/merged.bed",
            sim_method=SIM_METHODS, cohort=COHORTS, rep=REPS,
        ),
        expand(
            "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}"
            "/{effect_dist}/pheno.pheno",
            sim_method=SIM_METHODS, rep=REPS, trait=TRAIT_LABELS,
            h2=H2_LEVELS, p_causal=P_CAUSAL_LEVELS, effect_dist=EFFECT_DISTS,
        ),


rule all_pgs:
    """Full PGS evaluation across all simulation methods, panel sizes, ancestries, and PGS methods."""
    input:
        [
            (
                f"results/evaluation/{method}/{sim_method}/rep{rep}"
                f"/{panel_ancestry}/n{panel_n}"
                f"/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/metrics.tsv"
            )
            for sim_method         in SIM_METHODS
            for method             in PGS_METHODS
            for rep                in REPS
            for (panel_ancestry, panel_n) in VALID_PANEL_COMBOS
            for trait              in TRAIT_LABELS
            for h2                 in H2_LEVELS
            for p_causal           in P_CAUSAL_LEVELS
            for effect_dist        in EFFECT_DISTS
        ],


rule collect_all_metrics:
    """Aggregate all per-scenario metrics.tsv files into a single table."""
    input:
        metrics = rules.all_pgs.input,
        script  = "scripts/collect_metrics.R",
    output:
        "results/summary/all_metrics.tsv",
    log: "logs/collect_metrics.log"
    shell:
        """
        Rscript {input.script} \
            --results-dir results/evaluation \
            --out         {output} \
            2> {log}
        """


# ── Step 1a: Download 1KG Phase 3 genetic maps ────────────────────────────────

rule download_genetic_map:
    """Download per-chromosome 1KG genetic map from the Oxford stats server."""
    output: "resources/1kg/genetic_map_chr{chrom_n}_combined_b37.txt"
    params: url = lambda wc: KG3_URL + "/genetic_map_chr" + wc.chrom_n + "_combined_b37.txt"
    log: "logs/download/genetic_map_chr{chrom_n}.log"
    shell:
        """
        mkdir -p resources/1kg
        wget -q -O {output} "{params.url}" 2> {log}
        """


# ── Step 1b: Simulate genotypes with msprime (admixed population) ─────────────

rule simulate_msprime:
    """Simulate genotypes using the AmericanAdmixture_4B11 demographic model."""
    input:
        script        = "scripts/simulate_genotypes.py",
        hapmap3_sites = config["hapmap3_sites_file"],
    output:
        expand(
            "results/vcf/msprime/rep{{rep}}/{cohort}_{{chrom}}.vcf.gz",
            cohort=COHORTS,
        ),
    params:
        outdir = "results/vcf/msprime/rep{rep}",
        seed   = sim_seed,
    log: "logs/simulate/msprime_rep{rep}_{chrom}.log"
    resources:
        mem_mb = 16000,
    shell:
        """
        mkdir -p {params.outdir}
        python3 {input.script} \
            --chrom          {wildcards.chrom} \
            --n-gwas         {N_GWAS} \
            --n-target       {N_TARGET} \
            --n-panel        {N_PANEL_MAX} \
            --seed           {params.seed} \
            --hapmap3-sites  {input.hapmap3_sites} \
            --outdir         {params.outdir} \
            2>> {log}
        """


# ── Step 1b-ii: Simulate genotypes with HAPNEST ───────────────────────────────

rule simulate_hapnest:
    """Simulate genotypes using HAPNEST (LD-preserving haplotype copying from 1KG+HGDP)."""
    input:
        script    = "scripts/simulate_hapnest.py",
        **({} if HAPNEST_USE_DOCKER else {"container": HAPNEST_CONTAINER}),
    output:
        expand(
            "results/vcf/hapnest/rep{{rep}}/{cohort}_{{chrom}}.vcf.gz",
            cohort=COHORTS,
        ),
    params:
        outdir      = "results/vcf/hapnest/rep{rep}",
        seed        = sim_seed,
        superpop    = HAPNEST_SUPERPOP,
        data_dir    = HAPNEST_DATA_DIR,
        container   = HAPNEST_CONTAINER,
        docker_flag = "--use-docker" if HAPNEST_USE_DOCKER else "",
    log: "logs/simulate/hapnest_rep{rep}_{chrom}.log"
    resources:
        mem_mb = 4000,
    shell:
        """
        mkdir -p {params.outdir}
        python3 {input.script} \
            --chrom              {wildcards.chrom} \
            --n-gwas             {N_GWAS} \
            --n-target           {N_TARGET} \
            --n-panel            {N_PANEL_MAX} \
            --seed               {params.seed} \
            --superpopulation    {params.superpop} \
            --outdir             {params.outdir} \
            --hapnest-container  {params.container} \
            --hapnest-data-dir   {params.data_dir} \
            --threads            {threads} \
            --memory-mb          {resources.mem_mb} \
            {params.docker_flag} \
            2>> {log}
        """


# ── Step 1c: Download 1KG Phase 3 phased VCFs from EBI FTP ───────────────────

rule download_1kg_vcf:
    """Download a 1KG Phase 3 phased VCF from the EBI FTP server."""
    output:
        vcf = "resources/1kg/ALL.{chrom}.vcf.gz",
        tbi = "resources/1kg/ALL.{chrom}.vcf.gz.tbi",
    params:
        ebi_vcf = lambda wc: KG3_EBI_URL + "/ALL." + wc.chrom + ".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        ebi_tbi = lambda wc: KG3_EBI_URL + "/ALL." + wc.chrom + ".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",
    log: "logs/download/1kg_{chrom}.log"
    shell:
        """
        mkdir -p resources/1kg
        wget -q -O {output.vcf} "{params.ebi_vcf}" 2>  {log}
        wget -q -O {output.tbi} "{params.ebi_tbi}" 2>> {log}
        """


# ── Step 1d: Write per-super-population sample ID lists ───────────────────────

rule write_1kg_sample_ids:
    """Write a single-column sample-ID list for one 1KG super-population."""
    input:
        panel  = "resources/1kg/integrated_call_samples.panel",
        script = "scripts/write_1kg_sample_ids.py",
    output:
        ids = "resources/panels/{panel_ancestry}/sample_ids.txt",
    params:
        superpop = lambda wc: wc.panel_ancestry.replace("_1kg", ""),
    wildcard_constraints:
        panel_ancestry = r"EUR_1kg|AFR_1kg|AMR_1kg",
    log: "logs/panels/write_sample_ids_{panel_ancestry}.log"
    shell:
        """
        mkdir -p $(dirname {output.ids})
        python3 {input.script} \
            --panel    {input.panel} \
            --superpop {params.superpop} \
            --out      {output.ids} \
            2> {log}
        """


# ── Step 1e: Extract 1KG super-population and convert to PLINK ───────────────

rule extract_1kg_population:
    """Subset 1KG VCFs to one super-population and produce PLINK binary files."""
    input:
        vcfs = expand("resources/1kg/ALL.{chrom}.vcf.gz", chrom=CHROMS),
        tbis = expand("resources/1kg/ALL.{chrom}.vcf.gz.tbi", chrom=CHROMS),
        ids  = "resources/panels/{panel_ancestry}/sample_ids.txt",
    output:
        bed = "resources/panels/{panel_ancestry}/merged.bed",
        bim = "resources/panels/{panel_ancestry}/merged.bim",
        fam = "resources/panels/{panel_ancestry}/merged.fam",
    params:
        prefix = "resources/panels/{panel_ancestry}/merged",
        maf    = MAF,
    wildcard_constraints:
        panel_ancestry = r"EUR_1kg|AFR_1kg|AMR_1kg",
    log: "logs/panels/extract_{panel_ancestry}.log"
    resources: mem_mb = 8000
    shell:
        """
        mkdir -p $(dirname {output.bed})
        TMP_VCF=$(mktemp --suffix=.vcf.gz)
        # Subset to super-population samples, then convert to PLINK
        bcftools concat --allow-overlaps {input.vcfs} --output-type u 2> {log} | \
        bcftools view --samples-file {input.ids} --output-type z -o "$TMP_VCF" 2>> {log}
        plink2 \
            --vcf         "$TMP_VCF" \
            --make-bed \
            --maf         {params.maf} \
            --max-alleles 2 \
            --min-alleles 2 \
            --snps-only   just-acgt \
            --output-chr  chrM \
            --out         {params.prefix} \
            2>> {log}
        rm -f "$TMP_VCF"
        """


# ── Step 2: Filter each simulated VCF ─────────────────────────────────────────

rule filter_vcf:
    """Filter to biallelic SNPs with MAF >= threshold."""
    input:
        vcf = "results/vcf/{sim_method}/rep{rep}/{cohort}_{chrom}.vcf.gz",
    output:
        vcf = "results/vcf_filtered/{sim_method}/rep{rep}/{cohort}_{chrom}.vcf.gz",
        tbi = "results/vcf_filtered/{sim_method}/rep{rep}/{cohort}_{chrom}.vcf.gz.tbi",
    params:
        maf = MAF,
    log:
        "logs/filter/{sim_method}_rep{rep}_{cohort}_{chrom}.log",
    shell:
        """
        bcftools view \
            --min-af {params.maf}:minor \
            --max-alleles 2 \
            --min-alleles 2 \
            --type snps \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


# ── Step 3: Concatenate chromosomes and convert to PLINK binary ───────────────

rule merge_and_convert:
    """Concatenate per-chromosome VCFs and convert to PLINK2 binary format."""
    input:
        vcfs = expand(
            "results/vcf_filtered/{{sim_method}}/rep{{rep}}/{{cohort}}_{chrom}.vcf.gz",
            chrom=CHROMS,
        ),
    output:
        bed = "results/plink/{sim_method}/{cohort}/rep{rep}/merged.bed",
        bim = "results/plink/{sim_method}/{cohort}/rep{rep}/merged.bim",
        fam = "results/plink/{sim_method}/{cohort}/rep{rep}/merged.fam",
    params:
        prefix = "results/plink/{sim_method}/{cohort}/rep{rep}/merged",
    log:
        "logs/merge_convert/{sim_method}_rep{rep}_{cohort}.log",
    resources:
        mem_mb = 4000,
    shell:
        """
        TMP_VCF=$(mktemp --suffix=.vcf.gz)
        bcftools concat --allow-overlaps {input.vcfs} --output-type z -o "$TMP_VCF" 2> {log}
        plink2 \
            --vcf "$TMP_VCF" \
            --make-bed \
            --out {params.prefix} \
            2>> {log}
        rm -f "$TMP_VCF"
        """


# ── Step 4: Simulate phenotypes ───────────────────────────────────────────────

rule simulate_phenotypes:
    """Simulate phenotypes from GWAS population genotypes."""
    input:
        bed    = "results/plink/{sim_method}/gwas/rep{rep}/merged.bed",
        script = "scripts/simulate_phenotypes.R",
    output:
        pheno      = "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/pheno.pheno",
        causal     = "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/pheno.causal_snps.tsv",
        params_tsv = "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/pheno.sim_params.tsv",
    params:
        out_prefix = "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/pheno",
        trait_type  = lambda wc: trait_type(wc.trait),
        prevalence  = lambda wc: trait_prevalence(wc.trait),
        seed        = lambda wc: BASE_SEED + int(wc.rep),
    log:
        "logs/phenotypes/{sim_method}_rep{rep}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    resources:
        mem_mb = 4000,
    shell:
        """
        Rscript {input.script} \
            --bed         {input.bed} \
            --h2          {wildcards.h2} \
            --p-causal    {wildcards.p_causal} \
            --effect-dist {wildcards.effect_dist} \
            --trait-type  {params.trait_type} \
            --prevalence  {params.prevalence} \
            --seed        {params.seed} \
            --out         {params.out_prefix} \
            2> {log}
        """


# ── Step 5: Split GWAS individuals into training and test sets ────────────────

rule split_subjects:
    """Randomly split GWAS individuals into training and held-out test sets."""
    input:
        fam    = "results/plink/{sim_method}/gwas/rep{rep}/merged.fam",
        script = "scripts/split_subjects.R",
    output:
        train = "results/splits/{sim_method}/rep{rep}/train.txt",
        test  = "results/splits/{sim_method}/rep{rep}/test.txt",
    params:
        train_frac = TRAIN_FRAC,
        seed       = sim_seed,
    log:
        "logs/splits/{sim_method}_rep{rep}.log",
    shell:
        """
        Rscript {input.script} \
            --fam        {input.fam} \
            --train-frac {params.train_frac} \
            --seed       {params.seed} \
            --train      {output.train} \
            --test       {output.test} \
            2> {log}
        """


# ── Step 6a: Compute GWAS training-set PCA ────────────────────────────────────

rule run_pca_gwas:
    """Compute top PCs on the GWAS training set for use as GWAS covariates."""
    input:
        bed       = "results/plink/{sim_method}/gwas/rep{rep}/merged.bed",
        train_ids = "results/splits/{sim_method}/rep{rep}/train.txt",
    output:
        eigenvec = "results/pca/{sim_method}/gwas/rep{rep}/pcs.eigenvec",
        eigenval = "results/pca/{sim_method}/gwas/rep{rep}/pcs.eigenval",
    params:
        prefix     = "results/pca/{sim_method}/gwas/rep{rep}/pcs",
        bed_prefix = lambda wc, input: input.bed[:-4],
        n_pcs      = N_PCS,
    log: "logs/pca/{sim_method}_gwas_rep{rep}.log"
    resources: mem_mb = 4000
    shell:
        """
        mkdir -p $(dirname {output.eigenvec})
        plink2 \
            --bfile {params.bed_prefix} \
            --keep  {input.train_ids} \
            --pca   {params.n_pcs} \
            --out   {params.prefix} \
            2> {log}
        """


# ── Step 6b: Project test-set individuals onto training PCs ───────────────────

rule project_pca_test:
    """Project held-out test individuals onto training-set PC axes."""
    input:
        bed       = "results/plink/{sim_method}/gwas/rep{rep}/merged.bed",
        test_ids  = "results/splits/{sim_method}/rep{rep}/test.txt",
        eigenvec  = "results/pca/{sim_method}/gwas/rep{rep}/pcs.eigenvec",
    output:
        scores = "results/pca/{sim_method}/test/rep{rep}/pcs.sscore",
    params:
        bed_prefix   = lambda wc, input: input.bed[:-4],
        score_prefix = "results/pca/{sim_method}/test/rep{rep}/pcs",
    log: "logs/pca/{sim_method}_test_rep{rep}.log"
    resources: mem_mb = 4000
    shell:
        """
        mkdir -p $(dirname {output.scores})
        plink2 \
            --bfile  {params.bed_prefix} \
            --keep   {input.test_ids} \
            --score  {input.eigenvec} 1 2 header \
            --out    {params.score_prefix} \
            2> {log}
        """


# ── Step 7: Run GWAS on training subjects ─────────────────────────────────────

rule run_gwas:
    """Run GWAS on training-set individuals with plink2 --glm, adjusted for PCs."""
    input:
        bed       = "results/plink/{sim_method}/gwas/rep{rep}/merged.bed",
        pheno     = "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/pheno.pheno",
        train_ids = "results/splits/{sim_method}/rep{rep}/train.txt",
        pcs       = "results/pca/{sim_method}/gwas/rep{rep}/pcs.eigenvec",
    output:
        sumstats = "results/gwas/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/sumstats.tsv",
    params:
        bed_prefix  = lambda wc, input: input.bed[:-4],
        gwas_prefix = "results/gwas/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/gwas",
        binary_flag = lambda wc: "--1" if is_binary(wc) else "",
    log:
        "logs/gwas/{sim_method}_rep{rep}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    resources:
        mem_mb = 4000,
    shell:
        """
        mkdir -p $(dirname {output.sumstats})
        plink2 \
            --bfile   {params.bed_prefix} \
            --pheno   {input.pheno} \
            --pheno-name Y \
            --keep    {input.train_ids} \
            --covar   {input.pcs} \
            --covar-variance-standardize \
            --glm     no-x-sex hide-covar \
            {params.binary_flag} \
            --no-psam-pheno \
            --threads 1 \
            --out     {params.gwas_prefix} \
            2> {log}
        # Normalise output filename regardless of trait type
        OUTFILE=$(ls {params.gwas_prefix}.Y.glm.* 2>/dev/null | head -1)
        if [ -z "$OUTFILE" ]; then
            echo "ERROR: plink2 --glm produced no output" >> {log}; exit 1
        fi
        mv "$OUTFILE" {output.sumstats}
        """


# ── Step 8a: Precompute shared LD reference (panel FBM + SFBM) ────────────────
# Runs once per (sim_method, rep, panel_ancestry, panel_n); all run_ldpred2 jobs
# for that combination reuse the output.

rule prepare_ldpred2_ref:
    """Precompute LD reference panel (SFBM) for a given ancestry and panel size."""
    input:
        bed          = lambda wc: panel_bed(wc),
        bim          = lambda wc: panel_bed(wc)[:-4] + ".bim",
        fam          = lambda wc: panel_bed(wc)[:-4] + ".fam",
        ref_sumstats = _ref_sumstats,
        script       = "scripts/prepare_ldpred2_ref.R",
    output:
        sfbm         = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ld_sfbm.sbk",
        sfbm_rds     = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ld_sfbm.rds",
        matched_snps = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/matched_snps.tsv",
    params:
        out_dir  = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}",
        window   = 3000000,
        ncores   = 6,
        n_panel  = lambda wc: wc.panel_n,
        seed     = sim_seed,
    threads: 6
    log: "logs/ldpred2_ref/{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}.log"
    resources: mem_mb = 32000
    shell:
        """
        mkdir -p {params.out_dir}
        Rscript {input.script} \
            --panel-bed    {input.bed} \
            --ref-sumstats {input.ref_sumstats} \
            --out-dir      {params.out_dir} \
            --window       {params.window} \
            --ncores       {params.ncores} \
            --n-panel      {params.n_panel} \
            --seed         {params.seed} \
            2> {log}
        """


# ── Step 8b: Compute LDPred2-auto weights ─────────────────────────────────────

rule run_ldpred2:
    """Compute LDPred2-auto polygenic score weights using shared LD reference."""
    input:
        sumstats  = "results/gwas/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/sumstats.tsv",
        ref_sfbm  = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ld_sfbm.rds",
        ref_snps  = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/matched_snps.tsv",
        pheno     = "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/pheno.pheno",
        train_ids = "results/splits/{sim_method}/rep{rep}/train.txt",
        script    = "scripts/run_ldpred2.R",
    output:
        betas = "results/pgs_weights/ldpred2/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/betas.tsv",
    params:
        trait_type = lambda wc: trait_type(wc.trait),
        n_train    = _n_train,
        ref_dir    = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}",
        seed       = lambda wc: BASE_SEED + int(wc.rep),
    log:
        "logs/ldpred2/{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    resources:
        mem_mb = 4000,
    shell:
        """
        Rscript {input.script} \
            --sumstats   {input.sumstats} \
            --ref-dir    {params.ref_dir} \
            --pheno      {input.pheno} \
            --train-ids  {input.train_ids} \
            --trait-type {params.trait_type} \
            --n-train    {params.n_train} \
            --seed       {params.seed} \
            --out        {output.betas} \
            2> {log}
        """


# ── Step 8c: Compute PRS-CS weights ───────────────────────────────────────────

rule run_prscs:
    """Compute PRS-CS posterior effect-size weights."""
    input:
        sumstats = "results/gwas/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/sumstats.tsv",
        script   = "scripts/run_prscs.py",
    output:
        betas = "results/pgs_weights/prscs/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/betas.tsv",
    params:
        ref_dir = lambda wc: "resources/prscs_ref/" + (wc.panel_ancestry if wc.panel_ancestry != "oracle" else "matched_admixed"),
        n_train = _n_train,
        seed    = lambda wc: BASE_SEED + int(wc.rep),
    log:
        "logs/prscs/{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    resources:
        mem_mb = 4000,
    shell:
        """
        python3 {input.script} \
            --sumstats {input.sumstats} \
            --ref-dir  {params.ref_dir} \
            --n-gwas   {params.n_train} \
            --seed     {params.seed} \
            --out      {output.betas} \
            2> {log}
        """


# ── Step 9: Score test-set individuals ────────────────────────────────────────

rule score_test_set:
    """Apply PGS weights to held-out test individuals with plink2 --score."""
    input:
        bed      = "results/plink/{sim_method}/gwas/rep{rep}/merged.bed",
        betas    = "results/pgs_weights/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/betas.tsv",
        test_ids = "results/splits/{sim_method}/rep{rep}/test.txt",
    output:
        sscore = "results/pgs/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/scores.sscore",
    params:
        bed_prefix   = lambda wc, input: input.bed[:-4],
        score_prefix = "results/pgs/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/scores",
    log:
        "logs/score/{method}_{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    resources:
        mem_mb = 4000,
    shell:
        """
        mkdir -p $(dirname {output.sscore})
        plink2 \
            --bfile   {params.bed_prefix} \
            --keep    {input.test_ids} \
            --score   {input.betas} 1 2 3 header \
            --threads 1 \
            --out     {params.score_prefix} \
            2> {log}
        """


# ── Step 10: Evaluate prediction accuracy ─────────────────────────────────────

rule evaluate_pgs:
    """Compute R² (quantitative) or AUC (binary) for the test-set PGS."""
    input:
        scores   = "results/pgs/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/scores.sscore",
        pheno    = "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/pheno.pheno",
        test_ids = "results/splits/{sim_method}/rep{rep}/test.txt",
        pcs      = "results/pca/{sim_method}/test/rep{rep}/pcs.sscore",
        script   = "scripts/evaluate_pgs.R",
    output:
        metrics = "results/evaluation/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/metrics.tsv",
    params:
        trait_type = lambda wc: trait_type(wc.trait),
    log:
        "logs/evaluation/{method}_{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    shell:
        """
        Rscript {input.script} \
            --scores     {input.scores} \
            --pheno      {input.pheno} \
            --test-ids   {input.test_ids} \
            --covariates {input.pcs} \
            --trait-type {params.trait_type} \
            --out        {output.metrics} \
            2> {log}
        """
