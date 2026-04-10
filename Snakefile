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
KG38_EBI_URL = config.get(
    "kg38_ebi_url",
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/release/20181203_biallelic_SNV",
)

SIM_METHODS          = config.get("sim_methods", ["hapnest"])
HAPNEST_CONTAINER    = config.get("hapnest_container",
                                   "resources/hapnest/intervene-synthetic-data_latest.sif")
HAPNEST_DATA_DIR     = config.get("hapnest_data_dir", "resources/hapnest/data")
HAPNEST_SUPERPOP     = config.get("hapnest_superpopulation", "AMR")
HAPNEST_USE_DOCKER   = config.get("hapnest_use_docker", False)
HAPNEST_REPO         = config.get("hapnest_repo", "resources/hapnest/repo")

COHORTS = ["gwas", "target", "panel"]
REPS    = list(range(1, N_REPS + 1))

# ── Phenotype parameters ───────────────────────────────────────────────────────

H2_LEVELS       = config["h2_levels"]
P_CAUSAL_LEVELS = config["p_causal_levels"]
EFFECT_DISTS    = config["effect_dists"]
TRAIT_CONFIGS   = config["trait_configs"]
TRAIT_LABELS    = list(TRAIT_CONFIGS.keys())

# ── Panel & method parameters ──────────────────────────────────────────────────

# PRS-CS: maps each panel ancestry to its 1KG population code for ref download.
PRSCS_POP = {
    "matched_admixed": "amr",
    "EUR_1kg":         "eur",
    "AFR_1kg":         "afr",
    "AMR_1kg":         "amr",
    "oracle":          "amr",
    "gwas_subset":     "amr",
}
# PRScs.py checks os.path.basename(ref_dir) for '1kg' or 'ukbb' before loading
# the reference — directories that don't match either string leave ref_dict
# unbound and cause an UnboundLocalError.  Map every panel_ancestry to a ref
# directory whose basename contains '1kg'.
PRSCS_REF_DIR = {
    "matched_admixed": "AMR_1kg",
    "EUR_1kg":         "EUR_1kg",
    "AFR_1kg":         "AFR_1kg",
    "AMR_1kg":         "AMR_1kg",
    "oracle":          "AMR_1kg",
    "gwas_subset":     "AMR_1kg",
}
PRSCS_REF_URLS = config.get("prscs_ref_urls", {})

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
    if "gwas_subset" in PANEL_ANCESTRIES:
        for n in PANEL_SIZES:
            if n < N_GWAS:
                combos.append(("gwas_subset", n))
    return combos


VALID_PANEL_COMBOS = _valid_panel_combos()


# ── Wildcard constraints ───────────────────────────────────────────────────────

wildcard_constraints:
    sim_method     = r"msprime|hapgen2|hapnest|hapnest_public",
    panel_ancestry = r"matched_admixed|EUR_1kg|AFR_1kg|AMR_1kg|oracle|gwas_subset",
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
    elif a in ("oracle", "gwas_subset"):
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


rule all_hapnest:
    """All HAPNEST public-dataset VCF outputs (gwas/target/panel) for every chromosome."""
    input:
        expand(
            "results/vcf/hapnest_public/rep1/{cohort}_{chrom}.vcf.gz",
            cohort=COHORTS, chrom=CHROMS,
        ),


ALL_METRICS_FILES = [
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
]


rule all_pgs:
    """Full PGS evaluation across all simulation methods, panel sizes, ancestries, and PGS methods."""
    input:
        ALL_METRICS_FILES,
        "results/summary/all_metrics.tsv",


rule collect_all_metrics:
    """Aggregate all per-scenario metrics.tsv files into a single table."""
    input:
        metrics = ALL_METRICS_FILES,
        script  = "scripts/collect_metrics.R",
    output:
        "results/summary/all_metrics.tsv",
    log: "logs/collect_metrics.log"
    shell:
        """
        module load R/4.4.3-gcc-11.2.0-mkl
        Rscript {input.script} \
            --results-dir results/evaluation \
            --out         {output} \
            2> {log}
        """


# ── Step 1a-i: Download HapMap3 SNP positions ────────────────────────────────

rule download_hm3_ldpred2_rds:
    """Download the full LDpred2 HapMap3 map RDS from figshare.
    This file includes pos_hg38 alongside GRCh37 pos, used by extract_hm3_grch38_sites.
    Figshare item 13034123, file 25503788 (Privé et al. LDpred2 paper).
    """
    output:
        rds = "resources/map_hm3_ldpred2.rds",
    log: "logs/download/map_hm3_ldpred2.log"
    shell:
        """
        mkdir -p resources
        wget -q -O {output.rds} "https://ndownloader.figshare.com/files/25503788" 2> {log}
        """


rule download_hapmap3_sites:
    """Download HapMap3 Phase 3 SNP positions (GRCh37) used by simulate_msprime.
    Source: Privé et al. LDpred2 tutorial data (figshare).
    Note: only needed when sim_methods includes msprime.
    """
    output:
        tsv = config["hapmap3_sites_file"],
    log: "logs/download/hapmap3_sites.log"
    shell:
        """
        module load R/4.4.3-gcc-11.2.0-mkl
        mkdir -p $(dirname {output.tsv})
        Rscript -e "
          tmp <- tempfile(fileext = '.rds')
          download.file('https://figshare.com/ndownloader/files/36360325', tmp, quiet = TRUE)
          map <- readRDS(tmp)
          write.table(data.frame(chrom = map[['chr']], pos = map[['pos']]),
                      '{output.tsv}', sep = '\\t', quote = FALSE, row.names = FALSE)
          unlink(tmp)
        " 2> {log}
        """


rule extract_hm3_grch38_sites:
    """Extract HapMap3 SNP positions on GRCh38 from the LDpred2 map RDS.

    The map_hm3_ldpred2.rds published by Privé et al. includes a pos_hg38 column
    alongside the GRCh37 pos column. This rule writes chrom, pos_hg38, rsid, a0,
    a1 to a TSV used by run_prscs.py when processing GRCh38 data (hapnest_public).
    """
    input:
        rds = "resources/map_hm3_ldpred2.rds",
    output:
        tsv = "resources/hapmap3_sites_grch38.tsv",
    log: "logs/setup/extract_hm3_grch38.log"
    shell:
        """
        module load R/4.4.3-gcc-11.2.0-mkl
        Rscript -e "
          info <- readRDS('{input.rds}')
          out  <- data.frame(chrom = info[['chr']],
                             pos   = as.integer(info[['pos_hg38']]),
                             rsid  = info[['rsid']],
                             a0    = info[['a0']],
                             a1    = info[['a1']])
          out <- out[!is.na(out[['pos']]) & out[['pos']] > 0, ]
          out <- out[order(out[['chrom']], out[['pos']]), ]
          write.table(out, '{output.tsv}', sep='\\t', row.names=FALSE, quote=FALSE)
          cat(sprintf('[extract_hm3_grch38] wrote %d SNPs\\n', nrow(out)))
        " > {log} 2>&1
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


# ── Step 1b: Download HAPNEST public dataset (EBI BioStudies S-BSST936) ───────

EBI_HAPNEST_BASE = "https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/936/S-BSST936/Files"


rule download_hapnest_public_sample:
    """Download the population-label file for the HAPNEST public dataset."""
    output:
        sample = "resources/hapnest_public/synthetic_v1.sample",
    log: "logs/setup/download_hapnest_public_sample.log"
    shell:
        """
        mkdir -p resources/hapnest_public
        wget -q -O {output.sample} \
            {EBI_HAPNEST_BASE}/synthetic_v1.sample \
            2> {log}
        """


rule download_hapnest_public_plink:
    """Download PLINK binary files for one chromosome from the HAPNEST public dataset."""
    output:
        bed = "resources/hapnest_public/raw/synthetic_v1_chr-{chrom_n}.bed",
        bim = "resources/hapnest_public/raw/synthetic_v1_chr-{chrom_n}.bim",
        fam = "resources/hapnest_public/raw/synthetic_v1_chr-{chrom_n}.fam",
    log: "logs/setup/download_hapnest_public_chr{chrom_n}.log"
    resources:
        runtime = 480,   # large files — chr1 bed is ~125 GB
    shell:
        """
        mkdir -p resources/hapnest_public/raw
        wget -q -O {output.bed} \
            {EBI_HAPNEST_BASE}/genotypes/synthetic_v1_chr-{wildcards.chrom_n}.bed \
            2>  {log}
        wget -q -O {output.bim} \
            {EBI_HAPNEST_BASE}/genotypes/synthetic_v1_chr-{wildcards.chrom_n}.bim \
            2>> {log}
        wget -q -O {output.fam} \
            {EBI_HAPNEST_BASE}/genotypes/synthetic_v1_chr-{wildcards.chrom_n}.fam \
            2>> {log}
        """


rule partition_hapnest_public:
    """Filter to AMR subjects, randomly select 30 k, and split into cohorts.

    Reads the chr1 FAM file (sample IDs) and the population sample file,
    keeps only AMR individuals, draws n_gwas+n_target+n_panel at random,
    and writes per-cohort keep files (for plink2 --keep) and rename files
    (for bcftools reheader -s).  Run once; independent of chromosome.
    """
    input:
        fam    = "resources/hapnest_public/raw/synthetic_v1_chr-{}.fam".format(CHROMS[0].lstrip("chr")),
        sample = "resources/hapnest_public/synthetic_v1.sample",
        script = "scripts/partition_hapnest_public.py",
    output:
        keep_gwas   = "resources/hapnest_public/keep_gwas.txt",
        keep_target = "resources/hapnest_public/keep_target.txt",
        keep_panel  = "resources/hapnest_public/keep_panel.txt",
        rename_gwas   = "resources/hapnest_public/rename_gwas.txt",
        rename_target = "resources/hapnest_public/rename_target.txt",
        rename_panel  = "resources/hapnest_public/rename_panel.txt",
    log: "logs/setup/partition_hapnest_public.log"
    shell:
        """
        python3 {input.script} \
            --fam      {input.fam} \
            --sample   {input.sample} \
            --n-gwas   {N_GWAS} \
            --n-target {N_TARGET} \
            --n-panel  {N_PANEL_MAX} \
            --seed     {BASE_SEED} \
            --outdir   resources/hapnest_public \
            2> {log}
        """


rule extract_hapnest_public_vcf:
    """Extract one cohort from the public PLINK files and export as bgzipped VCF.

    After extraction and sample renaming, bcftools annotate adds rsIDs from the
    1KG GRCh38 lookup table so downstream tools (LDpred2, PRS-CS) can match
    variants by rsID rather than position.
    """
    input:
        bed      = lambda wc: "resources/hapnest_public/raw/synthetic_v1_chr-{}.bed".format(wc.chrom.lstrip("chr")),
        bim      = lambda wc: "resources/hapnest_public/raw/synthetic_v1_chr-{}.bim".format(wc.chrom.lstrip("chr")),
        fam      = lambda wc: "resources/hapnest_public/raw/synthetic_v1_chr-{}.fam".format(wc.chrom.lstrip("chr")),
        keep     = "resources/hapnest_public/keep_{cohort}.txt",
        rename   = "resources/hapnest_public/rename_{cohort}.txt",
        rsid_map = "resources/1kg_grch38/rsid_map_{chrom}.tsv.gz",
        rsid_tbi = "resources/1kg_grch38/rsid_map_{chrom}.tsv.gz.tbi",
    output:
        vcf = "results/vcf/hapnest_public/rep1/{cohort}_{chrom}.vcf.gz",
    params:
        plink_prefix = lambda wc: "resources/hapnest_public/raw/synthetic_v1_chr-{}".format(wc.chrom.lstrip("chr")),
        tmp_prefix   = "results/vcf/hapnest_public/rep1/_tmp_{cohort}_{chrom}",
    log: "logs/extract/hapnest_public_{cohort}_{chrom}.log"
    threads: 4
    resources:
        mem_mb  = 16000,
        runtime = 120,
    shell:
        """
        module load plink/2.0-alpha
        module load bcftools/1.19
        mkdir -p results/vcf/hapnest_public/rep1
        plink2 \
            --bfile   {params.plink_prefix} \
            --keep    {input.keep} \
            --export  vcf-4.2 bgz \
            --threads {threads} \
            --memory  {resources.mem_mb} \
            --out     {params.tmp_prefix} \
            2>> {log}
        bcftools reheader \
            --samples {input.rename} \
            {params.tmp_prefix}.vcf.gz \
        | bcftools annotate \
            --annotations {input.rsid_map} \
            --columns CHROM,POS,ID,REF,ALT \
            --output-type z \
            --output {output.vcf} \
            2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        rm -f {params.tmp_prefix}.vcf.gz {params.tmp_prefix}.log
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


# ── Step 1c-ii: Download 1KG GRCh38 VCFs (for hapnest_public rsID annotation) ─

rule download_1kg_grch38_vcf:
    """Download 1KG Phase 3 GRCh38 phased VCF (biallelic SNVs) for rsID annotation."""
    output:
        vcf = "resources/1kg_grch38/ALL.{chrom}.vcf.gz",
        tbi = "resources/1kg_grch38/ALL.{chrom}.vcf.gz.tbi",
    params:
        vcf_url = lambda wc: (
            KG38_EBI_URL + "/ALL." + wc.chrom +
            ".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
        ),
        tbi_url = lambda wc: (
            KG38_EBI_URL + "/ALL." + wc.chrom +
            ".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi"
        ),
    log: "logs/download/1kg_grch38_{chrom}.log"
    resources:
        runtime = 480,
    shell:
        """
        mkdir -p resources/1kg_grch38
        wget -q -O {output.vcf} "{params.vcf_url}" 2>  {log}
        wget -q -O {output.tbi} "{params.tbi_url}" 2>> {log}
        """


rule make_grch38_rsid_map:
    """Extract a (CHROM, POS, ID, REF, ALT) table from a 1KG GRCh38 VCF.

    Used by run_prscs.py to map hapnest_public GRCh38 positions to rsIDs before
    matching against the PRS-CS snpinfo_1kg_hm3 (GRCh37) by rsID.
    Output: bgzipped TSV, no header, columns: CHROM(no chr) POS ID REF ALT.
    """
    input:
        vcf = "resources/1kg_grch38/ALL.{chrom}.vcf.gz",
        tbi = "resources/1kg_grch38/ALL.{chrom}.vcf.gz.tbi",
    output:
        tsv = "resources/1kg_grch38/rsid_map_{chrom}.tsv.gz",
        tbi = "resources/1kg_grch38/rsid_map_{chrom}.tsv.gz.tbi",
    log: "logs/setup/grch38_rsid_map_{chrom}.log"
    shell:
        """
        module load bcftools/1.19
        bcftools query \
            --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
            {input.vcf} \
        | awk '{{gsub(/^chr/, "", $1); print}}' \
        | bgzip \
        > {output.tsv} \
        2> {log}
        tabix -s1 -b2 -e2 {output.tsv} 2>> {log}
        """


# ── Step 1c-ii: Download 1KG sample panel file ────────────────────────────────

rule download_1kg_panel:
    """Download the 1KG Phase 3 integrated sample panel file (population assignments)."""
    output:
        panel = "resources/1kg/integrated_call_samples.panel",
    params:
        url = KG3_EBI_URL + "/integrated_call_samples_v3.20130502.ALL.panel",
    log: "logs/download/1kg_panel.log"
    shell:
        """
        mkdir -p resources/1kg
        wget -q -O {output.panel} "{params.url}" 2> {log}
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
    """Subset 1KG GRCh38 VCFs to one super-population and produce PLINK binary files."""
    input:
        vcfs = expand("resources/1kg_grch38/ALL.{chrom}.vcf.gz", chrom=CHROMS),
        tbis = expand("resources/1kg_grch38/ALL.{chrom}.vcf.gz.tbi", chrom=CHROMS),
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
    resources: mem_mb = 16000
    shell:
        """
        module load bcftools/1.19 plink/2.0-alpha
        mkdir -p $(dirname {output.bed})
        TMP_VCF=$(mktemp --tmpdir=$(dirname {output.bed}) --suffix=.vcf.gz)
        # Subset to super-population samples, then convert to PLINK
        bcftools concat --allow-overlaps {input.vcfs} --output-type u 2> {log} | \
        bcftools view --samples-file {input.ids} --force-samples --output-type z -o "$TMP_VCF" 2>> {log}
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
        module load bcftools/1.19
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
        mem_mb = 16000,
    shell:
        """
        module load bcftools/1.19 plink/2.0-alpha
        TMP_VCF=$(mktemp --suffix=.vcf.gz)
        bcftools concat --allow-overlaps {input.vcfs} --output-type z -o "$TMP_VCF" 2> {log}
        plink2 \
            --vcf "$TMP_VCF" \
            --make-bed \
            --memory {resources.mem_mb} \
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
        mem_mb = 24000,
    shell:
        """
        module load R/4.4.3-gcc-11.2.0-mkl
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
        module load R/4.4.3-gcc-11.2.0-mkl
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
        eigenvec       = "results/pca/{sim_method}/gwas/rep{rep}/pcs.eigenvec",
        eigenval       = "results/pca/{sim_method}/gwas/rep{rep}/pcs.eigenval",
        eigenvec_allele = "results/pca/{sim_method}/gwas/rep{rep}/pcs.eigenvec.allele",
    params:
        prefix     = "results/pca/{sim_method}/gwas/rep{rep}/pcs",
        bed_prefix = lambda wc, input: input.bed[:-4],
        n_pcs      = N_PCS,
    log: "logs/pca/{sim_method}_gwas_rep{rep}.log"
    resources: mem_mb = 4000
    shell:
        """
        module load plink/2.0-alpha
        mkdir -p $(dirname {output.eigenvec})
        plink2 \
            --bfile {params.bed_prefix} \
            --keep  {input.train_ids} \
            --pca   {params.n_pcs} allele-wts \
            --out   {params.prefix} \
            2> {log}
        """


# ── Step 6b: Project test-set individuals onto training PCs ───────────────────

rule project_pca_test:
    """Project held-out test individuals onto training-set PC axes."""
    input:
        bed             = "results/plink/{sim_method}/gwas/rep{rep}/merged.bed",
        test_ids        = "results/splits/{sim_method}/rep{rep}/test.txt",
        eigenvec_allele = "results/pca/{sim_method}/gwas/rep{rep}/pcs.eigenvec.allele",
    output:
        scores = "results/pca/{sim_method}/test/rep{rep}/pcs.sscore",
    params:
        bed_prefix     = lambda wc, input: input.bed[:-4],
        score_prefix   = "results/pca/{sim_method}/test/rep{rep}/pcs",
        score_col_end  = 6 + N_PCS,
    log: "logs/pca/{sim_method}_test_rep{rep}.log"
    resources: mem_mb = 4000
    shell:
        """
        module load plink/2.0-alpha
        mkdir -p $(dirname {output.scores})
        plink2 \
            --bfile          {params.bed_prefix} \
            --keep           {input.test_ids} \
            --score          {input.eigenvec_allele} 2 6 header-read \
            --score-col-nums 7-{params.score_col_end} \
            --out            {params.score_prefix} \
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
        module load plink/2.0-alpha
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
        hm3_pos      = lambda wc: (
            "resources/hapmap3_sites_grch38.tsv"
            if wc.sim_method == "hapnest_public"
            else "resources/map_hm3_ldpred2.rds"
        ),
    output:
        sfbm         = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ld_sfbm.sbk",
        sfbm_rds     = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ld_sfbm.rds",
        matched_snps = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/matched_snps.tsv",
    params:
        out_dir      = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}",
        window       = 1000000,
        ncores       = 4,
        n_panel      = lambda wc: wc.panel_n,
        seed         = sim_seed,
        hm3_pos_arg  = lambda wc, input: (
            f"--hm3-positions {input.hm3_pos}"
            if wc.sim_method == "hapnest_public"
            else ""
        ),
    threads: 4
    log: "logs/ldpred2_ref/{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}.log"
    resources: mem_mb = 32000
    shell:
        """
        module load R/4.4.3-gcc-11.2.0-mkl
        mkdir -p {params.out_dir}
        Rscript {input.script} \
            --panel-bed    {input.bed} \
            --ref-sumstats {input.ref_sumstats} \
            --out-dir      {params.out_dir} \
            --window       {params.window} \
            --ncores       {params.ncores} \
            --n-panel      {params.n_panel} \
            --seed         {params.seed} \
            {params.hm3_pos_arg} \
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
        mem_mb = 24000,
    shell:
        """
        module load R/4.4.3-gcc-11.2.0-mkl
        Rscript {input.script} \
            --sumstats   {input.sumstats} \
            --ref-dir    {params.ref_dir} \
            --pheno      {input.pheno} \
            --train-ids  {input.train_ids} \
            --trait-type {params.trait_type} \
            --n-train    {params.n_train} \
            --seed       {params.seed} \
            --out        {output.betas} \
            > {log} 2>&1
        """


# ── Step 8c-ii: Install PRS-CS software ───────────────────────────────────────

rule install_prscs:
    """Clone the PRS-CS repository so PRScs.py is available locally."""
    output:
        exe = "resources/prscs/PRScs.py",
    log: "logs/setup/install_prscs.log"
    shell:
        """
        rm -rf resources/prscs
        mkdir -p resources/prscs
        wget -q -O - "https://github.com/getian107/PRScs/archive/refs/heads/master.tar.gz" \
            | tar -xz --strip-components=1 -C resources/prscs \
            2> {log}
        """


# ── Step 8c-i: Download PRS-CS 1KG LD reference panel ────────────────────────

rule download_prscs_ref:
    """Download and extract a PRS-CS 1KG LD reference panel for one panel ancestry.

    The tarball from getian107/PRScs contains ldblk_1kg_chr{N}.hdf5 files and
    snpinfo_1kg_hm3, all in a single top-level directory which is stripped so
    they land directly in resources/prscs_ref/{panel_ancestry}/.
    """
    output:
        sentinel = "resources/prscs_ref/{panel_ancestry}/snpinfo_1kg_hm3",
    params:
        outdir = "resources/prscs_ref/{panel_ancestry}",
        url    = lambda wc: PRSCS_REF_URLS[PRSCS_POP[wc.panel_ancestry]],
    wildcard_constraints:
        panel_ancestry = r"matched_admixed|EUR_1kg|AFR_1kg|AMR_1kg",
    log: "logs/setup/prscs_ref_{panel_ancestry}.log"
    resources:
        runtime = 120,
    shell:
        """
        mkdir -p {params.outdir}
        wget -q -O - "{params.url}" \
            | tar -xz --strip-components=1 -C {params.outdir} \
            2> {log}
        """


# ── Step 8c-i: Build custom PRS-CS LD reference from panel genotypes ──────────

rule build_prscs_ref_custom:
    """Build a PRS-CS HDF5 LD reference from the simulated/1KG panel BED file.

    LD block boundaries are found via the MDL criterion (Berisa & Pickrell 2016).
    The output subdirectory is named 'ref_1kg' so PRScs.py parse_ldblk() accepts
    it.  Runs once per (sim_method, rep, panel_ancestry, panel_n); all run_prscs
    jobs for that combination share the output.
    """
    input:
        bed       = lambda wc: panel_bed(wc),
        bim       = lambda wc: panel_bed(wc)[:-4] + ".bim",
        fam       = lambda wc: panel_bed(wc)[:-4] + ".fam",
        script    = "scripts/make_prscs_ref.py",
        hm3_grch38 = lambda wc: "resources/hapmap3_sites_grch38.tsv"
                                 if wc.sim_method == "hapnest_public" else [],
    output:
        sentinel = "results/prscs_custom_ref/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ref_1kg/snpinfo_1kg_hm3",
    params:
        bfile          = lambda wc, input: input.bed[:-4],
        out_dir        = "results/prscs_custom_ref/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ref_1kg",
        chroms         = " ".join(chrom_num(c) for c in CHROMS),
        n_panel        = lambda wc: wc.panel_n,
        seed           = sim_seed,
        max_block_size = 400,
        hm3_grch38_arg = lambda wc, input: f"--hm3-grch38 {input.hm3_grch38}"
                                            if wc.sim_method == "hapnest_public" else "",
    log: "logs/prscs_ref_custom/{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}.log"
    resources:
        mem_mb = 8000,
    shell:
        """
        module load python/3.12
        mkdir -p {params.out_dir}
        python3 {input.script} \
            --bfile          {params.bfile} \
            --out            {params.out_dir} \
            --chroms         {params.chroms} \
            --n-panel        {params.n_panel} \
            --seed           {params.seed} \
            --max-block-size {params.max_block_size} \
            {params.hm3_grch38_arg} \
            2> {log}
        """


# ── Step 8c: Compute PRS-CS weights ───────────────────────────────────────────

rule run_prscs:
    """Compute PRS-CS posterior effect-size weights."""
    input:
        sumstats   = "results/gwas/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/sumstats.tsv",
        bim        = "results/plink/{sim_method}/gwas/rep{rep}/merged.bim",
        script     = "scripts/run_prscs.py",
        prscs_exe  = "resources/prscs/PRScs.py",
        ref_data   = "results/prscs_custom_ref/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ref_1kg/snpinfo_1kg_hm3",
        hm3_grch38 = lambda wc: "resources/hapmap3_sites_grch38.tsv" if wc.sim_method == "hapnest_public" else [],
    output:
        betas = "results/pgs_weights/prscs/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/betas.tsv",
    params:
        ref_dir        = "results/prscs_custom_ref/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ref_1kg",
        bim_prefix     = "results/plink/{sim_method}/gwas/rep{rep}/merged",
        n_train        = _n_train,
        seed           = lambda wc: BASE_SEED + int(wc.rep),
        hm3_grch38_arg = lambda wc, input: (
            f"--hm3-grch38 {input.hm3_grch38}"
            if wc.sim_method == "hapnest_public" else ""
        ),
    log:
        "logs/prscs/{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    resources:
        mem_mb = 4000,
    shell:
        """
        module load python/3.12
        python3 {input.script} \
            --sumstats   {input.sumstats} \
            --ref-dir    {params.ref_dir} \
            --bim-prefix {params.bim_prefix} \
            --n-gwas     {params.n_train} \
            --seed       {params.seed} \
            --prscs-path {input.prscs_exe} \
            --out        {output.betas} \
            {params.hm3_grch38_arg} \
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
        module load plink/2.0-alpha
        mkdir -p $(dirname {output.sscore})
        TMP_BIM=$(mktemp --suffix=.bim)
        awk 'BEGIN{{OFS="\t"}} $2=="." {{chrom=$1; sub(/^chr/,"",chrom); $2=chrom":"$4":"$6":"$5}} {{print}}' \
            {params.bed_prefix}.bim > "$TMP_BIM"
        plink2 \
            --bed     {params.bed_prefix}.bed \
            --bim     "$TMP_BIM" \
            --fam     {params.bed_prefix}.fam \
            --keep    {input.test_ids} \
            --score   {input.betas} 1 2 3 header \
            --threads 1 \
            --out     {params.score_prefix} \
            2> {log}
        rm -f "$TMP_BIM"
        """


# ── Step 10a: Per-individual posterior PGS variance ───────────────────────────

rule score_pgs_variance:
    """Compute per-individual posterior PGS variance: Var(PGS_i) = w_i^T R w_i, w_i = sqrt(BETA_VAR) * G_i."""
    input:
        bed      = "results/plink/{sim_method}/gwas/rep{rep}/merged.bed",
        betas    = "results/pgs_weights/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/betas.tsv",
        test_ids = "results/splits/{sim_method}/rep{rep}/test.txt",
        ld_sfbm  = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/ld_sfbm.rds",
        ld_snps  = "results/ldpred2_work/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/matched_snps.tsv",
        script   = "scripts/score_pgs_variance.R",
    output:
        var_scores = "results/pgs/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/scores_var.tsv",
    log:
        "logs/score_var/{method}_{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    resources:
        mem_mb = 16000,
    shell:
        """
        module load R/4.4.3-gcc-11.2.0-mkl
        mkdir -p $(dirname {output.var_scores})
        Rscript {input.script} \
            --bed      {input.bed} \
            --betas    {input.betas} \
            --test-ids {input.test_ids} \
            --ld-sfbm  {input.ld_sfbm} \
            --ld-snps  {input.ld_snps} \
            --out      {output.var_scores} \
            2> {log}
        """


# ── Step 10b: Evaluate prediction accuracy ────────────────────────────────────

rule evaluate_pgs:
    """Compute R² (quantitative) or AUC (binary) for the test-set PGS."""
    input:
        scores     = "results/pgs/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/scores.sscore",
        var_scores = "results/pgs/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/scores_var.tsv",
        pheno      = "results/phenotypes/{sim_method}/rep{rep}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/pheno.pheno",
        test_ids   = "results/splits/{sim_method}/rep{rep}/test.txt",
        pcs        = "results/pca/{sim_method}/test/rep{rep}/pcs.sscore",
        script     = "scripts/evaluate_pgs.R",
    output:
        metrics = "results/evaluation/{method}/{sim_method}/rep{rep}/{panel_ancestry}/n{panel_n}/{trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/metrics.tsv",
    params:
        trait_type = lambda wc: trait_type(wc.trait),
    log:
        "logs/evaluation/{method}_{sim_method}_rep{rep}_{panel_ancestry}_n{panel_n}_{trait}_h2_{h2}_pc_{p_causal}_{effect_dist}.log",
    shell:
        """
        module load R/4.4.3-gcc-11.2.0-mkl
        Rscript {input.script} \
            --scores     {input.scores} \
            --var-scores {input.var_scores} \
            --pheno      {input.pheno} \
            --test-ids   {input.test_ids} \
            --covariates {input.pcs} \
            --trait-type {params.trait_type} \
            --out        {output.metrics} \
            2> {log}
        """
