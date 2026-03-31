#!/usr/bin/env python3
"""
Simulate genotypes using HAPNEST for one chromosome and replicate.

Runs HAPNEST via a Singularity container to generate n_gwas + n_target +
n_panel synthetic individuals, then splits the output into three per-cohort
bgzipped VCF files:

  gwas_{chrom}.vcf.gz    -- GWAS discovery cohort
  target_{chrom}.vcf.gz  -- Independent evaluation cohort
  panel_{chrom}.vcf.gz   -- Matched LD reference panel cohort

Sample IDs are renamed to {cohort}_{0-based index} to match the naming
convention used by simulate_genotypes.py (msprime).

Requirements (on PATH inside container): plink2, bcftools, tabix
Requirements (on host): singularity

Usage (via Snakemake rule):
  python3 simulate_hapnest.py \\
      --chrom          chr21 \\
      --n-gwas         15000 \\
      --n-target       5000 \\
      --n-panel        10000 \\
      --seed           1001 \\
      --superpopulation AMR \\
      --outdir         results/vcf/hapnest/rep1 \\
      --hapnest-container resources/hapnest/intervene-synthetic-data_latest.sif \\
      --hapnest-data-dir  resources/hapnest/data \\
      --threads        4
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path

import yaml


# ── Default HAPNEST population-genetic parameters (posterior estimates from
#    the HAPNEST publication, Wharrie et al. 2023, Bioinformatics) ─────────────

DEFAULT_RHO = {"AFR": 0.77, "AMR": 0.80, "EAS": 0.58,
               "EUR": 0.68, "CSA": 0.73, "MID": 0.65}
DEFAULT_NE  = {"AFR": 11900, "AMR": 10400, "EAS": 11700,
               "EUR": 11700, "CSA": 11500, "MID": 8100}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chrom", required=True,
                   help="Chromosome label, e.g. chr21")
    p.add_argument("--n-gwas", type=int, required=True,
                   help="Number of diploid GWAS individuals")
    p.add_argument("--n-target", type=int, required=True,
                   help="Number of diploid target individuals")
    p.add_argument("--n-panel", type=int, required=True,
                   help="Number of diploid LD-panel individuals")
    p.add_argument("--seed", type=int, required=True,
                   help="Random seed for HAPNEST")
    p.add_argument("--superpopulation", default="AMR",
                   choices=["AFR", "AMR", "EAS", "EUR", "CSA", "MID"],
                   help="1KG+HGDP superpopulation to use as reference "
                        "(default: AMR)")
    p.add_argument("--outdir", required=True,
                   help="Directory where output VCF files are written")
    p.add_argument("--hapnest-container", required=True,
                   help="Singularity: path to .sif file.  "
                        "Docker: image name (e.g. sophiewharrie/intervene-synthetic-data)")
    p.add_argument("--hapnest-data-dir", required=True,
                   help="Path to the HAPNEST data directory — must contain "
                        "inputs/ populated by the container 'fetch' command")
    p.add_argument("--use-docker", action="store_true",
                   help="Use Docker instead of Singularity to run HAPNEST")
    p.add_argument("--threads", type=int, default=4,
                   help="CPU threads passed to HAPNEST generate_geno "
                        "(default: 4)")
    p.add_argument("--memory-mb", type=int, default=30000,
                   help="Memory cap passed to HAPNEST in MB (default: 30000)")
    return p.parse_args()


def chrom_num(chrom: str) -> str:
    """Strip 'chr' prefix: 'chr21' → '21'."""
    return str(chrom).lstrip("chr")


def write_hapnest_config(args, config_path: Path, hapnest_outdir: str) -> None:
    """
    Write a HAPNEST config YAML for this chromosome / replicate.

    All file paths use the /data/ prefix because the host data directory is
    bind-mounted at /data/ inside the Singularity container.  HAPNEST replaces
    the literal string {chromosome} with the chromosome number at runtime.
    """
    chrom_n   = chrom_num(args.chrom)
    n_total   = args.n_gwas + args.n_target + args.n_panel
    pop       = args.superpopulation
    data_root = "/data"          # bind-mount point inside the container

    # File-path patterns — HAPNEST substitutes {chromosome} internally.
    inp_raw  = f"{data_root}/inputs/raw/1KG+HGDP"
    inp_proc = f"{data_root}/inputs/processed/1KG+HGDP"

    config = {
        "global_parameters": {
            "random_seed":    args.seed,
            "chromosome":     int(chrom_n),
            "superpopulation": pop,
            "memory":         args.memory_mb,
            "batchsize":      10000,
        },
        "filepaths": {
            "general": {
                "output_dir":    hapnest_outdir,
                "output_prefix": f"hapnest_chr-{{chromosome}}",
            },
            "genotype": {
                "vcf_input_raw":       f"{inp_raw}/1KG+HGDP.chr{{chromosome}}.hapmap.final.vcf.gz",
                "vcf_input_processed": f"{inp_proc}/1KG+HGDP.chr{{chromosome}}.hapmap.final.recode.vcf",
                "popfile_raw":         f"{inp_proc}/merged_pop_adjusted.tsv",
                "popfile_processed":   f"{inp_proc}/merged_pop.tsv",
                "genetic_mapfile":     f"{inp_raw}/genetic_maps/chr{{chromosome}}.interpolated_genetic_map",
                "genetic_distfile":    f"{inp_proc}/1KG+HGDP.chr{{chromosome}}.hapmap.distfile",
                "variant_list":        f"{inp_proc}/hapmap_variant_list_chr{{chromosome}}.txt",
                "rsid_list":           f"{inp_proc}/rsid_map_list_chr{{chromosome}}.txt",
                "vcf_metadata":        f"{inp_proc}/1KG+HGDP.chr{{chromosome}}.metadata",
                "mutation_mapfile":    f"{inp_raw}/mutation_maps/atlas.chr{{chromosome}}.csv",
                "mutation_agefile":    f"{inp_proc}/1KG+HGDP.chr{{chromosome}}.hapmap.agefile",
                "hap1_matrix":         f"{inp_proc}/1KG+HGDP.chr{{chromosome}}.hapmap.h1",
                "hap2_matrix":         f"{inp_proc}/1KG+HGDP.chr{{chromosome}}.hapmap.h2",
                "remove_list":         f"{inp_proc}/remove.txt",
            },
            "phenotype": {
                "causal_list": "",
                "reference":   "",
            },
            "software": {
                "plink":    "plink",
                "plink2":   "plink2",
                "king":     "king",
                "vcftools": "vcftools",
                "mapthin":  "mapthin",
                "phenoalg": "phenoalg",
            },
        },
        "genotype_data": {
            "samples": {
                "use_default": True,
                "default": {"nsamples": n_total},
            },
            "rho": DEFAULT_RHO,
            "Ne":  DEFAULT_NE,
        },
        # Phenotype and evaluation sections are required by the config schema
        # even when only genotypes are being generated.
        "phenotype_data": {
            "nPopulation":   1,
            "nTrait":        1,
            "ProportionGeno":  "0.5",
            "ProportionCovar": "0",
            "Prevalence":      "0.5",
            "TraitCorr":       1,
            "PopulationCorr":  [0],
            "Causality": {
                "UseCausalList": False,
                "Polygenicity":  0.005,
                "Pleiotropy":    1,
            },
        },
        "evaluation": {
            "metrics": {
                "aats":     False,
                "kinship":  False,
                "ld_corr":  False,
                "ld_decay": False,
                "maf":      False,
                "pca":      False,
                "gwas":     False,
            }
        },
    }

    with open(config_path, "w") as fh:
        yaml.dump(config, fh, default_flow_style=False, sort_keys=False)

    print(f"  Wrote HAPNEST config: {config_path}", file=sys.stderr)


def run_hapnest(args, config_path: Path) -> None:
    """Run HAPNEST generate_geno via Singularity or Docker."""
    data_dir   = Path(args.hapnest_data_dir).resolve()
    config_abs = config_path.resolve()

    # Config path as seen inside the container (/data/ bind-mount).
    config_in_container = "/data/" + str(config_abs.relative_to(data_dir))

    if args.use_docker:
        cmd = [
            "docker", "run", "--rm",
            "--volume", f"{data_dir}:/data/",
            args.hapnest_container,
            "generate_geno",
            str(args.threads),
            config_in_container,
        ]
    else:
        cmd = [
            "singularity", "exec",
            "--no-home",
            "--env", "JULIA_DEPOT_PATH=/root/.julia",
            "--bind", f"{data_dir}:/data/",
            args.hapnest_container,
            "generate_geno",
            str(args.threads),
            config_in_container,
        ]

    print(f"Running HAPNEST: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)


def split_and_convert(plink_prefix: str, args, outdir: Path) -> None:
    """
    Split the HAPNEST PLINK output into three cohorts and write VCFs.

    HAPNEST generates individuals in the order they were requested, so the
    FAM file rows map directly onto cohorts:
        rows [0,          n_gwas)            → gwas
        rows [n_gwas,     n_gwas+n_target)   → target
        rows [n_gwas+n_target, n_total)      → panel
    """
    fam_path = Path(plink_prefix + ".fam")
    if not fam_path.exists():
        raise FileNotFoundError(
            f"HAPNEST PLINK output not found: {fam_path}\n"
            "Check the HAPNEST log for errors."
        )

    # Read all sample IDs from the FAM file.
    all_samples = []
    with open(fam_path) as fh:
        for line in fh:
            parts = line.split()
            all_samples.append((parts[0], parts[1]))   # (FID, IID)

    n_total = args.n_gwas + args.n_target + args.n_panel
    if len(all_samples) != n_total:
        raise ValueError(
            f"Expected {n_total} samples in FAM file, found {len(all_samples)}."
        )

    cohort_slices = {
        "gwas":   all_samples[:args.n_gwas],
        "target": all_samples[args.n_gwas : args.n_gwas + args.n_target],
        "panel":  all_samples[args.n_gwas + args.n_target :],
    }

    chrom_n = chrom_num(args.chrom)

    for cohort_name, samples in cohort_slices.items():
        # ── 1. Write a keep file for plink2 ──────────────────────────────────
        keep_file = outdir / f"_tmp_{cohort_name}.keep"
        with open(keep_file, "w") as fh:
            for fid, iid in samples:
                fh.write(f"{fid}\t{iid}\n")

        # ── 2. Extract cohort → uncompressed VCF ─────────────────────────────
        tmp_prefix = outdir / f"_tmp_{cohort_name}"
        subprocess.run(
            [
                "plink2",
                "--bfile",  plink_prefix,
                "--keep",   str(keep_file),
                "--export", "vcf-4.2", "bgz",
                "--out",    str(tmp_prefix),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )

        # ── 3. Rename samples to {cohort}_{i} format ─────────────────────────
        #    bcftools reheader -s expects one new name per line, in FAM order.
        rename_file = outdir / f"_tmp_{cohort_name}.rename"
        with open(rename_file, "w") as fh:
            for i in range(len(samples)):
                fh.write(f"{cohort_name}_{i}\n")

        vcf_out = outdir / f"{cohort_name}_{args.chrom}.vcf.gz"
        subprocess.run(
            [
                "bcftools", "reheader",
                "--samples", str(rename_file),
                "--output",  str(vcf_out),
                str(tmp_prefix) + ".vcf.gz",
            ],
            check=True,
        )

        # ── 4. Index ──────────────────────────────────────────────────────────
        subprocess.run(["tabix", "-p", "vcf", str(vcf_out)], check=True)

        print(
            f"  Wrote {vcf_out}  ({len(samples)} individuals)",
            file=sys.stderr,
        )

        # ── 5. Clean up temporaries ───────────────────────────────────────────
        for tmp in [
            keep_file,
            rename_file,
            Path(str(tmp_prefix) + ".vcf.gz"),
            Path(str(tmp_prefix) + ".vcf.gz.tbi"),
            Path(str(tmp_prefix) + ".log"),
            Path(str(tmp_prefix) + ".pgen"),
            Path(str(tmp_prefix) + ".psam"),
            Path(str(tmp_prefix) + ".pvar"),
        ]:
            if tmp.exists():
                tmp.unlink()


def main():
    args   = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    chrom_n   = chrom_num(args.chrom)
    data_dir  = Path(args.hapnest_data_dir).resolve()

    # Config is written into the data dir so it's accessible inside the container.
    config_path = data_dir / f"config_rep{args.seed}_chr{chrom_n}.yaml"

    # HAPNEST output directory (as seen inside the container).
    hapnest_outdir_container = f"/data/outputs/rep{args.seed}/chr{chrom_n}"
    # … and on the host filesystem:
    hapnest_outdir_host      = data_dir / "outputs" / f"rep{args.seed}" / f"chr{chrom_n}"
    hapnest_outdir_host.mkdir(parents=True, exist_ok=True)

    print(
        f"HAPNEST simulation: {args.chrom}, {args.n_gwas + args.n_target + args.n_panel} "
        f"individuals ({args.superpopulation}), seed={args.seed}",
        file=sys.stderr,
    )

    # ── Step 1: Write config ──────────────────────────────────────────────────
    write_hapnest_config(args, config_path, hapnest_outdir_container)

    # ── Step 2: Run HAPNEST ───────────────────────────────────────────────────
    run_hapnest(args, config_path)

    # ── Step 3: Split output into cohorts and convert to per-cohort VCFs ─────
    plink_prefix = str(hapnest_outdir_host / f"hapnest_chr-{chrom_n}")
    split_and_convert(plink_prefix, args, outdir)

    # ── Step 4: Remove per-run config ─────────────────────────────────────────
    config_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
