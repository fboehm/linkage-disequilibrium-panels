#!/usr/bin/env python3
"""
Filter HAPNEST public dataset (S-BSST936) to AMR individuals and partition into cohorts.

Reads the PLINK FAM file (sample IDs, in the same row order as the dataset) and
the companion population-label file, keeps only AMR subjects, randomly draws
n_gwas + n_target + n_panel of them, then splits into cohorts.

For each cohort writes:
  keep_{cohort}.txt   — FID<tab>IID per line, for plink2 --keep
  rename_{cohort}.txt — one new sample name per line ({cohort}_0, {cohort}_1, ...),
                        for bcftools reheader -s

Usage (via Snakemake partition_hapnest_public rule):
  python3 partition_hapnest_public.py \\
      --fam      resources/hapnest_public/raw/synthetic_v1_chr-1.fam \\
      --sample   resources/hapnest_public/synthetic_v1.sample \\
      --n-gwas   15000 \\
      --n-target 5000  \\
      --n-panel  10000 \\
      --seed     1000  \\
      --outdir   resources/hapnest_public
"""

import argparse
import random
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fam",      required=True,
                   help="PLINK FAM file (any chromosome; all share the same individuals)")
    p.add_argument("--sample",   required=True,
                   help="Population-label file — one label per line, same order as FAM")
    p.add_argument("--n-gwas",   type=int, required=True)
    p.add_argument("--n-target", type=int, required=True)
    p.add_argument("--n-panel",  type=int, required=True)
    p.add_argument("--seed",     type=int, required=True,
                   help="Random seed for reproducible sampling")
    p.add_argument("--outdir",   required=True,
                   help="Directory where keep_*.txt and rename_*.txt are written")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read FAM: FID IID father mother sex pheno (whitespace-separated)
    samples = []
    with open(args.fam) as fh:
        for line in fh:
            parts = line.split()
            samples.append((parts[0], parts[1]))

    # Read population labels — one per line, same order as FAM rows
    labels = []
    with open(args.sample) as fh:
        for line in fh:
            label = line.strip()
            if label:
                labels.append(label)

    if len(samples) != len(labels):
        raise ValueError(
            f"FAM has {len(samples)} rows but sample file has {len(labels)} labels"
        )

    # Filter to AMR
    amr = [(fid, iid) for (fid, iid), lbl in zip(samples, labels) if lbl == "AMR"]
    n_total = args.n_gwas + args.n_target + args.n_panel
    if len(amr) < n_total:
        raise ValueError(
            f"Only {len(amr)} AMR individuals available; {n_total} requested"
        )

    print(f"  AMR individuals available: {len(amr)}", file=sys.stderr)

    # Randomly select n_total without replacement
    rng      = random.Random(args.seed)
    selected = rng.sample(amr, n_total)

    cohorts = {
        "gwas":   selected[: args.n_gwas],
        "target": selected[args.n_gwas : args.n_gwas + args.n_target],
        "panel":  selected[args.n_gwas + args.n_target :],
    }

    for cohort_name, cohort_samples in cohorts.items():
        keep_path = outdir / f"keep_{cohort_name}.txt"
        with open(keep_path, "w") as fh:
            for fid, iid in cohort_samples:
                fh.write(f"{fid}\t{iid}\n")

        rename_path = outdir / f"rename_{cohort_name}.txt"
        with open(rename_path, "w") as fh:
            for i in range(len(cohort_samples)):
                fh.write(f"{cohort_name}_{i}\n")

        print(f"  {cohort_name}: {len(cohort_samples)} individuals → "
              f"{keep_path.name}, {rename_path.name}", file=sys.stderr)


if __name__ == "__main__":
    main()
