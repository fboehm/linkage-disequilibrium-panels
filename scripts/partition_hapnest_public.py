#!/usr/bin/env python3
"""
Partition HAPNEST public dataset (S-BSST936) subjects into cohorts.

AMR individuals:
  Reads the PLINK FAM file and the companion population-label file, keeps only
  AMR subjects, randomly draws n_gwas + n_target + n_panel of them, then splits
  into cohorts.

Non-AMR ancestries (--non-amr-ancestries):
  For each listed ancestry label, writes a keep file containing ALL individuals
  of that ancestry.  Subsampling to a specific panel_n happens later in
  prepare_ldpred2_ref.R / build_prscs_ref_custom.py.

For each cohort / ancestry writes:
  keep_{name}.txt   — FID<tab>IID per line, for plink2 --keep
  rename_{name}.txt — one new sample name per line ({name}_0, {name}_1, ...),
                      for bcftools reheader -s

Usage (via Snakemake partition_hapnest_public rule):
  python3 partition_hapnest_public.py \\
      --fam                resources/hapnest_public/raw/synthetic_v1_chr-1.fam \\
      --sample             resources/hapnest_public/synthetic_v1.sample \\
      --n-gwas             100000 \\
      --n-target           5000   \\
      --n-panel            50000  \\
      --seed               1000   \\
      --outdir             resources/hapnest_public \\
      --non-amr-ancestries AFR EUR EAS CSA MID
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
    p.add_argument("--non-amr-ancestries", nargs="*", default=[],
                   metavar="POP",
                   help="Additional ancestry labels (e.g. AFR EUR EAS CSA MID) for "
                        "which keep/rename files are written using ALL available subjects")
    return p.parse_args()


def write_cohort(outdir: Path, name: str, cohort_samples: list) -> None:
    """Write keep_{name}.txt and rename_{name}.txt for a cohort."""
    keep_path = outdir / f"keep_{name}.txt"
    with open(keep_path, "w") as fh:
        for fid, iid in cohort_samples:
            fh.write(f"{fid}\t{iid}\n")

    rename_path = outdir / f"rename_{name}.txt"
    with open(rename_path, "w") as fh:
        for i in range(len(cohort_samples)):
            fh.write(f"{name}_{i}\n")

    print(f"  {name}: {len(cohort_samples)} individuals → "
          f"{keep_path.name}, {rename_path.name}", file=sys.stderr)


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

    # ── AMR: gwas / target / panel partition ──────────────────────────────────
    amr = [(fid, iid) for (fid, iid), lbl in zip(samples, labels) if lbl == "AMR"]
    n_total = args.n_gwas + args.n_target + args.n_panel
    if len(amr) < n_total:
        raise ValueError(
            f"Only {len(amr)} AMR individuals available; "
            f"{n_total} requested (n_gwas={args.n_gwas}, "
            f"n_target={args.n_target}, n_panel={args.n_panel})"
        )

    print(f"  AMR individuals available: {len(amr)}", file=sys.stderr)

    rng      = random.Random(args.seed)
    selected = rng.sample(amr, n_total)

    write_cohort(outdir, "gwas",   selected[: args.n_gwas])
    write_cohort(outdir, "target", selected[args.n_gwas : args.n_gwas + args.n_target])
    write_cohort(outdir, "panel",  selected[args.n_gwas + args.n_target :])

    # ── Non-AMR ancestries: all subjects as LD-panel candidates ───────────────
    for anc in args.non_amr_ancestries:
        pool = [(fid, iid) for (fid, iid), lbl in zip(samples, labels) if lbl == anc]
        if len(pool) == 0:
            print(f"  WARNING: no individuals found for ancestry '{anc}'", file=sys.stderr)
            continue
        print(f"  {anc} individuals available: {len(pool)}", file=sys.stderr)
        write_cohort(outdir, anc, pool)


if __name__ == "__main__":
    main()
