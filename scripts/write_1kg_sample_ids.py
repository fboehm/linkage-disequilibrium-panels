#!/usr/bin/env python3
"""
Write a single-column sample-ID list for one 1KG super-population.

The output file contains one sample ID per line (suitable for
bcftools --samples-file).

Usage:
    python3 scripts/write_1kg_sample_ids.py \
        --panel resources/1kg/integrated_call_samples.panel \
        --superpop EUR \
        --out resources/panels/EUR_1kg/sample_ids.txt
"""

import argparse
import csv
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel",    required=True,
                   help="Path to integrated_call_samples.panel")
    p.add_argument("--superpop", required=True,
                   help="Super-population code (EUR, AFR, AMR, EAS, SAS)")
    p.add_argument("--out",      required=True,
                   help="Output file: one sample ID per line")
    return p.parse_args()


def main():
    args = parse_args()

    samples = []
    with open(args.panel) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["super_pop"].strip() == args.superpop.strip():
                samples.append(row["sample"].strip())

    if not samples:
        print(f"ERROR: no samples found for super_pop={args.superpop}",
              file=sys.stderr)
        sys.exit(1)

    with open(args.out, "w") as fh:
        for s in samples:
            fh.write(s + "\n")

    print(f"Wrote {len(samples)} sample IDs for {args.superpop} → {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
