#!/usr/bin/env python3
"""
Filter 1000 Genomes Phase 3 Oxford-format haplotypes to a target
super-population and to HapMap3 SNP positions.

Inputs
------
--haps    : 1KG .haps.gz  (IMPUTE2/SHAPEIT2 format; 5 leading metadata columns
            then haplotype values: SNP_ID RSID POS A0 A1 hap1_ind1 hap2_ind1 ...)
--legend  : 1KG .legend.gz (header row: id position a0 a1 [type])
--sample  : Global 1KG .sample file (all populations).
            Expected format:
              Row 0 (header): ID_1 ID_2 missing sex population group
              Row 1 (types):  0    0    0       D   NA         NA
              Row 2+:         one row per individual
            where 'group' is the super-population code (AMR, EUR, AFR, SAS, EAS).
--hapmap3 : HapMap3 sites TSV with header; required columns: chrom  pos
            ('chrom' may or may not carry a 'chr' prefix)
--pop     : Super-population code to retain (AMR, EUR, AFR, SAS, or EAS)
--chrom   : Chromosome number WITHOUT 'chr' prefix (e.g. 1, 21)
--out-haps    : Output filtered .haps.gz
--out-legend  : Output filtered .legend.gz

The output .haps.gz retains the 5 metadata columns followed by haplotype columns
for the target population only, and only rows at HapMap3 positions.
The output .legend.gz retains the header and rows at HapMap3 positions.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

# Number of leading metadata columns in a 1KG Phase 3 .haps file
# (SNP_ID  RSID  POS  A0  A1)
HAPS_META_COLS = 5


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--haps",       required=True,
                   help="1KG per-chromosome .haps.gz")
    p.add_argument("--legend",     required=True,
                   help="1KG per-chromosome .legend.gz")
    p.add_argument("--sample",     required=True,
                   help="Global 1KG .sample file (all populations)")
    p.add_argument("--hapmap3",    required=True,
                   help="HapMap3 sites TSV (columns: chrom pos ...)")
    p.add_argument("--pop",        required=True,
                   help="Super-population code: AMR, EUR, AFR, SAS, or EAS")
    p.add_argument("--chrom",      required=True,
                   help="Chromosome number without 'chr' prefix (e.g. 1, 21)")
    p.add_argument("--out-haps",   required=True,
                   help="Output filtered .haps.gz")
    p.add_argument("--out-legend", required=True,
                   help="Output filtered .legend.gz")
    return p.parse_args()


def load_hapmap3_positions(hapmap3_file, chrom):
    """Return set of HapMap3 bp positions for *chrom* (no 'chr' prefix)."""
    positions = set()
    with open(hapmap3_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["chrom"].lstrip("chr") == chrom:
                positions.add(int(row["pos"]))
    return positions


def get_population_indices(sample_file, target_pop):
    """Return sorted 0-based individual indices for *target_pop* super-population.

    Skips the two-row preamble (header + type row) of the 1KG .sample file.
    Index 0 corresponds to the first data row (row 2 of the file).
    """
    indices = []
    with open(sample_file) as fh:
        for file_row, line in enumerate(fh):
            if file_row < 2:          # skip header and type rows
                continue
            parts = line.split()
            # Expected: ID_1 ID_2 missing sex population group
            if len(parts) < 6:
                print(
                    f"WARNING: sample row {file_row} has only {len(parts)} "
                    f"columns (expected ≥ 6); skipping",
                    file=sys.stderr,
                )
                continue
            group = parts[5]          # super-population column
            if group == target_pop:
                indices.append(file_row - 2)   # 0-based individual index
    return sorted(indices)


def open_maybe_gz(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def main():
    args = parse_args()

    # ── Load HapMap3 positions ────────────────────────────────────────────────
    print(f"Loading HapMap3 positions for chr{args.chrom} ...", file=sys.stderr)
    hm3_pos = load_hapmap3_positions(args.hapmap3, args.chrom)
    print(f"  {len(hm3_pos)} positions", file=sys.stderr)

    # ── Identify target-population individuals ────────────────────────────────
    print(f"Identifying {args.pop} individuals in {args.sample} ...", file=sys.stderr)
    pop_indices = get_population_indices(args.sample, args.pop)
    n_pop = len(pop_indices)
    print(f"  {n_pop} individuals", file=sys.stderr)

    if n_pop == 0:
        print(
            f"ERROR: no individuals found for population '{args.pop}'. "
            f"Check that the sample file has a 'group' column with this code.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Build list of haplotype column indices for the target population.
    # Individual i → columns HAPS_META_COLS + 2*i  and  HAPS_META_COLS + 2*i + 1
    hap_cols = []
    for idx in pop_indices:
        hap_cols.append(HAPS_META_COLS + 2 * idx)
        hap_cols.append(HAPS_META_COLS + 2 * idx + 1)

    # ── Pass 1: scan legend, build row mask and kept legend rows ─────────────
    print("Scanning legend for HapMap3 intersection ...", file=sys.stderr)
    row_mask = []       # True = keep this SNP
    kept_legend = []    # kept legend lines (strings)
    legend_header = None

    with open_maybe_gz(args.legend) as lh:
        for i, line in enumerate(lh):
            if i == 0:
                legend_header = line
                continue
            parts = line.split()
            pos = int(parts[1])     # 'position' column in legend
            keep = pos in hm3_pos
            row_mask.append(keep)
            if keep:
                kept_legend.append(line)

    n_kept = sum(row_mask)
    print(
        f"  {n_kept} / {len(row_mask)} SNPs retained "
        f"(HapMap3 ∩ chr{args.chrom})",
        file=sys.stderr,
    )

    if n_kept == 0:
        print(
            f"ERROR: no SNPs remain after HapMap3 intersection on chr{args.chrom}. "
            f"Check that 'chrom' values in the HapMap3 file match '{args.chrom}'.",
            file=sys.stderr,
        )
        sys.exit(1)

    # ── Write filtered legend ─────────────────────────────────────────────────
    Path(args.out_legend).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out_legend, "wt") as lout:
        lout.write(legend_header)
        for row in kept_legend:
            lout.write(row)

    # ── Pass 2: stream haps, filter rows and columns ──────────────────────────
    print("Filtering haplotype file ...", file=sys.stderr)
    Path(args.out_haps).parent.mkdir(parents=True, exist_ok=True)
    n_written = 0

    with open_maybe_gz(args.haps) as hin, gzip.open(args.out_haps, "wt") as hout:
        for snp_i, line in enumerate(hin):
            if snp_i >= len(row_mask):
                break
            if not row_mask[snp_i]:
                continue
            parts = line.rstrip().split()
            meta = parts[:HAPS_META_COLS]
            haps = [parts[j] for j in hap_cols]
            hout.write(" ".join(meta + haps) + "\n")
            n_written += 1

    print(
        f"Done. Wrote {n_written} SNPs × {n_pop} individuals "
        f"({2 * n_pop} haplotypes).",
        file=sys.stderr,
    )
    print(f"  haps   → {args.out_haps}", file=sys.stderr)
    print(f"  legend → {args.out_legend}", file=sys.stderr)


if __name__ == "__main__":
    main()
