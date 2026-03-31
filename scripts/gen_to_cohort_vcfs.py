#!/usr/bin/env python3
"""
Convert a HAPGEN2 .controls.gen output file to three per-cohort bgzipped VCF
files (gwas, target, panel).

Individuals are split by their sequential position in the .gen file:
  indices   0 .. n_gwas-1                      → gwas_{chrom}.vcf.gz
  indices   n_gwas .. n_gwas+n_target-1        → target_{chrom}.vcf.gz
  indices   n_gwas+n_target .. n_total-1       → panel_{chrom}.vcf.gz

Genotypes are hard-called as the argmax of the three probability columns
(P(0/0), P(0/1), P(1/1)) output by HAPGEN2.

HAPGEN2 .gen format (IMPUTE2 style):
  CHR  SNP_ID  POS  A0  A1  [P(0/0) P(0/1) P(1/1)]*N_individuals

HAPGEN2 .sample format:
  Row 0 (header):  ID_1  ID_2  missing  sex
  Row 1 (types):   0     0     0        D
  Row 2+:          one row per individual (controls)
"""

import argparse
import gzip
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gen",      required=True,
                   help="HAPGEN2 .controls.gen (plain or .gz)")
    p.add_argument("--sample",   required=True,
                   help="HAPGEN2 .controls.sample")
    p.add_argument("--chrom",    required=True,
                   help="Chromosome label for VCF CHROM column (e.g. chr1)")
    p.add_argument("--n-gwas",   type=int, required=True,
                   help="Number of GWAS individuals (first cohort)")
    p.add_argument("--n-target", type=int, required=True,
                   help="Number of target individuals (second cohort)")
    p.add_argument("--n-panel",  type=int, required=True,
                   help="Number of panel individuals (third cohort)")
    p.add_argument("--outdir",   required=True,
                   help="Directory for output VCF.gz files")
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────

def read_sample_ids(sample_file):
    """Parse individual IDs from HAPGEN2 .sample file.

    Skips the 2-row preamble (header + type row); returns a list of ID_1 values.
    """
    ids = []
    with open(sample_file) as fh:
        for i, line in enumerate(fh):
            if i < 2:
                continue
            parts = line.split()
            if parts:
                ids.append(parts[0])
    return ids


def call_gt(p_aa, p_ab, p_bb):
    """Hard-call a diploid genotype from three probability strings."""
    probs = (float(p_aa), float(p_ab), float(p_bb))
    best = probs.index(max(probs))
    return ("0/0", "0/1", "1/1")[best]


def vcf_header(chrom, sample_ids):
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={chrom}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(sample_ids),
    ]
    return "\n".join(lines) + "\n"


def open_gen(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    outdir   = Path(args.outdir)
    n_gwas   = args.n_gwas
    n_target = args.n_target
    n_panel  = args.n_panel
    n_total  = n_gwas + n_target + n_panel
    chrom    = args.chrom

    # Read sample IDs from HAPGEN2 output
    sample_ids = read_sample_ids(args.sample)
    if len(sample_ids) != n_total:
        print(
            f"ERROR: .sample has {len(sample_ids)} individuals, "
            f"expected {n_total} ({n_gwas}+{n_target}+{n_panel}).",
            file=sys.stderr,
        )
        sys.exit(1)

    gwas_ids   = sample_ids[:n_gwas]
    target_ids = sample_ids[n_gwas : n_gwas + n_target]
    panel_ids  = sample_ids[n_gwas + n_target :]

    # Open three output VCF streams simultaneously
    fh_gwas   = gzip.open(str(outdir / f"gwas_{chrom}.vcf.gz"),   "wt")
    fh_target = gzip.open(str(outdir / f"target_{chrom}.vcf.gz"), "wt")
    fh_panel  = gzip.open(str(outdir / f"panel_{chrom}.vcf.gz"),  "wt")

    fh_gwas.write(vcf_header(chrom, gwas_ids))
    fh_target.write(vcf_header(chrom, target_ids))
    fh_panel.write(vcf_header(chrom, panel_ids))

    # Detect the number of leading columns in the .gen file from the first line.
    # HAPGEN2 typically writes 5 (CHR SNP_ID POS A0 A1); some builds write 4.
    n_meta = None
    n_sites = 0

    with open_gen(args.gen) as gen_fh:
        for line in gen_fh:
            parts = line.rstrip().split()

            # Detect leading column count on first line
            if n_meta is None:
                total_cols = len(parts)
                rem = total_cols - 5
                if rem == 3 * n_total:
                    n_meta = 5
                elif rem - 3 == 3 * n_total:   # 4-column prefix
                    n_meta = 4
                else:
                    print(
                        f"ERROR: cannot determine .gen format. "
                        f"Line has {total_cols} columns; expected "
                        f"{5 + 3*n_total} (5-col prefix) or "
                        f"{4 + 3*n_total} (4-col prefix).",
                        file=sys.stderr,
                    )
                    sys.exit(1)

            if len(parts) < n_meta + 3 * n_total:
                print(
                    f"WARNING: skipping short line at site {n_sites + 1}",
                    file=sys.stderr,
                )
                continue

            snp_id = parts[1] if n_meta == 5 else parts[0]
            pos    = parts[2] if n_meta == 5 else parts[1]
            ref    = parts[3] if n_meta == 5 else parts[2]
            alt    = parts[4] if n_meta == 5 else parts[3]
            gcols  = parts[n_meta:]

            # Hard-call genotypes for all individuals
            gts = [
                call_gt(gcols[3*i], gcols[3*i+1], gcols[3*i+2])
                for i in range(n_total)
            ]

            fixed = f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t"

            fh_gwas.write(  fixed + "\t".join(gts[:n_gwas])                        + "\n")
            fh_target.write(fixed + "\t".join(gts[n_gwas:n_gwas+n_target])         + "\n")
            fh_panel.write( fixed + "\t".join(gts[n_gwas+n_target:])               + "\n")

            n_sites += 1

    fh_gwas.close()
    fh_target.close()
    fh_panel.close()

    for name, ids in [("gwas", gwas_ids), ("target", target_ids), ("panel", panel_ids)]:
        print(
            f"  Wrote {name}_{chrom}.vcf.gz  "
            f"({len(ids)} individuals, {n_sites} sites)",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
