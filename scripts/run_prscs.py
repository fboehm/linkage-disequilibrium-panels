#!/usr/bin/env python3
"""
Run PRS-CS to compute per-SNP posterior effect sizes.

Wraps the PRScs.py script (must be on PATH or set via --prscs-path) and
converts its output to the same SNP  A1  BETA  CHR  BP tab-delimited format
used by run_ldpred2.R, so that score_test_set (plink2 --score) can be
applied unchanged.

Inputs
------
--sumstats   PLINK2 --glm output (linear or logistic.hybrid), TEST==ADD rows
--ref-dir    Directory with PRS-CS LD block files (ldblk_1kg_*.hdf5 + snpinfo)
--n-gwas     GWAS training-set sample size
--seed       Random seed
--out        Output TSV: SNP  A1  BETA  CHR  BP

Dependencies
------------
PRScs (https://github.com/getian107/PRScs) must be installed.
Python packages: pandas, h5py (required by PRScs).
"""

import argparse
import os
import subprocess
import sys
import tempfile

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sumstats",   required=True,
                   help="PLINK2 GLM summary statistics file")
    p.add_argument("--ref-dir",    required=True,
                   help="PRS-CS LD block reference directory")
    p.add_argument("--n-gwas",     type=int, required=True,
                   help="GWAS training-set sample size")
    p.add_argument("--seed",       type=int, default=1,
                   help="Random seed [default: 1]")
    p.add_argument("--out",        required=True,
                   help="Output TSV: SNP  A1  BETA  CHR  BP")
    p.add_argument("--prscs-path", default="PRScs.py",
                   help="Path to PRScs.py [default: PRScs.py on PATH]")
    p.add_argument("--bim-prefix", required=True,
                   help="Path prefix of PLINK .bim file for the target dataset "
                        "(e.g. results/plink/hapnest_public/gwas/rep1/merged)")
    p.add_argument("--phi",        default=None,
                   help="Global shrinkage parameter phi (None = auto)")
    return p.parse_args()


def read_sumstats(path: str) -> pd.DataFrame:
    """Read PLINK2 --glm output and normalise to BETA/SE columns."""
    ss = pd.read_csv(path, sep="\t", comment=None, low_memory=False)
    ss.rename(columns={ss.columns[0]: "CHROM"}, inplace=True)
    ss = ss[ss["TEST"] == "ADD"].copy()

    if "BETA" in ss.columns:
        pass
    elif "OR" in ss.columns:
        import math
        ss["BETA"] = ss["OR"].apply(math.log)
        ss.rename(columns={"LOG(OR)_SE": "SE"}, inplace=True)
    else:
        raise ValueError("sumstats file has neither BETA nor OR column")

    ss = ss.dropna(subset=["BETA", "SE"])
    ss = ss[ss["SE"] > 0]
    return ss


def write_prscs_input(ss: pd.DataFrame, path: str) -> None:
    """Write PRS-CS sst file: SNP A1 A2 BETA SE."""
    out = ss[["ID", "A1", "REF", "BETA", "SE"]].copy()
    out.columns = ["SNP", "A1", "A2", "BETA", "SE"]
    out.to_csv(path, sep="\t", index=False)


def run_prscs(prscs_path, ref_dir, bim_prefix, sst_file, n_gwas, seed, out_prefix, phi,
              chroms):
    # PRScs treats --out_dir as a full path prefix, not a directory.
    # Output files are named {out_dir}_pst_eff_a1_b0.5_phi{phi}_chr{chrom}.txt
    chrom_str = ",".join(str(c).lstrip("chr") for c in sorted(chroms))
    cmd = [
        "python3", prscs_path,
        "--ref_dir",    ref_dir,
        "--bim_prefix", bim_prefix,
        "--sst_file",   sst_file,
        "--n_gwas",     str(n_gwas),
        "--seed",       str(seed),
        "--out_dir",    out_prefix,
        "--chrom",      chrom_str,
    ]
    if phi is not None:
        cmd += ["--phi", str(phi)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"PRScs.py exited with code {result.returncode}")
    print(result.stdout, file=sys.stderr)


def remap_to_rsids(ss: pd.DataFrame, ref_dir: str, bim_prefix: str,
                   tmpdir: str) -> tuple:
    """
    Remap positional SNP IDs to rsIDs using the PRS-CS snpinfo file.
    Matches on (CHR, BP).  Both the sumstats DataFrame and a temporary copy of
    the .bim are updated.  Returns (updated_ss, new_bim_prefix).
    """
    snpinfo_path = os.path.join(ref_dir, "snpinfo_1kg_hm3")
    snpinfo = pd.read_csv(snpinfo_path, sep="\t", low_memory=False)
    # columns: CHR  SNP  BP  A1  A2  MAF
    snpinfo["CHR"] = snpinfo["CHR"].astype(str).str.lstrip("chr")
    snpinfo["BP"]  = snpinfo["BP"].astype(int)
    lookup = dict(zip(zip(snpinfo["CHR"], snpinfo["BP"]), snpinfo["SNP"]))

    # ── remap sumstats ────────────────────────────────────────────────────────
    ss = ss.copy()
    chr_col = ss["CHROM"].astype(str).str.lstrip("chr")
    bp_col  = ss["POS"].astype(int)
    ss["_rsID"] = [lookup.get((c, b)) for c, b in zip(chr_col, bp_col)]
    n_mapped = ss["_rsID"].notna().sum()
    print(f"[run_prscs] {n_mapped}/{len(ss)} sumstats SNPs mapped to rsIDs",
          file=sys.stderr)
    ss = ss[ss["_rsID"].notna()].copy()
    ss["ID"] = ss["_rsID"]
    ss.drop(columns=["_rsID"], inplace=True)

    # ── remap bim ─────────────────────────────────────────────────────────────
    bim = pd.read_csv(f"{bim_prefix}.bim", sep="\t", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    chr_col = bim["CHR"].astype(str).str.lstrip("chr")
    bp_col  = bim["BP"].astype(int)
    bim["_rsID"] = [lookup.get((c, b)) for c, b in zip(chr_col, bp_col)]
    bim = bim[bim["_rsID"].notna()].copy()
    bim["SNP"] = bim["_rsID"]
    bim.drop(columns=["_rsID"], inplace=True)

    new_bim_prefix = os.path.join(tmpdir, "rsid_bim")
    bim.to_csv(f"{new_bim_prefix}.bim", sep="\t", header=False, index=False)
    print(f"[run_prscs] {len(bim)} bim SNPs mapped to rsIDs", file=sys.stderr)

    return ss, new_bim_prefix


def collect_prscs_output(out_prefix: str, ss: pd.DataFrame) -> pd.DataFrame:
    """Concatenate per-chromosome PRScs output files into a single data frame."""
    chroms = [str(c).lstrip("chr") for c in ss["CHROM"].unique()]
    parts = []
    for chrom in sorted(chroms, key=int):
        fname = f"{out_prefix}_pst_eff_a1_b0.5_phiauto_chr{chrom}.txt"
        if not os.path.exists(fname):
            print(f"WARNING: {fname} not found, skipping chr{chrom}",
                  file=sys.stderr)
            continue
        df = pd.read_csv(fname, sep="\t", header=None,
                         names=["CHR", "SNP", "BP", "A1", "A2", "BETA"])
        parts.append(df)
    if not parts:
        raise RuntimeError("No PRScs output files were produced.")
    return pd.concat(parts, ignore_index=True)


def main():
    args = parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    ss = read_sumstats(args.sumstats)
    print(f"[run_prscs] {len(ss)} SNPs in sumstats after filtering",
          file=sys.stderr)

    with tempfile.TemporaryDirectory() as tmpdir:
        ss, bim_prefix = remap_to_rsids(ss, args.ref_dir, args.bim_prefix, tmpdir)
        print(f"[run_prscs] {len(ss)} SNPs remaining after rsID mapping",
              file=sys.stderr)

        sst_file   = os.path.join(tmpdir, "sumstats_prscs.txt")
        out_prefix = os.path.join(tmpdir, "prscs_out")

        write_prscs_input(ss, sst_file)
        chroms = ss["CHROM"].unique()
        run_prscs(args.prscs_path, args.ref_dir, bim_prefix, sst_file,
                  args.n_gwas, args.seed, out_prefix, args.phi, chroms)

        results = collect_prscs_output(out_prefix, ss)

    out_df = results[["SNP", "A1", "BETA", "CHR", "BP"]]
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[run_prscs] Wrote {len(out_df)} SNP weights to {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
