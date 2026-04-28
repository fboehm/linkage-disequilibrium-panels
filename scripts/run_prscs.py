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
    p.add_argument("--hm3-grch38", dest="hm3_grch38", default=None, metavar="FILE",
                   help="TSV with HapMap3 SNP positions on GRCh38: chrom pos rsid a0 a1. "
                        "Used when input data is GRCh38 (e.g. hapnest_public). "
                        "Generate from map_hm3_ldpred2.rds via extract_hm3_grch38_sites rule.")
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
    cmd += ["--write_pst", "TRUE"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"PRScs.py exited with code {result.returncode}")
    print(result.stdout, file=sys.stderr)


def remap_to_rsids(ss: pd.DataFrame, ref_dir: str, bim_prefix: str,
                   tmpdir: str, hm3_grch38=None) -> tuple:
    """
    Remap positional SNP IDs to rsIDs so they match the PRS-CS snpinfo.

    Two modes:
    - GRCh37 (default): match bim/sumstats (CHR, BP) against snpinfo GRCh37 positions.
    - GRCh38: when hm3_grch38 is provided, map GRCh38 (CHR, POS) → rsID using the
      HapMap3 GRCh38 sites table (pos_hg38 from map_hm3_ldpred2.rds), then keep only
      rsIDs present in the snpinfo.

    Returns (updated_ss, new_bim_prefix, rsid_to_orig_id) where rsid_to_orig_id is a
    dict mapping rsID → original BIM SNP ID, used to restore original IDs in output.
    """
    snpinfo_path = os.path.join(ref_dir, "snpinfo_1kg_hm3")
    snpinfo = pd.read_csv(snpinfo_path, sep="\t", low_memory=False)
    # columns: CHR  SNP  BP  A1  A2  MAF

    if hm3_grch38:
        # Build (GRCh38_CHR, GRCh38_POS) → rsID from the HapMap3 GRCh38 sites table,
        # keeping only rsIDs that appear in the snpinfo (HM3 SNPs PRScs knows about).
        snpinfo_rsids = set(snpinfo["SNP"])
        hm3 = pd.read_csv(hm3_grch38, sep="\t")
        hm3 = hm3[hm3["rsid"].isin(snpinfo_rsids)]
        hm3["chrom"] = hm3["chrom"].astype(str)
        hm3["pos"]   = hm3["pos"].astype(int)
        lookup = dict(zip(zip(hm3["chrom"], hm3["pos"]), hm3["rsid"]))
    else:
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
    if n_mapped == 0 and len(ss) > 0:
        sample_ss = list(zip(chr_col[:5], bp_col[:5]))
        sample_lookup = list(lookup.items())[:5]
        print(f"[run_prscs] DEBUG: first 5 sumstats (CHR, POS): {sample_ss}",
              file=sys.stderr)
        print(f"[run_prscs] DEBUG: first 5 lookup keys (CHR, POS): {sample_lookup}",
              file=sys.stderr)
    ss = ss[ss["_rsID"].notna()].copy()
    ss["ID"] = ss["_rsID"]
    ss.drop(columns=["_rsID"], inplace=True)

    # Deduplicate by rsID: multi-allelic sites (two SNPs at same position) get
    # the same rsID; PRScs indexes by position and will crash on duplicates.
    n_before = len(ss)
    ss = ss.drop_duplicates(subset=["ID"])
    n_dropped = n_before - len(ss)
    if n_dropped > 0:
        print(f"[run_prscs] WARNING: dropped {n_dropped} duplicate rsID(s) from multi-allelic sites",
              file=sys.stderr)

    # ── remap bim ─────────────────────────────────────────────────────────────
    bim = pd.read_csv(f"{bim_prefix}.bim", sep="\t", header=None,
                      names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    chr_col = bim["CHR"].astype(str).str.lstrip("chr")
    bp_col  = bim["BP"].astype(int)
    bim["_rsID"] = [lookup.get((c, b)) for c, b in zip(chr_col, bp_col)]
    bim = bim[bim["_rsID"].notna()].copy()
    # Capture original positional IDs before overwriting SNP column
    orig_snp_ids = bim["SNP"].values.copy()
    bim["SNP"] = bim["_rsID"]
    bim["CHR"] = chr_col[bim.index]
    bim.drop(columns=["_rsID"], inplace=True)
    print(f"[run_prscs] {len(bim)} bim SNPs mapped to rsIDs", file=sys.stderr)

    # Deduplicate BIM by rsID: multi-allelic sites produce duplicate rsIDs which
    # cause an IndexError in PRScs parse_ldblk. Keep first occurrence.
    keep_mask = ~bim["SNP"].duplicated()
    n_bim_dropped = (~keep_mask).sum()
    if n_bim_dropped > 0:
        print(f"[run_prscs] WARNING: dropped {n_bim_dropped} duplicate rsID(s) from BIM multi-allelic sites",
              file=sys.stderr)
    orig_snp_ids = orig_snp_ids[keep_mask.values]
    bim = bim[keep_mask]

    rsid_to_orig_id = dict(zip(bim["SNP"], orig_snp_ids))

    new_bim_prefix = os.path.join(tmpdir, "rsid_bim")
    bim.to_csv(f"{new_bim_prefix}.bim", sep="\t", header=False, index=False)

    return ss, new_bim_prefix, rsid_to_orig_id


def collect_prscs_output(out_prefix: str, ss: pd.DataFrame) -> pd.DataFrame:
    """Concatenate per-chromosome PRScs output files into a single data frame.

    With --write_pst TRUE the output has 5 fixed columns (CHR SNP BP A1 A2)
    followed by N posterior-sample columns.  We compute the posterior mean
    (BETA) and posterior variance (BETA_VAR) across those samples.
    """
    chroms = [str(c).lstrip("chr") for c in ss["CHROM"].unique()]
    parts = []
    for chrom in sorted(chroms, key=int):
        fname = f"{out_prefix}_pst_eff_a1_b0.5_phiauto_chr{chrom}.txt"
        if not os.path.exists(fname):
            print(f"WARNING: {fname} not found, skipping chr{chrom}",
                  file=sys.stderr)
            continue
        df = pd.read_csv(fname, sep="\t", header=None)
        n_fixed = 5
        sample_cols = [f"_s{i}" for i in range(df.shape[1] - n_fixed)]
        df.columns = ["CHR", "SNP", "BP", "A1", "A2"] + sample_cols
        if sample_cols:
            df["BETA"]     = df[sample_cols].mean(axis=1)
            df["BETA_VAR"] = df[sample_cols].var(axis=1, ddof=1)
        else:
            df["BETA"]     = df.iloc[:, 5]
            df["BETA_VAR"] = float("nan")
        parts.append(df[["CHR", "SNP", "BP", "A1", "A2", "BETA", "BETA_VAR"]])
    if not parts:
        raise RuntimeError("No PRScs output files were produced.")
    return pd.concat(parts, ignore_index=True)


def main():
    args = parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    ss = read_sumstats(args.sumstats)
    n_snps_in_sumstats = len(ss)
    print(f"[run_prscs] {n_snps_in_sumstats} SNPs in sumstats after filtering",
          file=sys.stderr)

    with tempfile.TemporaryDirectory() as tmpdir:
        ss, bim_prefix, rsid_to_orig_id = remap_to_rsids(ss, args.ref_dir, args.bim_prefix,
                                                         tmpdir, hm3_grch38=args.hm3_grch38)
        n_snps_after_rsid_mapping = len(ss)
        print(f"[run_prscs] {n_snps_after_rsid_mapping} SNPs remaining after rsID mapping",
              file=sys.stderr)
        if len(ss) == 0:
            raise RuntimeError(
                "0 SNPs remaining after rsID mapping. "
                "Check that the rsID map files cover the correct genome build "
                "and chromosomes, and that sumstats CHR/POS match the reference."
            )

        sst_file   = os.path.join(tmpdir, "sumstats_prscs.txt")
        out_prefix = os.path.join(tmpdir, "prscs_out")

        write_prscs_input(ss, sst_file)
        chroms = ss["CHROM"].unique()
        run_prscs(args.prscs_path, args.ref_dir, bim_prefix, sst_file,
                  args.n_gwas, args.seed, out_prefix, args.phi, chroms)

        results = collect_prscs_output(out_prefix, ss)

    results["SNP"] = results["SNP"].map(rsid_to_orig_id).fillna(results["SNP"])
    out_df = results[["SNP", "A1", "BETA", "BETA_VAR", "CHR", "BP"]]
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[run_prscs] Wrote {len(out_df)} SNP weights to {args.out}",
          file=sys.stderr)

    conv_path = os.path.join(os.path.dirname(args.out), "convergence.tsv")
    n_snps_nonzero = int((out_df["BETA"].abs() > 1e-12).sum())
    beta_l2 = float((out_df["BETA"] ** 2).sum() ** 0.5)
    phi_value = float("nan") if args.phi is None else float(args.phi)
    conv_df = pd.DataFrame({
        "metric": ["n_snps_in_sumstats", "n_snps_after_rsid_mapping",
                   "n_snps_in_output", "n_snps_nonzero", "beta_l2",
                   "phi_value", "seed"],
        "value":  [n_snps_in_sumstats, n_snps_after_rsid_mapping,
                   len(out_df), n_snps_nonzero, beta_l2,
                   phi_value, args.seed],
    })
    conv_df.to_csv(conv_path, sep="\t", index=False)
    print(f"[run_prscs] Wrote convergence diagnostics to {conv_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
