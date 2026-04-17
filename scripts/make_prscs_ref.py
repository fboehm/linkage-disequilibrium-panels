#!/usr/bin/env python3
"""
Build a PRS-CS compatible HDF5 LD reference from a PLINK BED file.

LD block boundaries are found using the Minimum Description Length (MDL)
criterion from Berisa & Pickrell 2016 (Bioinformatics 32:283-285).

For a block of k SNPs with empirical correlation matrix R_k and n individuals:

    MDL(block) = n/2 * log|R_k| + k*(k-1)/4 * log(n)

A dynamic programming search over all valid partitions minimises the total MDL.
Block boundaries are found via Schur-complement incremental updates — the full
m×m correlation matrix is never stored, only a banded version up to
--max-block-size offsets wide.

For each requested chromosome, writes:
  {out}/ldblk_1kg_chr{N}.hdf5  – LD block matrices in PRS-CS format
  {out}/snpinfo_1kg_hm3         – SNP info table (CHR SNP BP A1 A2 MAF)

The output directory basename must contain '1kg' or 'ukbb' (PRScs.py
uses the basename to decide the HDF5 filename pattern inside parse_ldblk).

Usage (standalone):
  python3 scripts/make_prscs_ref.py \\
      --bfile         resources/panels/AMR_1kg/merged \\
      --out           results/prscs_custom_ref/AMR_1kg \\
      --chroms        22 \\
      --max-block-size 400

Usage (via Snakemake build_prscs_ref_custom rule):
  arguments are filled in from wildcards / params by the rule.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import h5py


# PLINK BED encoding: 2-bit code → dosage of A1 allele
# 00 = hom ref (A2/A2) → 0
# 01 = missing          → nan
# 10 = het              → 1
# 11 = hom alt (A1/A1)  → 2
_DECODE = np.array([0, np.nan, 1, 2], dtype=np.float32)


# ── Argument parsing ───────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bfile",          required=True,
                   help="PLINK bfile prefix (.bed/.bim/.fam)")
    p.add_argument("--out",            required=True,
                   help="Output directory (basename must contain '1kg' or 'ukbb')")
    p.add_argument("--chroms",         nargs="+", type=int, required=True,
                   help="Chromosome number(s) to process")
    p.add_argument("--n-panel",        type=int, default=None,
                   help="Individuals to subsample (default: use all)")
    p.add_argument("--seed",           type=int, default=42,
                   help="Random seed for subsampling [default: 42]")
    p.add_argument("--max-block-size", type=int, default=400,
                   help="Maximum SNPs per LD block; caps DP search width "
                        "[default: 400]")
    p.add_argument("--maf",            type=float, default=0.01,
                   help="Minor allele frequency filter [default: 0.01]")
    p.add_argument("--hm3-grch38",     default=None, metavar="FILE",
                   help="TSV with HapMap3 GRCh38 positions (columns: chrom pos rsid). "
                        "When provided, SNPs are restricted to HM3 sites and their "
                        "rsIDs are used in snpinfo and HDF5 snplists, so the custom "
                        "reference is compatible with run_prscs.py rsID remapping.")
    return p.parse_args()


# ── PLINK I/O ──────────────────────────────────────────────────────────────────

def read_fam(path):
    return [line.split()[1] for line in open(path)]


def read_bim(path):
    """Return list of dicts with global_idx, chrom, rsid, pos, a1 (ALT), a0 (REF)."""
    snps = []
    with open(path) as f:
        for i, line in enumerate(f):
            parts = line.split()
            snps.append(dict(global_idx=i, chrom=int(parts[0].lstrip("chr")),
                             rsid=parts[1], pos=int(parts[3]),
                             a1=parts[4], a0=parts[5]))
    return snps


def load_genotypes(bed_path, chrom_snps, ind_idx, n_total_inds):
    """Load genotype matrix (individuals × SNPs) for one chromosome.

    Reads SNPs one at a time by seeking to each row in the BED file.
    Missing values are mean-imputed per SNP.
    """
    bytes_per_snp = (n_total_inds + 3) // 4
    G = np.empty((len(ind_idx), len(chrom_snps)), dtype=np.float32)

    with open(bed_path, "rb") as f:
        magic = f.read(3)
        if magic != b'\x6c\x1b\x01':
            raise ValueError(
                f"{bed_path}: not a valid SNP-major PLINK BED file "
                f"(magic bytes: {magic.hex()})"
            )
        for j, snp in enumerate(chrom_snps):
            f.seek(3 + snp["global_idx"] * bytes_per_snp)
            raw      = np.frombuffer(f.read(bytes_per_snp), dtype=np.uint8)
            bits     = np.unpackbits(raw, bitorder="little")[: n_total_inds * 2]
            bits     = bits.reshape(n_total_inds, 2)
            codes    = bits[:, 0] + 2 * bits[:, 1]
            geno_all = _DECODE[codes]
            g        = geno_all[ind_idx]
            mean_g   = np.nanmean(g)
            g[np.isnan(g)] = 0.0 if np.isnan(mean_g) else mean_g
            G[:, j]  = g

    return G


# ── MDL block detection ────────────────────────────────────────────────────────

def standardize(G):
    """Return float64 zero-mean, unit-variance standardized genotype matrix.

    Uses in-place operations to avoid allocating intermediate arrays.
    Peak memory is one float32 copy (input) + one float64 copy (output).
    """
    G   = G.astype(np.float64)
    mu  = G.mean(axis=0)
    sd  = G.std(axis=0)
    sd[sd == 0] = 1.0
    G  -= mu
    G  /= sd
    return G


def precompute_banded_corr(G_std, max_block_size):
    """Compute banded correlation matrix without storing the full m×m matrix.

    R_band[i, d] = corr(SNP i, SNP i+d) for d in 1..max_block_size.

    Uses vectorised dot products: for each offset d,
        R_band[:m-d, d] = col-wise dot of G_std[:,i] and G_std[:,i+d] / n
    which is O(n * m) per offset and O(n * m * max_block_size) overall.
    """
    n, m    = G_std.shape
    R_band  = np.zeros((m, max_block_size + 1), dtype=np.float64)
    # d=0 is the diagonal (all 1s); we leave it at zero as it's unused.
    for d in range(1, min(max_block_size + 1, m)):
        R_band[:m - d, d] = np.einsum("ni,ni->i",
                                      G_std[:, :m - d],
                                      G_std[:, d:]) / n
    return R_band


def mdl_block_detection(R_band, n_ind, max_block_size):
    """Find optimal LD block boundaries via MDL dynamic programming.

    MDL cost for a block of k SNPs starting at position j:
        MDL(j, k) = n/2 * log|R_{j:j+k}| + k*(k-1)/4 * log(n)

    log|R| is accumulated incrementally using the Schur complement:
        log|R_{j:j+k}| = log|R_{j:j+k-1}| + log(1 - c^T R_{j:j+k-1}^{-1} c)
    where c = [R[j, j+k-1], R[j+1, j+k-1], ..., R[j+k-2, j+k-1]]
            = R_band[j:j+k-1, k-1:-1:-1]   (read from banded matrix)

    R^{-1} is also updated via the Schur complement (matrix inversion lemma):
        R_new^{-1} = block([[R_old^{-1} + vv^T/s, -v/s],
                             [         -v^T/s,      1/s]])
    where v = R_old^{-1} c  and  s = 1 - c^T v.

    Parameters
    ----------
    R_band : ndarray, shape (m, max_block_size+1)
        Banded correlation matrix; R_band[i, d] = corr(SNP i, SNP i+d).
    n_ind : int
        Number of individuals (used in the MDL penalty).
    max_block_size : int
        Maximum allowed block size (caps the DP search width).

    Returns
    -------
    block_starts : list of int
        0-based indices where each block begins (first entry is always 0).
    """
    m     = R_band.shape[0]
    log_n = np.log(n_ind)

    # ── Precompute MDL cost table ──────────────────────────────────────────────
    # cost[j, k] = MDL cost of block starting at SNP j with k SNPs.
    # k=1 → cost = 0 (single SNP has log|1×1|=0 and zero penalty).
    cost = np.full((m, max_block_size + 1), np.inf, dtype=np.float64)

    for j in range(m):
        cost[j, 1] = 0.0
        if j + 1 >= m:
            continue

        R_inv   = np.ones((1, 1), dtype=np.float64)
        log_det = 0.0

        for k in range(2, min(max_block_size + 1, m - j + 1)):
            # c = correlations of existing block j..j+k-2 with new SNP j+k-1
            # R_band[i, d] = corr(i, i+d), so corr(j+i, j+k-1) = R_band[j+i, k-1-i]
            rows    = np.arange(j, j + k - 1)
            offsets = np.arange(k - 1, 0, -1)        # k-1, k-2, ..., 1
            c       = R_band[rows, offsets]            # shape (k-1,)

            v = R_inv @ c
            s = max(1.0 - float(c @ v), 1e-12)        # Schur complement (clamped)
            log_det += np.log(s)

            # Schur complement update of R^{-1}
            R_inv_new                  = np.empty((k, k), dtype=np.float64)
            R_inv_new[:k-1, :k-1]     = R_inv + np.outer(v, v) / s
            R_inv_new[:k-1, k-1]      = -v / s
            R_inv_new[k-1, :k-1]      = -v / s
            R_inv_new[k-1,  k-1]      = 1.0 / s
            R_inv = R_inv_new

            penalty    = k * (k - 1) / 4 * log_n
            cost[j, k] = n_ind / 2 * log_det + penalty

        if j % 500 == 0:
            print(f"    MDL cost table: SNP {j}/{m} ...",
                  file=sys.stderr, flush=True)

    # ── Dynamic programming ────────────────────────────────────────────────────
    dp     = np.full(m + 1, np.inf, dtype=np.float64)
    parent = np.zeros(m + 1, dtype=np.int32)
    dp[0]  = 0.0

    for i in range(1, m + 1):
        for j in range(max(0, i - max_block_size), i):
            val = dp[j] + cost[j, i - j]
            if val < dp[i]:
                dp[i]     = val
                parent[i] = j

    # ── Traceback ──────────────────────────────────────────────────────────────
    block_starts = []
    i = m
    while i > 0:
        block_starts.append(int(parent[i]))
        i = parent[i]
    block_starts.reverse()

    return block_starts


def compute_ld_blocks_mdl(G, snps, max_block_size):
    """Detect LD blocks via MDL and return per-block correlation matrices.

    Returns
    -------
    ld_blocks  : list of float32 ndarrays (one correlation matrix per block)
    snp_blocks : list of lists of rsID strings (one list per block)
    """
    n_ind   = G.shape[0]
    m       = G.shape[1]
    G_std   = standardize(G)
    del G   # free float32 matrix; G_std (float64) is all we need from here on

    print(f"    Precomputing banded correlation (max offset={max_block_size}) ...",
          file=sys.stderr, flush=True)
    R_band = precompute_banded_corr(G_std, max_block_size)

    print(f"    Running MDL dynamic programming ({m} SNPs, n={n_ind}) ...",
          file=sys.stderr, flush=True)
    block_starts = mdl_block_detection(R_band, n_ind, max_block_size)
    block_starts.append(m)   # sentinel

    n_blocks = len(block_starts) - 1
    print(f"    MDL found {n_blocks} blocks "
          f"(mean size {m / n_blocks:.1f} SNPs)",
          file=sys.stderr, flush=True)

    # Compute correlation matrix for each block
    ld_blocks, snp_blocks = [], []
    for a, b in zip(block_starts[:-1], block_starts[1:]):
        G_blk = G_std[:, a:b]
        R_blk = (G_blk.T @ G_blk) / n_ind
        R_blk = (R_blk + R_blk.T) / 2
        np.clip(R_blk, -1.0, 1.0, out=R_blk)
        np.fill_diagonal(R_blk, 1.0)
        ld_blocks.append(R_blk.astype(np.float32))
        snp_blocks.append([snps[i]["rsid"] for i in range(a, b)])

    return ld_blocks, snp_blocks


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    basename = os.path.basename(args.out.rstrip("/\\"))
    if "1kg" not in basename and "ukbb" not in basename:
        print(
            f"WARNING: output directory basename '{basename}' does not contain "
            "'1kg' or 'ukbb'. PRScs.py will raise an error when loading this "
            "reference. Rename the directory or patch parse_ldblk().",
            file=sys.stderr,
        )

    # Build (chrom_str, pos) → rsID lookup from HM3 GRCh38 table if provided.
    # This replaces positional BIM IDs with proper rsIDs in the output, making
    # the custom reference compatible with the rsID remapping in run_prscs.py.
    hm3_lookup = {}
    if args.hm3_grch38:
        hm3 = pd.read_csv(args.hm3_grch38, sep="\t")
        hm3["chrom"] = hm3["chrom"].astype(str).str.lstrip("chr")
        hm3["pos"]   = hm3["pos"].astype(int)
        hm3_lookup   = dict(zip(zip(hm3["chrom"], hm3["pos"]), hm3["rsid"]))
        print(f"Loaded {len(hm3_lookup)} HM3 GRCh38 positions from {args.hm3_grch38}",
              file=sys.stderr)

    all_snps = read_bim(args.bfile + ".bim")
    iids     = read_fam(args.bfile + ".fam")
    n_total  = len(iids)

    rng = np.random.default_rng(args.seed)
    if args.n_panel is not None and args.n_panel < n_total:
        ind_idx = np.sort(rng.choice(n_total, size=args.n_panel, replace=False))
        print(f"Subsampled {args.n_panel} / {n_total} individuals (seed={args.seed})",
              file=sys.stderr)
    else:
        ind_idx = np.arange(n_total)
        print(f"Using all {n_total} individuals", file=sys.stderr)

    snpinfo_rows = []

    for chrom in sorted(args.chroms):
        chrom_snps = [s for s in all_snps if s["chrom"] == chrom]
        if not chrom_snps:
            print(f"chr{chrom}: no SNPs in BIM — skipping", file=sys.stderr)
            continue

        # If an HM3 lookup was provided, restrict to HM3 SNPs and replace BIM IDs
        # with rsIDs BEFORE loading genotypes — avoids loading the full chromosome
        # into memory (e.g. 100k SNPs → ~5k HM3 SNPs).
        if hm3_lookup:
            hm3_snps = []
            for s in chrom_snps:
                rsid = hm3_lookup.get((str(chrom), s["pos"]))
                if rsid:
                    hm3_snps.append(dict(s, rsid=rsid))
            print(f"chr{chrom}: {len(hm3_snps)} / {len(chrom_snps)} SNPs have HM3 rsIDs",
                  file=sys.stderr, flush=True)
            if not hm3_snps:
                print(f"chr{chrom}: no HM3 SNPs found — skipping", file=sys.stderr)
                continue
            chrom_snps = hm3_snps

        print(f"chr{chrom}: loading {len(chrom_snps)} SNPs ...",
              file=sys.stderr, flush=True)
        G = load_genotypes(args.bfile + ".bed", chrom_snps, ind_idx, n_total)

        # MAF filter
        alt_freq = G.mean(axis=0) / 2.0
        maf      = np.minimum(alt_freq, 1.0 - alt_freq)
        keep     = maf >= args.maf
        if not keep.any():
            print(f"chr{chrom}: no SNPs pass MAF >= {args.maf} — skipping",
                  file=sys.stderr)
            continue

        G             = G[:, keep]
        chrom_snps_k  = [s for s, k in zip(chrom_snps, keep) if k]
        maf_k         = maf[keep]
        print(f"chr{chrom}: {keep.sum()} / {len(keep)} SNPs pass MAF >= {args.maf}",
              file=sys.stderr, flush=True)

        ld_blocks, snp_blocks = compute_ld_blocks_mdl(
            G, chrom_snps_k, args.max_block_size
        )

        hdf5_path = os.path.join(args.out, f"ldblk_1kg_chr{chrom}.hdf5")
        with h5py.File(hdf5_path, "w") as hf:
            for b_i, (ldblk, snplist) in enumerate(zip(ld_blocks, snp_blocks), 1):
                grp = hf.create_group(f"blk_{b_i}")
                grp.create_dataset("ldblk",   data=ldblk)
                grp.create_dataset("snplist",
                                   data=[s.encode("UTF-8") for s in snplist])
        print(f"chr{chrom}: wrote {hdf5_path} ({len(ld_blocks)} blocks)",
              file=sys.stderr)

        for snp, maf_val in zip(chrom_snps_k, maf_k):
            snpinfo_rows.append((chrom, snp["rsid"], snp["pos"],
                                 snp["a1"], snp["a0"], float(maf_val)))

    snpinfo_path = os.path.join(args.out, "snpinfo_1kg_hm3")
    with open(snpinfo_path, "w") as f:
        f.write("CHR\tSNP\tBP\tA1\tA2\tMAF\n")
        for row in snpinfo_rows:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t"
                    f"{row[4]}\t{row[5]:.6f}\n")
    print(f"Wrote {snpinfo_path} ({len(snpinfo_rows)} SNPs)", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
