#!/usr/bin/env python3
"""
Simulate genotypes with msprime and write per-cohort bgzipped VCF files.

Replaces HAPGEN2 for environments where the HAPGEN2 binary (a 2011 static
Linux binary) crashes due to missing vsyscall support (e.g. WSL2).

Individuals are split sequentially:
  indices   0 .. n_gwas-1                  -> gwas_{chrom}.vcf.gz
  indices   n_gwas .. n_gwas+n_target-1    -> target_{chrom}.vcf.gz
  indices   n_gwas+n_target .. n_total-1   -> panel_{chrom}.vcf.gz

A constant-size coalescent model (Ne = hapgen2_ne) with the per-chromosome
Oxford IMPUTE2 genetic map is used for recombination rates, giving realistic
linkage-disequilibrium structure.
"""

import argparse
import gzip
import sys
from pathlib import Path

import msprime


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genmap",   required=True,
                   help="Oxford IMPUTE2 genetic map: position  rate(cM/Mb)  map(cM)")
    p.add_argument("--chrom",    required=True,
                   help="Chromosome label for VCF CHROM column (e.g. chr21)")
    p.add_argument("--n-gwas",   type=int, required=True)
    p.add_argument("--n-target", type=int, required=True)
    p.add_argument("--n-panel",  type=int, required=True)
    p.add_argument("--ne",       type=float, default=11418.0,
                   help="Effective population size (default: 11418, 1KG AMR)")
    p.add_argument("--mu",       type=float, default=1.25e-8,
                   help="Mutation rate per bp per generation (default: 1.25e-8)")
    p.add_argument("--seed",     type=int, required=True,
                   help="Base random seed (ancestry uses seed, mutations uses seed+1)")
    p.add_argument("--outdir",   required=True,
                   help="Output directory for VCF.gz files")
    return p.parse_args()


# ── Genetic map ───────────────────────────────────────────────────────────────

def load_rate_map(path):
    """Read Oxford IMPUTE2 genetic map and return an msprime.RateMap.

    File format (with header):
        position  COMBINED_rate(cM/Mb)  Genetic_Map(cM)

    The hapmap convention is that rate[i] applies from position[i] to
    position[i+1].  We extrapolate rate[0] back to position 0 and
    rate[-1] forward to sequence_length.
    """
    positions = []
    rates_per_bp = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts or not parts[0][0].isdigit():
                continue  # skip header / blank lines
            positions.append(float(parts[0]))
            rates_per_bp.append(float(parts[1]) * 1e-8)  # cM/Mb -> per bp

    if not positions:
        raise ValueError(f"No data rows found in genetic map: {path}")

    # Build breakpoints and per-interval rates for msprime.RateMap.
    # Breakpoints: [0, p0, p1, ..., pN-1, seq_len]  (N+2 values)
    # Rates:       [r0, r0, r1, ..., rN-1]           (N+1 values)
    #   -> rate r0 extrapolated to [0, p0]
    #   -> rate ri applies in [pi, p(i+1)]  (hapmap convention)
    #   -> rate rN-1 extrapolated past pN-1
    seq_len = positions[-1] + 1.0
    breakpoints = [0.0] + positions + [seq_len]
    rates = [rates_per_bp[0]] + rates_per_bp   # prepend r0 for the [0, p0] interval

    return msprime.RateMap(position=breakpoints, rate=rates)


# ── VCF helpers ───────────────────────────────────────────────────────────────

def vcf_header(chrom, sample_ids):
    return (
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={chrom}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(sample_ids) + "\n"
    )


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args    = parse_args()
    outdir  = Path(args.outdir)
    n_gwas  = args.n_gwas
    n_target = args.n_target
    n_panel = args.n_panel
    n_total = n_gwas + n_target + n_panel
    chrom   = args.chrom
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Genetic map ───────────────────────────────────────────────────────────
    print(f"Loading genetic map: {args.genmap}", file=sys.stderr, flush=True)
    rate_map = load_rate_map(args.genmap)
    print(f"  sequence length: {rate_map.sequence_length:,.0f} bp", file=sys.stderr, flush=True)

    # ── Demographic model (constant-size, matching HAPGEN2 -Ne) ───────────────
    demography = msprime.Demography()
    demography.add_population(name="AMR", initial_size=args.ne)

    # ── Ancestry simulation ───────────────────────────────────────────────────
    print(
        f"Simulating ancestry: {n_total} diploid individuals  "
        f"Ne={args.ne:.0f}  seed={args.seed}",
        file=sys.stderr, flush=True,
    )
    ts = msprime.sim_ancestry(
        samples=n_total,
        demography=demography,
        recombination_rate=rate_map,
        sequence_length=rate_map.sequence_length,
        random_seed=args.seed,
    )

    # ── Mutation simulation ───────────────────────────────────────────────────
    print(f"Adding mutations: mu={args.mu:.2e}  seed={args.seed + 1}", file=sys.stderr, flush=True)
    ts = msprime.sim_mutations(ts, rate=args.mu, random_seed=args.seed + 1)
    print(f"  {ts.num_sites} variant sites simulated", file=sys.stderr, flush=True)

    # ── Sample IDs ────────────────────────────────────────────────────────────
    sample_ids = [f"ind{i + 1}" for i in range(n_total)]
    gwas_ids   = sample_ids[:n_gwas]
    target_ids = sample_ids[n_gwas:n_gwas + n_target]
    panel_ids  = sample_ids[n_gwas + n_target:]

    # ── Open output VCF streams ───────────────────────────────────────────────
    fh_gwas   = gzip.open(outdir / f"gwas_{chrom}.vcf.gz",   "wt")
    fh_target = gzip.open(outdir / f"target_{chrom}.vcf.gz", "wt")
    fh_panel  = gzip.open(outdir / f"panel_{chrom}.vcf.gz",  "wt")

    fh_gwas.write(vcf_header(chrom, gwas_ids))
    fh_target.write(vcf_header(chrom, target_ids))
    fh_panel.write(vcf_header(chrom, panel_ids))

    # ── Stream variants ───────────────────────────────────────────────────────
    print("Writing VCF files ...", file=sys.stderr, flush=True)
    n_written = 0
    for var in ts.variants():
        alleles = var.alleles
        # Keep only biallelic SNPs with single-character alleles
        if len(alleles) != 2:
            continue
        ref, alt = alleles
        if alt is None or len(ref) != 1 or len(alt) != 1:
            continue

        pos    = int(var.site.position) + 1   # 0-based tskit -> 1-based VCF
        snp_id = f"{chrom}:{pos}"

        # var.genotypes is a flat haploid array of length 2*n_total
        gh  = var.genotypes
        gts = [f"{gh[2*i]}/{gh[2*i+1]}" for i in range(n_total)]

        fixed = f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t"
        fh_gwas.write(  fixed + "\t".join(gts[:n_gwas])                   + "\n")
        fh_target.write(fixed + "\t".join(gts[n_gwas:n_gwas + n_target])  + "\n")
        fh_panel.write( fixed + "\t".join(gts[n_gwas + n_target:])        + "\n")
        n_written += 1

    fh_gwas.close()
    fh_target.close()
    fh_panel.close()

    for name, ids in [("gwas", gwas_ids), ("target", target_ids), ("panel", panel_ids)]:
        print(
            f"  Wrote {name}_{chrom}.vcf.gz  "
            f"({len(ids)} individuals, {n_written} sites)",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
