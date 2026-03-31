#!/usr/bin/env python3
"""
Simulate admixed genotypes using msprime + stdpopsim AmericanAdmixture_4B11.

Outputs three gzipped VCF files per chromosome:
  gwas_{chrom}.vcf.gz   -- GWAS training population  (n_gwas diploid individuals)
  target_{chrom}.vcf.gz -- Independent evaluation population (n_target)
  panel_{chrom}.vcf.gz  -- Matched-admixed LD reference panel (n_panel)

The admixed population is the last population in the AmericanAdmixture_4B11
demographic model (Browning et al. 2011), which combines African, European,
and Native American ancestral components to approximate Hispanics/Latinos.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import stdpopsim


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chrom", required=True,
                   help="Chromosome label as used in stdpopsim, e.g. chr1")
    p.add_argument("--n-gwas", type=int, default=15000,
                   help="Number of diploid GWAS individuals (default: 15000)")
    p.add_argument("--n-target", type=int, default=5000,
                   help="Number of diploid target individuals (default: 5000)")
    p.add_argument("--n-panel", type=int, default=10000,
                   help="Number of diploid matched-panel individuals (default: 10000)")
    p.add_argument("--seed", type=int, required=True,
                   help="Random seed for the simulation")
    p.add_argument("--outdir", required=True,
                   help="Directory where output VCF files are written")
    p.add_argument("--length-multiplier", type=float, default=1.0,
                   help="Fraction of chromosome length to simulate. "
                        "Use a value < 1 (e.g. 0.1) for quick test runs.")
    p.add_argument("--hapmap3-sites", default=None,
                   help="Path to a tab-delimited file (with header) containing "
                        "HapMap3 SNP positions. Required columns: 'chrom' (integer, "
                        "no 'chr' prefix) and 'pos' (1-based GRCh37 bp). When "
                        "provided, random mutations are replaced by mutations "
                        "injected at exactly these positions.")
    return p.parse_args()


def get_contig(species, chrom, length_multiplier):
    """Load contig with HapMapII recombination map, falling back to flat rate."""
    try:
        return species.get_contig(
            chrom,
            genetic_map="HapMapII_GRCh37",
            length_multiplier=length_multiplier,
        )
    except Exception as exc:
        print(
            f"Warning: could not load HapMapII genetic map ({exc}). "
            "Falling back to flat recombination rate.",
            file=sys.stderr,
        )
        return species.get_contig(chrom, length_multiplier=length_multiplier)


def drop_multiroot_sites(ts):
    """Remove sites that fall in trees with more than one root.

    Incomplete coalescence leaves some genomic intervals with a forest of
    disconnected subtrees rather than a single tree.  tskit's write_vcf
    accesses tree.root (singular) internally, raising ValueError for such
    intervals.  Dropping the affected sites is safe: they cannot be assigned
    unambiguous ancestral/derived alleles anyway.
    """
    sites_to_drop = [
        site.id
        for tree in ts.trees()
        if tree.num_roots > 1
        for site in tree.sites()
    ]
    if sites_to_drop:
        print(
            f"  Dropping {len(sites_to_drop)} site(s) in multi-root trees",
            file=sys.stderr,
        )
        return ts.delete_sites(sites_to_drop)
    return ts


def write_cohort_vcf(ts, ind_ids, cohort_name, chrom, outdir):
    """Simplify the tree sequence to a cohort's individuals and write a VCF."""
    # Collect haplotype node IDs for these individuals
    nodes = []
    for i in ind_ids:
        nodes.extend(ts.individual(i).nodes.tolist())

    sub_ts = ts.simplify(samples=nodes)
    sub_ts = drop_multiroot_sites(sub_ts)

    # Name individuals consistently: {cohort}_{0-based index}
    ind_names = [f"{cohort_name}_{i}" for i in range(sub_ts.num_individuals)]

    # stdpopsim/msprime uses numeric contig IDs internally; strip "chr" prefix
    contig_id = chrom.lstrip("chr")

    vcf_path = outdir / f"{cohort_name}_{chrom}.vcf.gz"
    with gzip.open(str(vcf_path), "wt") as fh:
        sub_ts.write_vcf(fh, individual_names=ind_names, contig_id=contig_id)

    print(
        f"  Wrote {vcf_path}  "
        f"({sub_ts.num_individuals} individuals, {sub_ts.num_sites} sites)",
        file=sys.stderr,
    )


def load_hapmap3_positions(sites_file, chrom):
    """Return sorted list of HapMap3 bp positions for *chrom* from *sites_file*.

    The file must be tab-delimited with a header row containing at least the
    columns 'chrom' and 'pos'.  The 'chrom' column may carry a 'chr' prefix
    or not; both are accepted.
    """
    chrom_num = chrom.lstrip("chr")
    positions = []
    with open(sites_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["chrom"].lstrip("chr") == chrom_num:
                positions.append(int(row["pos"]))
    return sorted(set(positions))


def inject_hapmap3_mutations(ts, positions, seed):
    """Replace all mutations in *ts* with ones placed at *positions*.

    For each requested position that falls within the simulated sequence, a
    mutation is added to a branch chosen at random proportional to branch
    length.  This preserves the demographic structure: the probability that
    a mutation lands on a given lineage equals the fraction of total tree
    height on that lineage, so derived-allele frequencies reflect the
    simulated history.

    Alleles are set to ancestral='A' / derived='T' (arbitrary biallelic
    coding; the downstream MAF filter keeps only polymorphic sites anyway).
    """
    rng = np.random.default_rng(seed)
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()

    seq_len = ts.sequence_length
    n_added = 0

    for pos in positions:
        pos_f = float(pos)
        if pos_f >= seq_len:
            continue                          # outside the simulated region

        tree = ts.at(pos_f)
        nodes = [n for n in tree.nodes() if n not in set(tree.roots)]
        if not nodes:
            continue

        branch_lengths = np.array([tree.branch_length(n) for n in nodes],
                                   dtype=float)
        total = branch_lengths.sum()
        if total == 0:
            continue

        branch_lengths /= total
        mut_node = int(rng.choice(nodes, p=branch_lengths))

        site_idx = tables.sites.add_row(position=pos_f, ancestral_state="A")
        tables.mutations.add_row(site=site_idx, node=mut_node, derived_state="T")
        n_added += 1

    tables.sort()
    print(
        f"  Injected mutations at {n_added} / {len(positions)} requested "
        f"HapMap3 sites (sequence length = {seq_len:.0f} bp)",
        file=sys.stderr,
    )
    return tables.tree_sequence()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Load demographic model ────────────────────────────────────────────────
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("AmericanAdmixture_4B11")

    print("AmericanAdmixture_4B11 populations:", file=sys.stderr)
    for pop in model.populations:
        print(f"  name={pop.name!r}  description={pop.description!r}",
              file=sys.stderr)

    # The admixed population (ADMIX) is conventionally the last in this model.
    # Verify by checking the description if needed.
    admix_pop = model.populations[-1].name
    print(f"Using admixed population: {admix_pop!r}", file=sys.stderr)

    # ── Simulate ──────────────────────────────────────────────────────────────
    n_total = args.n_gwas + args.n_target + args.n_panel
    contig = get_contig(species, args.chrom, args.length_multiplier)

    engine = stdpopsim.get_engine("msprime")
    ts = engine.simulate(
        demographic_model=model,
        contig=contig,
        samples={admix_pop: n_total},
        seed=args.seed,
    )

    print(
        f"Simulation complete: {ts.num_individuals} diploid individuals, "
        f"{ts.num_sites} variant sites, "
        f"sequence length {ts.sequence_length:.0f} bp",
        file=sys.stderr,
    )

    # ── Optionally replace mutations with HapMap3-anchored sites ─────────────
    if args.hapmap3_sites:
        print(f"Loading HapMap3 positions from: {args.hapmap3_sites}",
              file=sys.stderr)
        positions = load_hapmap3_positions(args.hapmap3_sites, args.chrom)
        print(f"  {len(positions)} positions on {args.chrom}", file=sys.stderr)
        ts = inject_hapmap3_mutations(ts, positions, seed=args.seed)

    # ── Split individuals into cohorts ────────────────────────────────────────
    # Individuals are assigned in simulation order: GWAS first, then target,
    # then matched panel. This is deterministic given the seed.
    cohorts = {
        "gwas":   list(range(0, args.n_gwas)),
        "target": list(range(args.n_gwas, args.n_gwas + args.n_target)),
        "panel":  list(range(args.n_gwas + args.n_target, n_total)),
    }

    for cohort_name, ind_ids in cohorts.items():
        write_cohort_vcf(ts, ind_ids, cohort_name, args.chrom, outdir)


if __name__ == "__main__":
    main()
