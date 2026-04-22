#!/bin/bash
# Run Snakemake on the login node — it submits HAPNEST genotype jobs to SLURM.
#
# Usage:
#   screen -S ld-panels        # or: tmux new -s ld-panels
#   bash run_slurm.sh
#   Ctrl-A D                   # detach (screen) / Ctrl-B D (tmux)
#
# Re-attach later with:
#   screen -r ld-panels        # or: tmux attach -t ld-panels

set -euo pipefail

# ── Environment ────────────────────────────────────────────────────────────────
# Initialize conda before loading any modules so that module load snakemake
# can call conda deactivate/activate without hitting an uninitialised shell.
source "$(conda info --base)/etc/profile.d/conda.sh"

# Update the module name below if needed: run `module avail snakemake` to check.
module load snakemake

cd /scratch/jacks.local/frederick.boehm/linkage-disequilibrium-panels

mkdir -p logs/slurm

# ── Run Snakemake ──────────────────────────────────────────────────────────────
# Target: all HAPNEST VCF outputs (gwas/target/panel) for every chromosome and
# replicate.  Change 'all_hapnest' to 'all_pgs' to run the full PGS pipeline.

export TMPDIR=/scratch/jacks.local/frederick.boehm/tmp
mkdir -p $TMPDIR

#snakemake all_pgs \
results/evaluation/prscs/hapnest_public/rep1/hapnest_AMR/n500/quantitative/h2_0.4/pc_0.001/gaussian/metrics.tsv \
    --profile profile/slurm \
    --config sim_methods="[hapnest_public]" \
    --rerun-incomplete \
    --keep-going
