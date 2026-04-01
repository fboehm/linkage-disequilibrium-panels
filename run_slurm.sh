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
# Update the module name below if needed: run `module avail snakemake` to check.
module load snakemake

# Initialize conda for non-interactive shells (needed when running as a script).
source "$(conda info --base)/etc/profile.d/conda.sh"

cd /scratch/jacks.local/frederick.boehm/linkage-disequilibrium-panels

mkdir -p logs/slurm

# ── Run Snakemake ──────────────────────────────────────────────────────────────
# Target: all HAPNEST VCF outputs (gwas/target/panel) for every chromosome and
# replicate.  Change 'all_hapnest' to 'all_pgs' to run the full PGS pipeline.

snakemake all_pgs \
    --profile profile/slurm \
    --config sim_methods="[hapnest]" \
    --rerun-incomplete \
    --keep-going
