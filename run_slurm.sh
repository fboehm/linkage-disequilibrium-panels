#!/bin/bash
# Run Snakemake on the login node — it submits jobs to SLURM itself.
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

cd /home/jacks.local/frederick.boehm/linkage-disequilibrium-panels

mkdir -p logs/slurm

# ── Run Snakemake ──────────────────────────────────────────────────────────────
# Change 'all_pgs' to 'all' to stop after genotype + phenotype generation.

snakemake all_pgs \
    --profile profile/slurm \
    --use-envmodules \
    --rerun-incomplete \
    --keep-going
