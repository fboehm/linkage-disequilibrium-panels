#!/bin/bash
#SBATCH --job-name=ld-panels-snakemake
#SBATCH --partition=compute
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/slurm/snakemake_%j.out
#SBATCH --error=logs/slurm/snakemake_%j.err

# ── Environment ────────────────────────────────────────────────────────────────
# Update the module name below if needed: run `module avail snakemake` to check.
module load snakemake

set -euo pipefail

cd /home/jacks.local/frederick.boehm/linkage-disequilibrium-panels

mkdir -p logs/slurm

# ── Run Snakemake ──────────────────────────────────────────────────────────────
# Requires: pip install snakemake-executor-plugin-slurm
#
# Change 'all_pgs' to 'all' to stop after genotype + phenotype generation.

snakemake all_pgs \
    --profile profile/slurm \
    --rerun-incomplete \
    --keep-going
