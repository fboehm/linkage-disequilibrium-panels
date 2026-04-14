#!/usr/bin/env bash
# cleanup_obsolete_results.sh
#
# Remove results directories that are obsolete after the study redesign:
#   - sim_methods other than hapnest_public (msprime, hapgen2, hapnest)
#   - panel ancestries other than hapnest_* (1KG panels, matched_admixed,
#     oracle, gwas_subset)
#
# Run from the root of the repository on the cluster:
#   bash scripts/cleanup_obsolete_results.sh
#
# Pass --dry-run to print what would be deleted without removing anything.

set -euo pipefail

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=true
  echo "[dry-run] No files will be deleted."
fi

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

remove() {
  local path="$1"
  if [[ -e "$path" ]]; then
    if $DRY_RUN; then
      echo "[dry-run] would remove: $path"
    else
      echo "Removing: $path"
      rm -rf "$path"
    fi
  fi
}

# ── 1. Obsolete sim_method top-level trees ────────────────────────────────────

OBSOLETE_SIM_METHODS=(msprime hapgen2 hapnest)

RESULT_ROOTS=(
  results/vcf
  results/vcf_filtered
  results/plink
  results/phenotypes
  results/splits
  results/pca
  results/gwas
  results/ldpred2_work
  results/prscs_custom_ref
)

for root in "${RESULT_ROOTS[@]}"; do
  for sm in "${OBSOLETE_SIM_METHODS[@]}"; do
    remove "$root/$sm"
  done
done

# pgs_weights and evaluation have an extra method-level directory
for method in ldpred2 prscs; do
  for sm in "${OBSOLETE_SIM_METHODS[@]}"; do
    remove "results/pgs_weights/$method/$sm"
    remove "results/evaluation/$method/$sm"
  done
done

# ── 2. Obsolete panel ancestry subtrees under hapnest_public ─────────────────

OBSOLETE_PANEL_ANCESTRIES=(
  AFR_1kg
  AMR_1kg
  EUR_1kg
  matched_admixed
  oracle
  gwas_subset
)

# Directories that contain panel_ancestry as a subdirectory level
PANEL_ROOTS=(
  results/ldpred2_work/hapnest_public
  results/prscs_custom_ref/hapnest_public
)

for root in "${PANEL_ROOTS[@]}"; do
  if [[ -d "$root" ]]; then
    for rep_dir in "$root"/rep*/; do
      [[ -d "$rep_dir" ]] || continue
      for pa in "${OBSOLETE_PANEL_ANCESTRIES[@]}"; do
        remove "${rep_dir}${pa}"
      done
    done
  fi
done

# pgs_weights: results/pgs_weights/{method}/hapnest_public/rep*/{panel_ancestry}
# evaluation:  results/evaluation/{method}/hapnest_public/rep*/{panel_ancestry}
for method in ldpred2 prscs; do
  for subdir in \
    "results/pgs_weights/$method/hapnest_public" \
    "results/evaluation/$method/hapnest_public"
  do
    if [[ -d "$subdir" ]]; then
      for rep_dir in "$subdir"/rep*/; do
        [[ -d "$rep_dir" ]] || continue
        for pa in "${OBSOLETE_PANEL_ANCESTRIES[@]}"; do
          remove "${rep_dir}${pa}"
        done
      done
    fi
  done
done

# ── 3. Obsolete 1KG panel genotype resources ─────────────────────────────────

for pa in AFR_1kg AMR_1kg EUR_1kg; do
  remove "resources/panels/$pa"
done

# ── Done ──────────────────────────────────────────────────────────────────────

if $DRY_RUN; then
  echo "[dry-run] Done. Re-run without --dry-run to delete."
else
  echo "Done."
fi
