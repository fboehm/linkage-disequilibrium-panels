#!/usr/bin/env Rscript
# plot_convergence.R
#
# Diagnostic plots from the convergence metrics emitted by run_ldpred2.R and
# run_prscs.py. All values are averaged across scenarios (trait, h2, p_causal,
# effect_dist) within each (panel_ancestry, panel_n) cell.
#
# Usage:
#   Rscript scripts/plot_convergence.R \
#       results/all_metrics.tsv \
#       results/plots/convergence

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

args   <- commandArgs(trailingOnly = TRUE)
input  <- if (length(args) >= 1) args[1] else "results/all_metrics.tsv"
outdir <- if (length(args) >= 2) args[2] else "results/plots/convergence"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11))

ancestry_colours <- c(
  hapnest_AMR = "#E69F00",
  hapnest_AFR = "#56B4E9",
  hapnest_EUR = "#009E73",
  hapnest_EAS = "#0072B2",
  hapnest_CSA = "#D55E00",
  hapnest_MID = "#CC79A7"
)

d <- read.delim(input, stringsAsFactors = FALSE)
d <- d[d$panel_ancestry %in% names(ancestry_colours), ]

panel_breaks <- sort(unique(d$panel_n))

save_plot <- function(p, name, w = 10, h = 6) {
  path <- file.path(outdir, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 150)
  message("Wrote ", path)
}

# ── 1. LDpred2 chain-convergence rate vs panel size ──────────────────────────

conv_wide <- d |>
  filter(method == "ldpred2",
         metric %in% c("n_chains_converged", "n_chains_total")) |>
  pivot_wider(names_from = metric, values_from = value) |>
  mutate(conv_rate = n_chains_converged / n_chains_total) |>
  group_by(panel_ancestry, panel_n) |>
  summarise(mean_rate = mean(conv_rate, na.rm = TRUE),
            n_runs    = sum(!is.na(conv_rate)),
            .groups   = "drop")

if (nrow(conv_wide) > 0) {
  p1 <- ggplot(conv_wide, aes(panel_n, mean_rate,
                              colour = panel_ancestry,
                              group  = panel_ancestry)) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_log10(breaks = panel_breaks, labels = comma) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    scale_colour_manual(values = ancestry_colours) +
    labs(title  = "LDpred2: fraction of converged chains vs LD-panel size",
         x      = "LD-panel size N (log scale)",
         y      = "Mean fraction of chains converged",
         colour = "Panel ancestry")
  save_plot(p1, "01_ldpred2_chain_convergence")
}

# ── 2. LDpred2 fallback-path frequency vs panel size ─────────────────────────

fallback <- d |>
  filter(method == "ldpred2", metric == "chains_used_fallback") |>
  group_by(panel_ancestry, panel_n) |>
  summarise(fallback_rate = mean(value, na.rm = TRUE),
            n_runs        = sum(!is.na(value)),
            .groups       = "drop")

if (nrow(fallback) > 0) {
  p2 <- ggplot(fallback, aes(panel_n, fallback_rate,
                             colour = panel_ancestry,
                             group  = panel_ancestry)) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_log10(breaks = panel_breaks, labels = comma) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    scale_colour_manual(values = ancestry_colours) +
    labs(title    = paste("LDpred2: fraction of runs hitting the",
                          "all-chains-diverged fallback"),
         subtitle = "1.0 = all runs at this (ancestry, N) had zero converged chains",
         x        = "LD-panel size N (log scale)",
         y        = "Fraction of runs",
         colour   = "Panel ancestry")
  save_plot(p2, "02_ldpred2_fallback_rate")
}

# ── 3. Sparsity of the fitted PGS vs panel size (both methods) ──────────────

sparsity <- d |>
  filter(metric %in% c("n_snps_nonzero", "n_snps_total",
                       "n_snps_in_output")) |>
  pivot_wider(names_from = metric, values_from = value) |>
  mutate(
    denom    = ifelse(is.na(n_snps_total), n_snps_in_output, n_snps_total),
    nz_frac  = n_snps_nonzero / denom
  ) |>
  filter(!is.na(nz_frac), is.finite(nz_frac)) |>
  group_by(method, panel_ancestry, panel_n) |>
  summarise(mean_nz_frac = mean(nz_frac, na.rm = TRUE),
            .groups      = "drop")

if (nrow(sparsity) > 0) {
  p3 <- ggplot(sparsity, aes(panel_n, mean_nz_frac,
                             colour = panel_ancestry,
                             group  = panel_ancestry)) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~method) +
    scale_x_log10(breaks = panel_breaks, labels = comma) +
    scale_y_continuous(labels = percent) +
    scale_colour_manual(values = ancestry_colours) +
    labs(title  = "Fraction of SNPs with non-zero effect (|β| > 1e-12) vs panel size",
         x      = "LD-panel size N (log scale)",
         y      = "Mean fraction of non-zero β",
         colour = "Panel ancestry")
  save_plot(p3, "03_nonzero_betas")
}
