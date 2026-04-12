#!/usr/bin/env Rscript
# explore_results.R
#
# Exploratory plots of results/summary/all_metrics.tsv
#
# Usage:
#   Rscript scripts/explore_results.R \
#       results/summary/all_metrics.tsv \
#       results/summary/plots

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

args   <- commandArgs(trailingOnly = TRUE)
input  <- if (length(args) >= 1) args[1] else "results/summary/all_metrics.tsv"
outdir <- if (length(args) >= 2) args[2] else "results/summary/plots"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(input, stringsAsFactors = FALSE)

# ── helpers ────────────────────────────────────────────────────────────────────

save_plot <- function(p, name, w = 10, h = 7) {
  path <- file.path(outdir, paste0(name, ".png"))
  ggsave(path, p, width = w, height = h, dpi = 150)
  message("Wrote ", path)
}

theme_set(theme_bw(base_size = 11))

# Make panel_n and h2 ordered factors for clean axis ordering
d <- d |>
  mutate(
    panel_n_f = factor(panel_n, levels = sort(unique(panel_n))),
    h2_f      = factor(h2,      levels = sort(unique(h2))),
    method    = factor(method,  levels = c("ldpred2", "prscs"))
  )

ancestry_colours <- c(
  AFR_1kg         = "#E41A1C",
  AMR_1kg         = "#FF7F00",
  EUR_1kg         = "#377EB8",
  matched_admixed = "#4DAF4A",
  oracle          = "#984EA3",
  gwas_subset     = "#A65628"
)

# ── 1. AUC vs panel size (binary traits, averaged over replicates & p_causal) ──

auc <- d |>
  filter(metric == "AUC_adjusted", trait %in% c("binary_prev10", "binary_prev30")) |>
  group_by(method, panel_ancestry, panel_n_f, panel_n, h2_f, trait, effect_dist) |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val   = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups  = "drop")

p1 <- ggplot(auc, aes(panel_n, mean_val,
                      colour = panel_ancestry,
                      linetype = method,
                      shape = method)) +
  geom_line() +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                width = 0, alpha = 0.4) +
  facet_grid(trait + h2_f ~ effect_dist,
             labeller = label_both) +
  scale_x_log10(breaks = c(100, 500, 1000, 5000, 15000),
                labels  = comma) +
  scale_colour_manual(values = ancestry_colours) +
  labs(title  = "AUC (adjusted) vs LD-panel size",
       x      = "Panel N (log scale)",
       y      = "Mean AUC (± 1 SE)",
       colour = "Panel ancestry",
       linetype = "Method", shape = "Method") +
  theme(legend.position = "right",
        strip.text       = element_text(size = 7))

save_plot(p1, "01_auc_vs_panelN", w = 14, h = 18)

# ── 2. R² vs panel size (quantitative trait) ───────────────────────────────────

r2 <- d |>
  filter(metric == "R2_full", trait == "quantitative") |>
  group_by(method, panel_ancestry, panel_n_f, panel_n, h2_f, effect_dist) |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val   = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups  = "drop")

p2 <- ggplot(r2, aes(panel_n, mean_val,
                     colour = panel_ancestry,
                     linetype = method,
                     shape = method)) +
  geom_line() +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                width = 0, alpha = 0.4) +
  facet_grid(h2_f ~ effect_dist,
             labeller = label_both) +
  scale_x_log10(breaks = c(100, 500, 1000, 5000, 15000),
                labels  = comma) +
  scale_colour_manual(values = ancestry_colours) +
  labs(title  = "R² (full model) vs LD-panel size — quantitative trait",
       x      = "Panel N (log scale)",
       y      = "Mean R² (± 1 SE)",
       colour = "Panel ancestry",
       linetype = "Method", shape = "Method") +
  theme(legend.position = "right")

save_plot(p2, "02_r2_vs_panelN", w = 12, h = 10)

# ── 3. Method comparison: ldpred2 vs prscs (scatter, one point per scenario) ──

wide_method <- d |>
  filter(metric %in% c("AUC_adjusted", "R2_full")) |>
  tidyr::pivot_wider(
    id_cols     = c(sim_method, rep, panel_ancestry, panel_n, trait, h2, p_causal, effect_dist, metric),
    names_from  = method,
    values_from = value
  )

p3 <- ggplot(wide_method |> filter(!is.na(ldpred2), !is.na(prscs)),
             aes(ldpred2, prscs, colour = panel_ancestry)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.4, size = 1) +
  facet_grid(trait ~ metric, scales = "free",
             labeller = label_both) +
  scale_colour_manual(values = ancestry_colours) +
  labs(title  = "LDpred2 vs PRScs — per-scenario accuracy",
       x      = "LDpred2",
       y      = "PRScs",
       colour = "Panel ancestry") +
  theme(legend.position = "bottom")

save_plot(p3, "03_method_scatter", w = 10, h = 10)

# ── 4. Accuracy by ancestry (box plots, best large-panel scenarios only) ───────

large <- d |> filter(panel_n >= 5000, metric %in% c("AUC_adjusted", "R2_full"))

p4 <- ggplot(large, aes(panel_ancestry, value,
                        fill = method)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
  facet_grid(metric ~ trait, scales = "free_y",
             labeller = label_both) +
  scale_fill_manual(values = c(ldpred2 = "#1F78B4", prscs = "#33A02C")) +
  labs(title  = "Accuracy by panel ancestry (panel N ≥ 5 000)",
       x      = NULL, y      = "Metric value",
       fill   = "Method") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(p4, "04_accuracy_by_ancestry", w = 12, h = 8)

# ── 5. Effect of heritability on R²  ──────────────────────────────────────────

r2_h2 <- d |>
  filter(metric == "R2_full", trait == "quantitative", panel_n >= 1000) |>
  group_by(method, panel_ancestry, h2_f, effect_dist) |>
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

p5 <- ggplot(r2_h2, aes(h2_f, mean_val,
                        fill = panel_ancestry)) +
  geom_col(position = position_dodge(0.85)) +
  facet_grid(effect_dist ~ method,
             labeller = label_both) +
  scale_fill_manual(values = ancestry_colours) +
  labs(title  = "Mean R² by heritability (panel N ≥ 1 000, quantitative)",
       x      = "h²", y = "Mean R²",
       fill   = "Panel ancestry") +
  theme(legend.position = "right")

save_plot(p5, "05_r2_by_h2", w = 10, h = 8)

# ── 6. Incremental R² (PC-adjusted) vs panel size ─────────────────────────────

r2_inc <- d |>
  filter(metric == "R2_incremental", trait == "quantitative") |>
  group_by(method, panel_ancestry, panel_n, h2_f) |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val   = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups  = "drop")

p6 <- ggplot(r2_inc, aes(panel_n, mean_val,
                         colour = panel_ancestry,
                         linetype = method)) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~h2_f, labeller = label_both) +
  scale_x_log10(breaks = c(100, 500, 1000, 5000, 15000),
                labels  = comma) +
  scale_colour_manual(values = ancestry_colours) +
  labs(title  = "Incremental R² (PC-adjusted) vs panel size — quantitative",
       x      = "Panel N (log scale)", y = "Mean incremental R²",
       colour = "Panel ancestry", linetype = "Method")

save_plot(p6, "06_incremental_r2_vs_panelN", w = 12, h = 8)

# ── 7. PGS posterior variance vs panel size ───────────────────────────────────
# matched_admixed, quantitative trait, gaussian effect distribution
# colour = h2, linetype/shape = p_causal, facet = method

var_d <- d |>
  filter(
    metric          == "pgs_var_mean",
    panel_ancestry  == "matched_admixed",
    trait           == "quantitative",
    effect_dist     == "gaussian",
    method          == "prscs"
  ) |>
  mutate(p_causal_f = factor(p_causal)) |>
  group_by(panel_n, h2_f, p_causal_f) |>
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

p7 <- ggplot(var_d, aes(panel_n, mean_val,
                        colour   = h2_f,
                        linetype = p_causal_f,
                        shape    = p_causal_f,
                        group    = interaction(h2_f, p_causal_f))) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_log10(breaks = c(100, 500, 1000, 5000, 10000),
                labels  = comma) +
  labs(title    = "PGS posterior variance vs panel size — PRScs\n(matched_admixed, quantitative, gaussian)",
       x        = "Panel N (log scale)",
       y        = "Mean PGS variance",
       colour   = "h²",
       linetype = "p_causal",
       shape    = "p_causal") +
  theme(legend.position = "right")

save_plot(p7, "07_pgs_variance_vs_panelN", w = 12, h = 6)

# ── 8. Accuracy vs panel size: matched_admixed only ───────────────────────────

admixed_vals <- d |>
  filter(
    panel_ancestry == "matched_admixed",
    metric %in% c("R2_full", "AUC_adjusted")
  ) |>
  group_by(method, panel_n, metric, trait, h2_f, effect_dist) |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val   = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups  = "drop")

p8_base <- function(metric_sel, trait_sel, ylabel) {
  ggplot(
    admixed_vals |> filter(metric == metric_sel, trait == trait_sel),
    aes(panel_n, mean_val, colour = method, shape = method)
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_grid(h2_f ~ effect_dist, labeller = label_both) +
    scale_x_log10(breaks = c(100, 500, 1000, 5000, 10000),
                  labels  = comma) +
    scale_colour_manual(values = c(ldpred2 = "#1F78B4", prscs = "#33A02C")) +
    labs(x      = "Panel N (log scale)",
         y      = ylabel,
         colour = "Method", shape = "Method") +
    theme(legend.position = "right",
          strip.text       = element_text(size = 8))
}

p8a <- p8_base("R2_full", "quantitative", "Mean R²") +
  ggtitle("R² vs panel size: matched_admixed — quantitative")

p8b <- p8_base("AUC_adjusted", "binary_prev10", "Mean AUC") +
  ggtitle("AUC vs panel size: matched_admixed — binary (prev = 0.10)")

p8c <- p8_base("AUC_adjusted", "binary_prev30", "Mean AUC") +
  ggtitle("AUC vs panel size: matched_admixed — binary (prev = 0.30)")

save_plot(p8a, "08a_r2_admixed",        w = 12, h = 10)
save_plot(p8b, "08b_auc_admixed_prev10", w = 12, h = 10)
save_plot(p8c, "08c_auc_admixed_prev30", w = 12, h = 10)

# ── 9. R² vs panel size — all settings (quantitative) ────────────────────────
# One figure per effect_dist; rows = h2, cols = p_causal; no averaging over
# p_causal so every simulation setting gets its own panel.

r2_all <- d |>
  filter(metric == "R2_full", trait == "quantitative") |>
  mutate(p_causal_f = factor(p_causal, levels = sort(unique(p_causal)))) |>
  group_by(method, panel_ancestry, panel_n_f, panel_n, h2_f, p_causal_f, effect_dist) |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val   = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups  = "drop")

for (ed in sort(unique(r2_all$effect_dist))) {
  p9 <- ggplot(r2_all |> filter(effect_dist == ed),
               aes(panel_n, mean_val,
                   colour   = panel_ancestry,
                   linetype = method,
                   shape    = method)) +
    geom_line() +
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  width = 0, alpha = 0.4) +
    facet_grid(h2_f ~ p_causal_f,
               labeller = labeller(h2_f      = label_both,
                                   p_causal_f = label_both)) +
    scale_x_log10(breaks = c(100, 500, 1000, 5000, 15000),
                  labels  = comma) +
    scale_colour_manual(values = ancestry_colours) +
    labs(title    = sprintf("R\u00b2 (full model) vs LD-panel size — quantitative, %s", ed),
         subtitle = "rows = h\u00b2 | columns = p_causal",
         x        = "Panel N (log scale)",
         y        = "Mean R\u00b2 (\u00b1 1 SE)",
         colour   = "Panel ancestry",
         linetype = "Method", shape = "Method") +
    theme(legend.position = "right",
          strip.text       = element_text(size = 7))

  save_plot(p9, sprintf("09_r2_all_settings_%s", ed), w = 16, h = 12)
}

message("Done — ", length(list.files(outdir, "*.png")), " plots written to ", outdir)
