#!/usr/bin/env Rscript
# Plot PGS accuracy vs panel sample size, one line per panel ancestry,
# one PNG per (trait, effect_dist, h2, method) combination.
#
# Accuracy metric: R2_incremental for quantitative trait,
#                  AUC_incremental for binary traits.
# Within each PNG: colour = panel ancestry, facet = p_causal.
#
# Usage:
#   Rscript scripts/plot_accuracy_by_ancestry.R \
#       results/summary/all_metrics.tsv \
#       results/summary/plots

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

args   <- commandArgs(trailingOnly = TRUE)
input  <- if (length(args) >= 1) args[1] else "results/summary/all_metrics.tsv"
outdir <- if (length(args) >= 2) args[2] else "results/summary/plots"

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
d <- d[d$h2 >= 0.2, ]

acc <- d |>
  filter(
    (startsWith(trait, "binary") & metric == "AUC_incremental") |
    (trait == "quantitative"     & metric == "R2_incremental")
  ) |>
  filter(panel_ancestry %in% names(ancestry_colours)) |>
  mutate(method = factor(method, levels = c("ldpred2", "prscs")))

if (nrow(acc) == 0) stop("No matching rows in ", input)

agg <- acc |>
  group_by(trait, effect_dist, h2, p_causal,
           method, panel_ancestry, panel_n) |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            n_rep    = sum(!is.na(value)),
            .groups  = "drop")

panel_breaks <- sort(unique(agg$panel_n))

combos <- agg |>
  distinct(trait, effect_dist, h2, method) |>
  arrange(trait, effect_dist, h2, method)

cat(sprintf("[plot] %d (trait, effect_dist, h2, method) combos to render\n",
            nrow(combos)))

for (i in seq_len(nrow(combos))) {
  s   <- combos[i, ]
  sub <- agg |>
    filter(trait       == s$trait,
           effect_dist == s$effect_dist,
           h2          == s$h2,
           method      == s$method)

  y_label <- if (s$trait == "quantitative")
               expression(Incremental ~ R^2)
             else
               "Incremental AUC"

  title <- sprintf("%s  |  %s  |  h²=%g  |  %s",
                   s$trait, s$effect_dist, s$h2, s$method)

  p <- ggplot(sub, aes(panel_n, mean_val,
                       colour = panel_ancestry,
                       group  = panel_ancestry)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ p_causal, labeller = label_both) +
    scale_x_log10(breaks = panel_breaks, labels = comma) +
    scale_colour_manual(values = ancestry_colours) +
    labs(x = "LD-panel sample size",
         y = y_label,
         colour = "Panel ancestry",
         title  = title)

  out <- file.path(outdir,
                   sprintf("13_accuracy_vs_panelN_by_ancestry_%s_%s_h2-%g_%s.png",
                           s$trait, s$effect_dist, s$h2, s$method))
  ggsave(out, p, width = 8, height = 5.5, dpi = 150)
  cat(sprintf("[plot] wrote %s\n", out))
}
