#!/usr/bin/env Rscript
# Plot quantitative-trait R2_incremental vs panel sample size,
# one line per h2, for the AMR panel only.
# Separate PNGs per (method, effect_dist); facet by p_causal.
#
# Usage:
#   Rscript scripts/plot_quant_r2_by_h2.R \
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

d <- read.delim(input, stringsAsFactors = FALSE)
d <- d[d$h2 >= 0.2, ]

acc <- d |>
  filter(trait == "quantitative",
         metric == "R2_incremental",
         panel_ancestry == "hapnest_AMR") |>
  mutate(h2 = factor(h2, levels = sort(unique(h2))))

if (nrow(acc) == 0) stop("No matching rows in ", input)

agg <- acc |>
  group_by(method, effect_dist, p_causal, h2, panel_n) |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            n_rep    = sum(!is.na(value)),
            .groups  = "drop")

panel_breaks <- sort(unique(agg$panel_n))

combos <- agg |>
  distinct(method, effect_dist) |>
  arrange(method, effect_dist)

for (i in seq_len(nrow(combos))) {
  m  <- combos$method[i]
  ed <- combos$effect_dist[i]
  sub <- agg |> filter(method == m, effect_dist == ed)

  p <- ggplot(sub, aes(panel_n, mean_val,
                       colour = h2, group = h2)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ p_causal, labeller = label_both) +
    scale_x_log10(breaks = panel_breaks, labels = comma) +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    labs(x = "LD-panel sample size",
         y = expression(Incremental ~ R^2),
         colour = expression(h^2),
         title = sprintf("Quantitative R² vs panel N  |  %s  |  %s  |  AMR panel",
                         m, ed))

  out <- file.path(outdir,
                   sprintf("11_quant_r2_vs_panelN_by_h2_AMR_%s_%s.png", m, ed))
  ggsave(out, p, width = 7, height = 5, dpi = 150)
  cat(sprintf("[plot] wrote %s\n", out))
}
