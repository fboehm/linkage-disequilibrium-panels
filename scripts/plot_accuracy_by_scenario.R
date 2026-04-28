#!/usr/bin/env Rscript
# plot_accuracy_by_scenario.R
#
# One PNG per (trait, h2, p_causal, effect_dist) scenario.
#   x = LD-panel size (log scale)
#   y = accuracy (R2_incremental for quantitative, AUC_incremental for binary)
#   one trace per panel ancestry; methods (ldpred2 / prscs) shown via linetype
#
# Usage:
#   Rscript scripts/plot_accuracy_by_scenario.R \
#       results/all_metrics.tsv \
#       results/plots/by_scenario

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

args   <- commandArgs(trailingOnly = TRUE)
input  <- if (length(args) >= 1) args[1] else "results/all_metrics.tsv"
outdir <- if (length(args) >= 2) args[2] else "results/plots/by_scenario"

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

# Pick one accuracy metric per trait family.
acc <- d |>
  filter(
    (startsWith(trait, "binary") & metric == "AUC_incremental") |
    (trait == "quantitative"     & metric == "R2_incremental")
  ) |>
  mutate(method = factor(method, levels = c("ldpred2", "prscs")))

if (nrow(acc) == 0) stop("No matching accuracy rows in ", input)

# Average over replicates within each (scenario × method × ancestry × panel_n).
agg <- acc |>
  group_by(trait, h2, p_causal, effect_dist,
           method, panel_ancestry, panel_n) |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            se_val   = sd(value,   na.rm = TRUE) / sqrt(sum(!is.na(value))),
            n_rep    = sum(!is.na(value)),
            .groups  = "drop")

panel_breaks <- sort(unique(agg$panel_n))

scenarios <- agg |>
  distinct(trait, h2, p_causal, effect_dist) |>
  arrange(trait, h2, p_causal, effect_dist)

cat(sprintf("[plot] %d scenarios to render\n", nrow(scenarios)))

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  sub <- agg |>
    filter(trait       == s$trait,
           h2          == s$h2,
           p_causal    == s$p_causal,
           effect_dist == s$effect_dist)

  y_label <- if (s$trait == "quantitative") "R²" else "AUC"

  title <- sprintf("%s | h²=%g | p_causal=%g | %s",
                   s$trait, s$h2, s$p_causal, s$effect_dist)

  p <- ggplot(sub, aes(panel_n, mean_val,
                       colour   = panel_ancestry,
                       linetype = method,
                       shape    = method,
                       group    = interaction(panel_ancestry, method))) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_val - se_val,
                      ymax = mean_val + se_val),
                  width = 0, alpha = 0.4) +
    scale_x_log10(breaks = panel_breaks, labels = comma) +
    scale_colour_manual(values = ancestry_colours) +
    labs(title    = title,
         x        = "LD-panel size N (log scale)",
         y        = y_label,
         colour   = "Panel ancestry",
         linetype = "Method",
         shape    = "Method")

  fname <- sprintf("%s_h2-%s_pc-%s_%s.png",
                   s$trait,
                   sub("\\.", "p", format(s$h2,       trim = TRUE)),
                   sub("\\.", "p", format(s$p_causal, trim = TRUE)),
                   s$effect_dist)
  path <- file.path(outdir, fname)
  ggsave(path, p, width = 9, height = 6, dpi = 150)
  message("Wrote ", path)
}
