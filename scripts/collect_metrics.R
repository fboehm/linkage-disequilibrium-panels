#!/usr/bin/env Rscript
# collect_metrics.R
#
# Walk the results/evaluation directory tree and aggregate all per-scenario
# metrics.tsv files into a single wide table.
#
# Expected path structure:
#   results/evaluation/{method}/rep{rep}/{panel_ancestry}/n{panel_n}/
#       {trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/metrics.tsv
#
# Usage:
#   Rscript scripts/collect_metrics.R \
#       --results-dir results/evaluation \
#       --out         results/summary/all_metrics.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--results-dir", type = "character",
              default = "results/evaluation",
              help = "Root of evaluation directory tree [default: %default]"),
  make_option("--out", type = "character",
              help = "Output TSV path")
)

opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$out))
  stop("--out is required")

# ── Discover all metrics files ─────────────────────────────────────────────────

files <- list.files(opt$`results-dir`,
                    pattern   = "^metrics\\.tsv$",
                    recursive = TRUE,
                    full.names = TRUE)

cat(sprintf("[collect_metrics] Found %d metrics.tsv files\n", length(files)))
if (length(files) == 0L) {
  warning("No metrics.tsv files found; writing empty output.")
  fwrite(data.table(), opt$out, sep = "\t")
  quit(status = 0L)
}

# ── Parse path components ──────────────────────────────────────────────────────
# Path relative to results-dir:
#   {method}/rep{rep}/{panel_ancestry}/n{panel_n}/
#   {trait}/h2_{h2}/pc_{p_causal}/{effect_dist}/metrics.tsv

parse_path <- function(full_path, root) {
  rel <- sub(paste0("^", normalizePath(root, mustWork = FALSE), "/"), "",
             normalizePath(full_path, mustWork = FALSE))
  parts <- strsplit(rel, "/", fixed = TRUE)[[1L]]
  # parts: method / rep{N} / ancestry / n{N} / trait / h2_{V} / pc_{V} / dist / metrics.tsv
  if (length(parts) != 9L) return(NULL)
  list(
    method         = parts[1L],
    rep            = as.integer(sub("rep", "", parts[2L])),
    panel_ancestry = parts[3L],
    panel_n        = as.integer(sub("n", "", parts[4L])),
    trait          = parts[5L],
    h2             = as.numeric(sub("h2_", "", parts[6L])),
    p_causal       = as.numeric(sub("pc_", "", parts[7L])),
    effect_dist    = parts[8L]
  )
}

rows <- lapply(files, function(f) {
  meta <- parse_path(f, opt$`results-dir`)
  if (is.null(meta)) {
    warning("Unexpected path structure, skipping: ", f)
    return(NULL)
  }
  m <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0L) return(NULL)
  cbind(as.data.table(meta), m)
})

rows <- Filter(Negate(is.null), rows)
if (length(rows) == 0L) {
  warning("All metrics files were empty or unparseable.")
  fwrite(data.table(), opt$out, sep = "\t")
  quit(status = 0L)
}

out <- rbindlist(rows, fill = TRUE)
setcolorder(out, c("method", "rep", "panel_ancestry", "panel_n",
                   "trait", "h2", "p_causal", "effect_dist"))
setorder(out, method, panel_ancestry, panel_n, rep, trait, h2, p_causal, effect_dist)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
fwrite(out, opt$out, sep = "\t")
cat(sprintf("[collect_metrics] Wrote %d rows to %s\n", nrow(out), opt$out))
