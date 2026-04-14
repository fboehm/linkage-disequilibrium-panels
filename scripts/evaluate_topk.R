#!/usr/bin/env Rscript
# evaluate_topk.R
#
# For each threshold k, compute the fraction of top-k% PGS subjects who
# are truly in the top-k% of the trait distribution (precision at k),
# and the enrichment over random expectation (precision_at_k / k).
#
# Only meaningful for quantitative traits. Exits with an informative
# message and writes an empty file for binary traits.
#
# Usage:
#   Rscript scripts/evaluate_topk.R \
#       --scores     results/pgs/.../scores.sscore \
#       --pheno      results/phenotypes/.../pheno.pheno \
#       --test-ids   results/splits/rep1/test.txt \
#       --trait-type quantitative \
#       --k-levels   0.01,0.05,0.10,0.20,0.50 \
#       --out        results/evaluation/.../topk_metrics.tsv

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
  make_option("--scores",     type = "character",
              help = "PLINK2 .sscore file"),
  make_option("--pheno",      type = "character",
              help = "Phenotype file with columns FID IID Y"),
  make_option("--test-ids",   type = "character",
              help = "Test-set keep-file (FID IID, no header)"),
  make_option("--trait-type", type = "character",
              help = "quantitative | binary"),
  make_option("--k-levels",   type = "character",
              default = "0.01,0.05,0.10,0.20,0.50",
              help = "Comma-separated top-k thresholds in (0,1) [default: %default]"),
  make_option("--out",        type = "character",
              help = "Output TSV path")
)

opt <- parse_args(OptionParser(option_list = opt_list))

missing_args <- Filter(is.null, list(
  scores     = opt$scores,
  pheno      = opt$pheno,
  `test-ids` = opt$`test-ids`,
  `trait-type` = opt$`trait-type`,
  out        = opt$out
))
if (length(missing_args) > 0L)
  stop("Required arguments missing: ", paste(names(missing_args), collapse = ", "))

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)

# ── Guard: binary traits ───────────────────────────────────────────────────────

if (opt$`trait-type` != "quantitative") {
  cat("[evaluate_topk] Skipping binary trait — top-k not defined for 0/1 outcomes.\n")
  write.table(
    data.frame(metric = character(), value = numeric(), n = integer(),
               k = numeric(), n_top = integer(), n_overlap = integer()),
    opt$out, sep = "\t", row.names = FALSE, quote = FALSE
  )
  quit(status = 0L)
}

# ── Parse k levels ─────────────────────────────────────────────────────────────

k_levels <- as.numeric(strsplit(opt$`k-levels`, ",")[[1L]])
if (any(is.na(k_levels)) || any(k_levels <= 0) || any(k_levels >= 1))
  stop("--k-levels must be comma-separated values strictly in (0, 1)")

# ── Read inputs ────────────────────────────────────────────────────────────────

scores <- read.table(opt$scores, header = TRUE, check.names = FALSE,
                     comment.char = "")
colnames(scores)[1L] <- sub("^#", "", colnames(scores)[1L])
score_col <- intersect(c("SCORE1_AVG", "SCORE1_SUM", "SCORE1"), colnames(scores))
if (length(score_col) == 0L)
  stop("Cannot find score column. Columns: ",
       paste(colnames(scores), collapse = ", "))
colnames(scores)[colnames(scores) == score_col[1L]] <- "SCORE1_AVG"
scores <- scores[, c("FID", "IID", "SCORE1_AVG")]

pheno    <- read.table(opt$pheno,      header = TRUE)
test_ids <- read.table(opt$`test-ids`, header = FALSE,
                       col.names = c("FID", "IID"))

# ── Merge to test set ──────────────────────────────────────────────────────────

df <- merge(scores, pheno,    by = c("FID", "IID"))
df <- merge(df,     test_ids, by = c("FID", "IID"))
df <- df[!is.na(df$Y) & !is.na(df$SCORE1_AVG), ]
n  <- nrow(df)

if (n == 0L)
  stop("No individuals remain after merging scores, phenotypes, and test IDs.")

cat(sprintf("[evaluate_topk] %d test individuals after merge\n", n))

# ── Compute precision and enrichment at each k ─────────────────────────────────

rows <- do.call(rbind, lapply(k_levels, function(k) {
  n_top         <- max(1L, round(k * n))
  top_pgs_idx   <- order(df$SCORE1_AVG, decreasing = TRUE)[seq_len(n_top)]
  top_trait_idx <- order(df$Y,          decreasing = TRUE)[seq_len(n_top)]
  n_overlap     <- length(intersect(top_pgs_idx, top_trait_idx))
  precision     <- n_overlap / n_top
  enrichment    <- precision / k   # fold-enrichment over random expectation

  cat(sprintf(
    "  k = %.2f  n_top = %d  n_overlap = %d  precision = %.4f  enrichment = %.4f\n",
    k, n_top, n_overlap, precision, enrichment
  ))

  data.frame(
    metric    = c(sprintf("precision_at_k%.2f", k),
                  sprintf("enrichment_at_k%.2f", k)),
    value     = c(precision, enrichment),
    n         = n,
    k         = k,
    n_top     = n_top,
    n_overlap = c(n_overlap, n_overlap)
  )
}))

# ── Write output ───────────────────────────────────────────────────────────────

write.table(rows, opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("[evaluate_topk] Written to %s\n", opt$out))
