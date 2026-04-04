#!/usr/bin/env Rscript
# evaluate_pgs.R
#
# Compute prediction accuracy of a polygenic score in the held-out test set.
#
# Metrics:
#   quantitative traits  →  R² (proportion of phenotypic variance explained)
#   binary traits        →  AUC (area under the ROC curve)
#
# Usage:
#   Rscript scripts/evaluate_pgs.R \
#       --scores     results/pgs/rep1/.../scores.sscore \
#       --pheno      results/phenotypes/rep1/.../pheno.pheno \
#       --test-ids   results/splits/rep1/test.txt \
#       --trait-type quantitative \
#       --out        results/evaluation/rep1/.../metrics.tsv

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
  make_option("--scores",      type = "character",
              help = "PLINK2 .sscore file"),
  make_option("--pheno",       type = "character",
              help = "Phenotype file with columns FID IID Y"),
  make_option("--test-ids",    type = "character",
              help = "Test-set keep-file (FID IID, no header)"),
  make_option("--trait-type",  type = "character",
              help = "quantitative | binary"),
  make_option("--covariates",  type = "character", default = NULL,
              help = "PLINK2 .eigenvec file of PCs (optional). When provided, "
                     "incremental metrics (adjusted for PCs) are also reported."),
  make_option("--var-scores",  type = "character", default = NULL,
              help = "TSV with FID IID PGS_VAR (from score_pgs_variance.R). "
                     "When provided, adds mean/median per-individual PGS "
                     "posterior variance to the metrics output."),
  make_option("--out",         type = "character",
              help = "Output TSV of accuracy metrics")
)

opt <- parse_args(OptionParser(option_list = opt_list))

missing_args <- Filter(is.null, list(
  scores      = opt$scores,
  pheno       = opt$pheno,
  `test-ids`  = opt$`test-ids`,
  `trait-type`= opt$`trait-type`,
  out         = opt$out
))
if (length(missing_args) > 0)
  stop("Required arguments missing: ", paste(names(missing_args), collapse = ", "))

# ── Read data ─────────────────────────────────────────────────────────────────

# plink2 .sscore: #FID IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM SCORE1_AVG
scores <- read.table(opt$scores, header = TRUE, check.names = FALSE)
colnames(scores)[1] <- "FID"           # strip leading '#'
scores <- scores[, c("FID", "IID", "SCORE1_AVG"), drop = FALSE]

pheno <- read.table(opt$pheno, header = TRUE)   # FID IID Y

test_ids <- read.table(opt$`test-ids`, header = FALSE,
                       col.names = c("FID", "IID"))

# ── Restrict to test-set individuals and merge ────────────────────────────────

df <- merge(scores, pheno,    by = c("FID", "IID"))
df <- merge(df,     test_ids, by = c("FID", "IID"))  # inner join → test set only
df <- df[!is.na(df$Y) & !is.na(df$SCORE1_AVG), ]

# Optionally merge principal components
pc_cols <- character(0)
if (!is.null(opt$covariates)) {
  pcs <- read.table(opt$covariates, header = TRUE, check.names = FALSE)
  colnames(pcs)[1] <- "FID"   # strip leading '#'
  pc_cols <- grep("^PC", colnames(pcs), value = TRUE)
  df <- merge(df, pcs[, c("FID", "IID", pc_cols)], by = c("FID", "IID"))
  df <- df[complete.cases(df[, pc_cols]), ]
}

if (nrow(df) == 0L)
  stop("No individuals remain after merging scores, phenotypes, and test IDs.")

cat(sprintf("[evaluate_pgs] %d test individuals after merge\n", nrow(df)))

# ── Compute metric ────────────────────────────────────────────────────────────

wilcoxon_auc <- function(scores_vec, labels) {
  n1 <- sum(labels == 1L)
  n0 <- sum(labels == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  ranks    <- rank(scores_vec)
  rank_sum <- sum(ranks[labels == 1L])
  (rank_sum - n1 * (n1 + 1) / 2) / (n1 * n0)
}

if (opt$`trait-type` == "quantitative") {

  fit_null <- if (length(pc_cols) > 0L)
    lm(as.formula(paste("Y ~", paste(pc_cols, collapse = " + "))), data = df)
  else
    lm(Y ~ 1, data = df)

  fit_full <- if (length(pc_cols) > 0L)
    lm(as.formula(paste("Y ~ SCORE1_AVG +", paste(pc_cols, collapse = " + "))), data = df)
  else
    lm(Y ~ SCORE1_AVG, data = df)

  r2_null <- summary(fit_null)$r.squared
  r2_full <- summary(fit_full)$r.squared
  r2_incr <- r2_full - r2_null

  metrics <- data.frame(
    metric       = c("R2_full", "R2_null", "R2_incremental"),
    value        = c(r2_full,   r2_null,   r2_incr),
    n            = nrow(df)
  )
  cat(sprintf("[evaluate_pgs] R²_full = %.4f  R²_null = %.4f  R²_incr = %.4f  (n = %d)\n",
              r2_full, r2_null, r2_incr, nrow(df)))

} else {

  # Marginal AUC of the raw PGS (rank-based, no external packages)
  auc_raw <- wilcoxon_auc(df$SCORE1_AVG, df$Y)

  # PC-adjusted AUC: residualise PGS on PCs, then compute AUC of residuals
  if (length(pc_cols) > 0L) {
    pgs_resid <- residuals(
      lm(as.formula(paste("SCORE1_AVG ~", paste(pc_cols, collapse = " + "))), data = df)
    )
    auc_adj <- wilcoxon_auc(pgs_resid, df$Y)
  } else {
    auc_adj <- auc_raw
  }

  n1 <- sum(df$Y == 1L)
  n0 <- sum(df$Y == 0L)
  metrics <- data.frame(
    metric     = c("AUC_raw", "AUC_adjusted"),
    value      = c(auc_raw,   auc_adj),
    n          = nrow(df),
    n_cases    = n1,
    n_controls = n0
  )
  cat(sprintf("[evaluate_pgs] AUC_raw = %.4f  AUC_adj = %.4f  (n_cases = %d, n_controls = %d)\n",
              auc_raw, auc_adj, n1, n0))
}

# ── Per-individual PGS variance (optional) ───────────────────────────────────

if (!is.null(opt$`var-scores`)) {
  var_df <- read.table(opt$`var-scores`, header = TRUE)
  var_df <- merge(var_df, test_ids, by = c("FID", "IID"))
  var_df <- var_df[!is.na(var_df$PGS_VAR), ]
  if (nrow(var_df) > 0L) {
    pgs_var_mean   <- mean(var_df$PGS_VAR)
    pgs_var_median <- median(var_df$PGS_VAR)
    cat(sprintf("[evaluate_pgs] PGS posterior variance: mean = %.4g  median = %.4g\n",
                pgs_var_mean, pgs_var_median))
    var_rows <- data.frame(
      metric = c("pgs_var_mean", "pgs_var_median"),
      value  = c(pgs_var_mean,   pgs_var_median),
      n      = nrow(var_df)
    )
    # Align columns with the main metrics frame before rbinding
    for (col in setdiff(names(metrics), names(var_rows)))
      var_rows[[col]] <- NA
    metrics <- rbind(metrics, var_rows[, names(metrics)])
  }
}

# ── Write output ──────────────────────────────────────────────────────────────

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.table(metrics,
            file      = opt$out,
            sep       = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote     = FALSE)
cat(sprintf("[evaluate_pgs] Metrics written to: %s\n", opt$out))
