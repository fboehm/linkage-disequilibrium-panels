#!/usr/bin/env Rscript
# run_ldpred2.R
#
# Compute LDPred2-auto polygenic score weights from GWAS summary statistics
# using a pre-computed LD reference (panel FBM + SFBM from prepare_ldpred2_ref.R).
#
# Inputs
#   --sumstats   : PLINK2 --glm output (linear or logistic.hybrid), TEST==ADD rows
#   --ref-dir    : directory with shared panel.rds, ld_sfbm.rds, matched_snps.tsv
#   --pheno      : phenotype file (FID IID Y) — used to count cases/controls
#   --train-ids  : training-set keep-file (FID IID, no header)
#   --trait-type : quantitative | binary
#   --n-train    : number of individuals in the GWAS training set
#   --seed       : random seed
#   --out        : output TSV with per-SNP columns: SNP, A1, BETA
#
# Dependencies: bigsnpr, bigsparser, optparse, data.table

suppressPackageStartupMessages({
  library(bigsnpr)
  library(bigsparser)
  library(optparse)
  library(data.table)
})

# ── Argument parsing ──────────────────────────────────────────────────────────

opt_list <- list(
  make_option("--sumstats",   type = "character",
              help = "PLINK2 GLM summary-statistics file"),
  make_option("--ref-dir",    type = "character",
              help = "Directory with shared LD reference (panel.rds, ld_sfbm.rds, matched_snps.tsv)"),
  make_option("--pheno",      type = "character",
              help = "Phenotype file (FID IID Y) for case/control counting"),
  make_option("--train-ids",  type = "character",
              help = "Training-set keep-file (FID IID, no header)"),
  make_option("--trait-type", type = "character", default = "quantitative",
              help = "quantitative | binary [default: %default]"),
  make_option("--n-train",    type = "integer",
              help = "Total training-set sample size"),
  make_option("--seed",       type = "integer",   default = 1L,
              help = "Random seed [default: %default]"),
  make_option("--ncores",     type = "integer",   default = 1L,
              help = "Number of cores for LDSC and LDpred2-auto [default: %default]"),
  make_option("--out",        type = "character",
              help = "Output TSV: SNP  A1  BETA")
)

opt <- parse_args(OptionParser(option_list = opt_list))

missing_args <- Filter(is.null, list(
  sumstats  = opt$sumstats,
  `ref-dir` = opt$`ref-dir`,
  `n-train` = opt$`n-train`,
  out       = opt$out
))
if (length(missing_args) > 0)
  stop("Required arguments missing: ", paste(names(missing_args), collapse = ", "))

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

ref_dir <- opt$`ref-dir`

# ── Load pre-computed LD reference ────────────────────────────────────────────

cat(sprintf("[ldpred2] Loading shared LD reference from: %s\n", ref_dir))

# Pre-computed sparse LD matrix
sfbm_rds <- file.path(ref_dir, "ld_sfbm.rds")
if (!file.exists(sfbm_rds))
  stop("ld_sfbm.rds not found in --ref-dir: ", ref_dir)
corr <- readRDS(sfbm_rds)
cat(sprintf("[ldpred2] LD matrix: %d × %d (loaded from cache)\n",
            nrow(corr), ncol(corr)))

# Pre-computed matched-SNP table
msnps_file <- file.path(ref_dir, "matched_snps.tsv")
if (!file.exists(msnps_file))
  stop("matched_snps.tsv not found in --ref-dir: ", ref_dir)
matched_snps <- fread(msnps_file)
rename_map <- c("_NUM_ID_" = "NUM_ID", "_FLIP_" = "FLIP")
old_nm <- intersect(names(rename_map), names(matched_snps))
if (length(old_nm) > 0) setnames(matched_snps, old_nm, rename_map[old_nm])
# Keep only needed columns
if (!"sfbm_row" %in% names(matched_snps))
  stop("matched_snps.tsv is missing the sfbm_row column.")

cat(sprintf("[ldpred2] Matched SNP reference: %d SNPs\n", nrow(matched_snps)))

# ── Read and normalise GWAS summary statistics ────────────────────────────────

cat(sprintf("[ldpred2] Reading summary statistics: %s\n", opt$sumstats))
ss <- fread(opt$sumstats, data.table = TRUE)
setnames(ss, 1L, "CHROM")
ss <- ss[TEST == "ADD"]

# plink2 --glm linear  → columns BETA, SE
# plink2 --glm logistic hybrid → columns OR, LOG(OR)_SE
if ("BETA" %in% names(ss)) {
  ss <- ss[!is.na(BETA) & !is.na(SE) & SE > 0]
} else if ("OR" %in% names(ss)) {
  setnames(ss, "LOG(OR)_SE", "SE")
  ss[, BETA := log(OR)]
  ss <- ss[!is.na(BETA) & !is.na(SE) & SE > 0 & is.finite(BETA)]
} else {
  stop("[ldpred2] Cannot find BETA or OR column in summary statistics.")
}
ss[, CHROM := as.integer(sub("^chr", "", CHROM))]
ss[, POS   := as.integer(POS)]

cat(sprintf("[ldpred2] Summary statistics: %d SNPs after filtering\n", nrow(ss)))

# ── Restrict to HM3 SNPs ──────────────────────────────────────────────────────
# matched_snps defines the HM3 reference set; inner-join to keep only those SNPs.
hm3_pos <- matched_snps[, .(chr, pos)]
ss <- ss[hm3_pos, on = c("CHROM" = "chr", "POS" = "pos"), nomatch = NULL]
cat(sprintf("[ldpred2] After HM3 restriction: %d SNPs\n", nrow(ss)))

# ── Effective sample size ─────────────────────────────────────────────────────

if (opt$`trait-type` == "binary") {
  if (is.null(opt$pheno) || is.null(opt$`train-ids`))
    stop("--pheno and --train-ids are required for binary traits")

  pheno_df    <- read.table(opt$pheno,      header = TRUE)
  train_df    <- read.table(opt$`train-ids`, header = FALSE,
                             col.names = c("FID", "IID"))
  train_pheno <- merge(pheno_df, train_df, by = c("FID", "IID"))
  n_cases     <- sum(train_pheno$Y == 1L, na.rm = TRUE)
  n_controls  <- sum(train_pheno$Y == 0L, na.rm = TRUE)

  if (n_cases == 0L || n_controls == 0L)
    stop(sprintf("Training set has %d cases and %d controls.", n_cases, n_controls))

  n_eff <- 4 / (1 / n_cases + 1 / n_controls)
  cat(sprintf("[ldpred2] Binary trait: %d cases + %d controls  n_eff = %.1f\n",
              n_cases, n_controls, n_eff))
} else {
  n_eff <- opt$`n-train`
}
cat(sprintf("[ldpred2] Effective N: %.1f\n", n_eff))

# ── Align summary stats to pre-computed matched-SNP order ────────────────────
# The SFBM rows/columns correspond to matched_snps in sfbm_row order.
# We merge the current sumstats into this reference order.

ss_key <- ss[, .(chr = CHROM, pos = POS, ID, REF, A1, BETA, SE)]

# Merge on chr+pos (works with positional IDs and cross-build data)
merged <- ss_key[matched_snps, on = c("chr", "pos"), nomatch = NA]

# Deduplicate by sfbm_row: multi-allelic sites (two GWAS SNPs at same position)
# produce one extra row per duplicated position; keep the first match.
if (anyDuplicated(merged$sfbm_row)) {
  n_dup <- sum(duplicated(merged$sfbm_row))
  cat(sprintf("[ldpred2] WARNING: dropping %d duplicate sfbm_row(s) from multi-allelic sites\n", n_dup))
  merged <- merged[!duplicated(merged$sfbm_row), ]
}

# Determine per-SNP beta after potential allele flip
# If FLIP == TRUE (stored as 1 in TSV), the panel a0 == ss A1, so negate beta
if ("FLIP" %in% names(merged)) {
  merged[, beta_aligned := fifelse(FLIP == 1L, -BETA, BETA)]
} else {
  merged[, beta_aligned := BETA]
}

# Drop SNPs not found in current sumstats (set beta=0, large SE = no effect)
n_missing <- sum(is.na(merged$BETA))
if (n_missing > 0) {
  cat(sprintf("[ldpred2] WARNING: %d / %d reference SNPs not found in current sumstats; filling beta=0\n",
              n_missing, nrow(matched_snps)))
  merged[is.na(BETA), beta_aligned := 0]
  merged[is.na(SE),   SE           := 1e6]
}

# Sort by sfbm_row to match SFBM column order
setorder(merged, sfbm_row)

n_matched <- nrow(merged)
if (n_matched != nrow(corr))
  stop(sprintf("[ldpred2] Alignment error: df_beta has %d rows but SFBM has %d.",
               n_matched, nrow(corr)))

cat(sprintf("[ldpred2] Aligned %d SNPs to SFBM\n", n_matched))

# Assemble df_beta in SFBM row order
info_snp <- data.frame(
  chr     = merged$chr,
  pos     = merged$pos,
  a0      = matched_snps$a0,
  a1      = matched_snps$a1,
  rsid    = matched_snps$rsid,
  beta    = merged$beta_aligned,
  beta_se = merged$SE,
  n_eff   = n_eff,
  stringsAsFactors = FALSE
)

# ── LD score regression for h2 initialisation ────────────────────────────────

cat("[ldpred2] Running LDSC for h2 initialisation ...\n")
ldsc_res <- tryCatch(
  snp_ldsc2(corr, info_snp, ncores = opt$ncores),
  error = function(e) {
    message("[ldpred2] snp_ldsc2 failed (", conditionMessage(e), "); using h2_init = 0.1")
    list(h2 = 0.1)
  }
)
h2_init <- max(ldsc_res[["h2"]], 1e-4)
cat(sprintf("[ldpred2] LDSC h2 estimate: %.4f  (used as init)\n", h2_init))

# ── LDPred2-auto ──────────────────────────────────────────────────────────────

cat("[ldpred2] Running LDPred2-auto ...\n")

multi_auto <- snp_ldpred2_auto(
  corr            = corr,
  df_beta         = info_snp,
  h2_init         = h2_init,
  vec_p_init      = seq_log(1e-4, 0.9, length.out = 30L),
  burn_in         = 500L,
  num_iter        = 500L,
  report_step     = 20L,   # 500 / 20 = 25 stored samples per chain
  allow_jump_sign = FALSE,
  shrink_corr     = 0.95,
  ncores          = opt$ncores
)

# Filter to converged chains
converged <- sapply(multi_auto, function(auto)
  !any(is.na(auto$beta_est)) &&
  auto$h2_est  > 0.001 &&
  auto$h2_est  < 1.5
)
cat(sprintf("[ldpred2] %d / %d chains converged\n",
            sum(converged), length(converged)))

if (!any(converged)) converged[] <- TRUE

beta_mat <- sapply(multi_auto[converged], `[[`, "beta_est")
betas_final <- if (is.matrix(beta_mat)) rowMeans(beta_mat) else beta_mat

# Posterior variance from pooled per-iteration MCMC draws across converged chains.
# sample_beta is stored on the standardised scale (beta * sd); divide variance by sd^2.
sample_mats <- lapply(multi_auto[converged], `[[`, "sample_beta")
has_samples <- length(sample_mats) > 0L &&
  all(vapply(sample_mats,
             function(m) !is.null(m) && ncol(m) > 0L,
             logical(1)))

if (has_samples) {
  all_samples <- do.call(cbind, sample_mats)
  n_samp      <- ncol(all_samples)
  row_mean    <- Matrix::rowMeans(all_samples)
  row_sqmean  <- Matrix::rowMeans(all_samples * all_samples)
  var_std     <- (row_sqmean - row_mean^2) * n_samp / max(n_samp - 1L, 1L)
  sd_snp      <- 1 / sqrt(info_snp$n_eff * info_snp$beta_se^2 + info_snp$beta^2)
  betas_var   <- as.numeric(var_std) / sd_snp^2
  cat(sprintf("[ldpred2] Posterior variance from %d MCMC samples across %d converged chains\n",
              n_samp, sum(converged)))
} else {
  cat("[ldpred2] WARNING: sample_beta unavailable; falling back to between-chain variance\n")
  betas_var <- if (is.matrix(beta_mat))
    matrixStats::rowVars(beta_mat)
  else
    rep(NA_real_, length(beta_mat))
}

betas_final[is.na(betas_final)] <- 0
betas_var[is.na(betas_var)]     <- 0

h2_final <- mean(sapply(multi_auto[converged], `[[`, "h2_est"))
p_final  <- mean(sapply(multi_auto[converged], `[[`, "p_est"))
cat(sprintf("[ldpred2] Final estimates: h2 = %.4f  p_causal = %.5f\n",
            h2_final, p_final))

# ── Write per-SNP weights ─────────────────────────────────────────────────────

out_df <- data.frame(
  SNP      = paste0("chr", info_snp$chr, ":", info_snp$pos, ":",
                    info_snp$a0,  ":", info_snp$a1),
  A1       = info_snp$a1,
  BETA     = betas_final,
  BETA_VAR = betas_var,
  CHR      = info_snp$chr,
  BP       = info_snp$pos,
  stringsAsFactors = FALSE
)
write.table(out_df,
            file      = opt$out,
            sep       = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote     = FALSE)

cat(sprintf("[ldpred2] Wrote %d SNP weights to: %s\n", nrow(out_df), opt$out))
