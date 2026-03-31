#!/usr/bin/env Rscript
# simulate_phenotypes.R
#
# Simulate phenotypes from a PLINK binary GWAS genotype file.
#
# Effect-size distributions:
#   gaussian   – ALL p SNPs have N(0, h²/p) effects (infinitesimal / purely polygenic)
#   spikeslab  – p_causal fraction get N(0, h²/(p·p_causal)) effects, rest = 0
#   scaledt    – p_causal fraction get scaled-t(df=4) effects, rest = 0
#
# Trait types:
#   quantitative – continuous Y on the liability scale
#   binary       – liability-threshold model; Y = 1 if liability > qnorm(1 - prevalence)
#
# Outputs (at --out prefix):
#   <prefix>.pheno            – tab-separated FID / IID / Y (PLINK2-compatible)
#   <prefix>.causal_snps.tsv  – per-causal-SNP true effect sizes and allele statistics
#   <prefix>.sim_params.tsv   – single-row summary of realised simulation parameters
#
# Dependencies: BEDMatrix, optparse
#
# Usage example:
#   Rscript scripts/simulate_phenotypes.R \
#       --bed         results/plink/gwas/rep1/merged.bed \
#       --h2          0.20 \
#       --p-causal    0.01 \
#       --effect-dist spikeslab \
#       --trait-type  quantitative \
#       --seed        42 \
#       --out         results/phenotypes/rep1/quantitative/h2_0.20/pc_0.01/spikeslab/pheno

suppressPackageStartupMessages({
  library(BEDMatrix)
  library(optparse)
})

# ── Argument parsing ──────────────────────────────────────────────────────────

opt_list <- list(
  make_option("--bed",         type = "character",
              help = "Path to PLINK .bed file (prefix or full path with .bed extension)"),
  make_option("--h2",          type = "double",
              help = "SNP heritability on the liability scale, in (0, 1)"),
  make_option("--p-causal",    type = "double",   default = 0.01,
              help = "Proportion of causal SNPs; ignored for --effect-dist gaussian [default: %default]"),
  make_option("--effect-dist", type = "character", default = "spikeslab",
              help = "Effect-size distribution: gaussian | spikeslab | scaledt [default: %default]"),
  make_option("--trait-type",  type = "character", default = "quantitative",
              help = "Trait type: quantitative | binary [default: %default]"),
  make_option("--prevalence",  type = "double",   default = 0.10,
              help = "Disease prevalence for binary traits [default: %default]"),
  make_option("--seed",        type = "integer",
              help = "Random seed (integer)"),
  make_option("--out",         type = "character",
              help = "Output file prefix (directory will be created if needed)")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# ── Validation ────────────────────────────────────────────────────────────────

missing_args <- Filter(is.null, list(bed = opt$bed, h2 = opt$h2,
                                     seed = opt$seed, out = opt$out))
if (length(missing_args) > 0) {
  stop("Required arguments missing: ", paste(names(missing_args), collapse = ", "))
}
stopifnot(
  "h2 must be in (0, 1)"          = opt$h2 > 0 && opt$h2 < 1,
  "p-causal must be in (0, 1]"    = opt$`p-causal` > 0 && opt$`p-causal` <= 1,
  "effect-dist not recognised"    = opt$`effect-dist` %in% c("gaussian", "spikeslab", "scaledt"),
  "trait-type not recognised"     = opt$`trait-type` %in% c("quantitative", "binary"),
  "prevalence must be in (0, 1)"  = opt$prevalence > 0 && opt$prevalence < 1
)

set.seed(opt$seed)

# ── Load genotype matrix ──────────────────────────────────────────────────────

bed_prefix <- sub("\\.bed$", "", opt$bed)
cat(sprintf("[simulate_phenotypes] Loading genotype matrix from: %s\n", bed_prefix))

G_full <- BEDMatrix(bed_prefix, simple_names = TRUE)   # n × p, memory-mapped
n      <- nrow(G_full)
p      <- ncol(G_full)
cat(sprintf("  %d individuals × %d SNPs\n", n, p))

# ── Select causal SNPs ────────────────────────────────────────────────────────

if (opt$`effect-dist` == "gaussian") {
  # Infinitesimal: every SNP contributes
  n_causal   <- p
  causal_idx <- seq_len(p)
} else {
  n_causal   <- max(1L, round(opt$`p-causal` * p))
  causal_idx <- sort(sample.int(p, n_causal))
}
cat(sprintf("  Causal SNPs: %d (%.4f of %d)\n", n_causal, n_causal / p, p))

# ── Two-pass computation over causal SNP columns ──────────────────────────────
# Pass 1: column means and SDs (needed to standardise genotypes)
# Pass 2: accumulate G_std * beta into the genetic-value vector g
#
# Chunked to avoid loading all causal genotype columns into RAM at once.
# Chunk of 5 000 SNPs × 15 000 individuals × 4 bytes ≈ 300 MB.

CHUNK <- 5000L

cat("  Pass 1: computing column means and SDs ...\n")
col_means <- numeric(n_causal)
col_sds   <- numeric(n_causal)

for (k in seq(1L, n_causal, by = CHUNK)) {
  k2     <- min(k + CHUNK - 1L, n_causal)
  gchunk <- G_full[, causal_idx[k:k2], drop = FALSE]
  col_means[k:k2] <- colMeans(gchunk, na.rm = TRUE)
  col_sds[k:k2]   <- apply(gchunk, 2L, sd, na.rm = TRUE)
}
col_sds[col_sds == 0] <- 1   # constant column (monomorphic after rounding): zero effect

# ── Draw effect sizes ─────────────────────────────────────────────────────────

# Approximate prior scale so E[Var(G*beta)] ≈ h2.
# Under independence: Var(G_std * beta) = sum(beta_j^2) ≈ n_causal * sigma_b^2
# → sigma_b = sqrt(h2 / n_causal).
sigma_b <- sqrt(opt$h2 / n_causal)

betas <- switch(
  opt$`effect-dist`,
  gaussian  = rnorm(n_causal, mean = 0, sd = sigma_b),
  spikeslab = rnorm(n_causal, mean = 0, sd = sigma_b),
  scaledt   = {
    df <- 4
    # Standardise so that Var(beta_j) = sigma_b^2 (requires df > 2)
    sigma_b * rt(n_causal, df = df) / sqrt(df / (df - 2))
  }
)

# ── Pass 2: compute genetic values ────────────────────────────────────────────

cat("  Pass 2: computing genetic values ...\n")
g <- numeric(n)

for (k in seq(1L, n_causal, by = CHUNK)) {
  k2     <- min(k + CHUNK - 1L, n_causal)
  gchunk <- G_full[, causal_idx[k:k2], drop = FALSE]

  # Impute missing genotypes with column mean (rare after MAF filtering)
  for (j in seq_len(k2 - k + 1L)) {
    mis <- is.na(gchunk[, j])
    if (any(mis)) gchunk[mis, j] <- col_means[k + j - 1L]
  }

  # Standardise and accumulate contribution
  gstd <- sweep(
    sweep(gchunk, 2L, col_means[k:k2], "-"),
    2L, col_sds[k:k2], "/"
  )
  g <- g + as.numeric(gstd %*% betas[k:k2])
}

# ── Scale residual variance to hit target h2 exactly ─────────────────────────

var_g <- var(g)
if (var_g < .Machine$double.eps * 100) {
  stop("Genetic variance is effectively zero. Check genotype file and --p-causal.")
}
sigma2_e  <- var_g * (1 - opt$h2) / opt$h2
epsilon   <- rnorm(n, mean = 0, sd = sqrt(sigma2_e))
liability <- g + epsilon
h2_realised <- var_g / (var_g + sigma2_e)   # equals opt$h2 by construction

# ── Generate phenotype ────────────────────────────────────────────────────────

if (opt$`trait-type` == "quantitative") {
  Y        <- liability
  n_cases  <- NA_integer_
  prev_obs <- NA_real_
} else {
  # Liability-threshold model: liability on N(0,1) scale after rescaling
  lia_scaled <- (liability - mean(liability)) / sd(liability)
  threshold  <- qnorm(1 - opt$prevalence)
  Y          <- as.integer(lia_scaled > threshold)
  n_cases    <- sum(Y)
  prev_obs   <- mean(Y)
}

# ── Write outputs ─────────────────────────────────────────────────────────────

out_dir <- dirname(opt$out)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Phenotype file
fam <- read.table(
  paste0(bed_prefix, ".fam"), header = FALSE,
  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
)
pheno_df <- data.frame(FID = fam$FID, IID = fam$IID, Y = Y)
write.table(pheno_df,
            file      = paste0(opt$out, ".pheno"),
            sep       = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote     = FALSE)

# 2. Causal SNP effects
bim <- read.table(
  paste0(bed_prefix, ".bim"), header = FALSE,
  col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2")
)
causal_df <- data.frame(
  SNP      = bim$SNP[causal_idx],
  CHR      = bim$CHR[causal_idx],
  BP       = bim$BP[causal_idx],
  A1       = bim$A1[causal_idx],
  freq_A1  = col_means / 2,          # approximate allele frequency
  beta_std = betas                   # effect on standardised genotype scale
)
write.table(causal_df,
            file      = paste0(opt$out, ".causal_snps.tsv"),
            sep       = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote     = FALSE)

# 3. Simulation parameter summary
params_df <- data.frame(
  n             = n,
  p             = p,
  n_causal      = n_causal,
  p_causal      = n_causal / p,
  effect_dist   = opt$`effect-dist`,
  trait_type    = opt$`trait-type`,
  h2_target     = opt$h2,
  h2_realised   = h2_realised,
  prevalence_target  = if (opt$`trait-type` == "binary") opt$prevalence else NA_real_,
  prevalence_realised= if (opt$`trait-type` == "binary") prev_obs       else NA_real_,
  n_cases       = if (opt$`trait-type` == "binary") n_cases             else NA_integer_,
  seed          = opt$seed
)
write.table(params_df,
            file      = paste0(opt$out, ".sim_params.tsv"),
            sep       = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote     = FALSE)

# ── Console summary ───────────────────────────────────────────────────────────

cat(sprintf(
  "[simulate_phenotypes] Done.\n  trait=%s  effect_dist=%s\n  h2 target=%.3f  realised=%.3f\n  n_causal=%d  p_causal=%.4f\n",
  opt$`trait-type`, opt$`effect-dist`,
  opt$h2, h2_realised,
  n_causal, n_causal / p
))
if (opt$`trait-type` == "binary") {
  cat(sprintf("  prevalence target=%.3f  realised=%.3f  n_cases=%d / %d\n",
              opt$prevalence, prev_obs, n_cases, n))
}
