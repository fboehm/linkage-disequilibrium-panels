#!/usr/bin/env Rscript
# sanity_check_sample_beta_scale.R
#
# Confirms on which scale `snp_ldpred2_auto`'s `sample_beta` matrix is stored.
#
#   Hypothesis A (assumed in run_ldpred2.R):
#     sample_beta is on the *standardised* scale (beta * sd).
#     => rowMeans(sample_beta) ≈ beta_est * sd
#
#   Hypothesis B:
#     sample_beta is on the *unstandardised* scale (same as beta_est).
#     => rowMeans(sample_beta) ≈ beta_est
#
# With report_step = 1 and a single chain, rowMeans(sample_beta) should equal
# the within-chain posterior mean by construction, so whichever equality holds
# tells us the scale unambiguously.

suppressPackageStartupMessages({
  library(bigsnpr)
  library(Matrix)
})

set.seed(1)

# ── Toy data from bigsnpr's bundled example ──────────────────────────────────
bigsnp  <- snp_attachExtdata()
G       <- bigsnp$genotypes
CHR     <- bigsnp$map$chromosome
POS     <- bigsnp$map$physical.pos
y       <- bigsnp$fam$affection - 1   # 0/1 phenotype as continuous
N       <- length(y)

# Simple GWAS (linear) to get marginal betas + SEs
gwas    <- big_univLinReg(G, y)
df_beta <- data.frame(
  beta    = gwas$estim,
  beta_se = gwas$std.err,
  n_eff   = N
)

# Drop any SNPs with zero SE (monomorphic) to keep snp_ldpred2_auto happy
keep    <- df_beta$beta_se > 0 & is.finite(df_beta$beta) & is.finite(df_beta$beta_se)
df_beta <- df_beta[keep, ]
ind_kep <- which(keep)

# ── Build a small in-memory LD matrix on chromosome 1 only ───────────────────
chr1    <- which(CHR[ind_kep] == 1L)
ind_col <- ind_kep[chr1]
df_beta <- df_beta[chr1, ]

corr0   <- snp_cor(G, ind.col = ind_col, size = 500, ncores = 1)
corr    <- bigsparser::as_SFBM(corr0)

cat(sprintf("[sanity] SNPs in test: %d\n", nrow(df_beta)))

# ── Run LDpred2-auto with report_step = 1 on one chain ───────────────────────
# Note: beta_est is a Rao-Blackwellised cumulative average, while sample_beta
# stores raw per-iteration draws. Both estimate the same posterior mean, so
# rowMeans(sample_beta) ≈ beta_est (on the standardised scale) up to MC noise.
# We use a long run to push MC noise down and then check the slope.
BURN_IN  <- 500L
NUM_ITER <- 5000L

multi_auto <- snp_ldpred2_auto(
  corr        = corr,
  df_beta     = df_beta,
  h2_init     = 0.1,
  vec_p_init  = 0.01,
  burn_in     = BURN_IN,
  num_iter    = NUM_ITER,
  report_step = 1L,
  ncores      = 1
)
auto <- multi_auto[[1]]

stopifnot(!is.null(auto$sample_beta))
stopifnot(ncol(auto$sample_beta) == NUM_ITER)

# ── Compare scales ───────────────────────────────────────────────────────────
# beta_est is divided by sd inside snp_ldpred2_auto (see its source),
# with sd = 1 / sqrt(n_eff * beta_se^2 + beta^2).
sd_snp <- 1 / sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)

sample_mean_std   <- as.numeric(Matrix::rowMeans(auto$sample_beta))
expected_if_std   <- auto$beta_est * sd_snp   # Hypothesis A
expected_if_unstd <- auto$beta_est             # Hypothesis B

# Restrict comparison to SNPs with non-trivial posterior mass
nz <- which(abs(auto$beta_est) > 1e-8)

# Slope of rowMeans(sample_beta) on each candidate expectation.
# Under the correct hypothesis, slope should be ≈ 1 with tight residuals.
fit_slope <- function(y, x) {
  ok <- is.finite(x) & is.finite(y)
  coef(lm(y[ok] ~ 0 + x[ok]))[[1]]
}

slope_A <- fit_slope(sample_mean_std[nz], expected_if_std[nz])
slope_B <- fit_slope(sample_mean_std[nz], expected_if_unstd[nz])

rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))
corr_coef <- function(a, b) suppressWarnings(cor(a, b, use = "complete.obs"))

cat("\n[sanity] Comparing rowMeans(sample_beta) against beta_est on both scales\n")
cat(sprintf("         (using %d post-burn-in MCMC draws on %d non-zero SNPs)\n\n",
            NUM_ITER, length(nz)))
cat(sprintf("  Hypothesis A (standardised):   slope = %.5f   cor = %.6f   RMSE = %.3e\n",
            slope_A,
            corr_coef(sample_mean_std[nz], expected_if_std[nz]),
            rmse(sample_mean_std[nz], expected_if_std[nz])))
cat(sprintf("  Hypothesis B (unstandardised): slope = %.5f   cor = %.6f   RMSE = %.3e\n",
            slope_B,
            corr_coef(sample_mean_std[nz], expected_if_unstd[nz]),
            rmse(sample_mean_std[nz], expected_if_unstd[nz])))

# Verdict by slope-to-1 closeness
slope_tol <- 0.02
dev_A <- abs(slope_A - 1)
dev_B <- abs(slope_B - 1)

if (dev_A < slope_tol && dev_A < dev_B) {
  cat(sprintf("\n[sanity] PASS: sample_beta is on the STANDARDISED scale (slope %.4f ≈ 1).\n",
              slope_A))
  cat("         run_ldpred2.R is correct to divide var_std by sd^2.\n")
} else if (dev_B < slope_tol && dev_B < dev_A) {
  cat(sprintf("\n[sanity] FAIL: sample_beta is on the UNSTANDARDISED scale (slope %.4f ≈ 1).\n",
              slope_B))
  cat("         run_ldpred2.R should NOT divide var_std by sd^2.\n")
  quit(status = 1)
} else {
  cat("\n[sanity] INCONCLUSIVE: neither slope is within tolerance of 1. Inspect:\n")
  print(head(data.frame(
    sample_mean   = sample_mean_std[nz],
    beta_est      = auto$beta_est[nz],
    beta_est_x_sd = expected_if_std[nz],
    sd            = sd_snp[nz]
  )))
  quit(status = 2)
}
