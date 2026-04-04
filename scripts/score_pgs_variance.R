#!/usr/bin/env Rscript
# score_pgs_variance.R
#
# Compute per-individual posterior variance of the PGS:
#
#   Var(PGS_i) = sum_j  BETA_VAR_j * G_ij^2
#
# where BETA_VAR_j is the posterior variance of the effect size at SNP j
# (variance across LDpred2-auto MCMC chains, or across PRS-CS posterior
# samples), and G_ij is the allele dosage (0/1/2) for individual i at SNP j.
#
# This uses the independence approximation (posterior betas treated as
# uncorrelated across SNPs), which is consistent with the modelling
# assumptions of both LDpred2-auto and PRS-CS.
#
# Note: BETA_VAR is invariant to allele-flip, so no strand harmonisation
# is needed.  G_ij^2 is computed on the raw dosage coding in the BED file.
#
# Usage:
#   Rscript scripts/score_pgs_variance.R \
#       --bed      results/plink/<sim_method>/gwas/rep<r>/merged.bed \
#       --betas    results/pgs_weights/<method>/.../betas.tsv \
#       --test-ids results/splits/<sim_method>/rep<r>/test.txt \
#       --out      results/pgs/<method>/.../scores_var.tsv \
#       --ncores   1

suppressPackageStartupMessages({
  library(bigsnpr)
  library(data.table)
  library(optparse)
})

opt_list <- list(
  make_option("--bed",      type = "character",
              help = "Path to merged PLINK .bed file"),
  make_option("--betas",    type = "character",
              help = "betas.tsv with SNP, A1, BETA, BETA_VAR columns"),
  make_option("--test-ids", type = "character",
              help = "Test-set keep-file (FID IID, no header)"),
  make_option("--ncores",   type = "integer",  default = 1L,
              help = "Cores for big_apply [default: 1]"),
  make_option("--out",      type = "character",
              help = "Output TSV: FID IID PGS_VAR")
)

opt <- parse_args(OptionParser(option_list = opt_list))

missing_args <- Filter(is.null, list(
  bed       = opt$bed,
  betas     = opt$betas,
  `test-ids`= opt$`test-ids`,
  out       = opt$out
))
if (length(missing_args) > 0)
  stop("Required arguments missing: ", paste(names(missing_args), collapse = ", "))

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)

# ── Load panel BED into FBM ──────────────────────────────────────────────────
bk_path  <- file.path(tempdir(), "pgs_var_panel")
rds_file <- paste0(bk_path, ".rds")
if (!file.exists(rds_file))
  snp_readBed(opt$bed, backingfile = bk_path)

obj <- snp_attach(rds_file)
G   <- obj$genotypes
map <- setDT(obj$map)
setnames(map,
         c("chromosome", "marker.ID", "genetic.dist",
           "physical.pos", "allele1", "allele2"),
         c("chr", "rsid", "cm", "pos", "a1", "a0"))

fam <- setDT(obj$fam)
# plink .fam columns: family.ID sample.ID paternal maternal sex affection
setnames(fam, 1:2, c("FID", "IID"))

cat(sprintf("[score_pgs_variance] Panel: %d individuals × %d SNPs\n",
            nrow(G), ncol(G)))

# ── Identify test individuals ────────────────────────────────────────────────
test_ids <- fread(opt$`test-ids`, header = FALSE, col.names = c("FID", "IID"))
ind_row  <- which(paste(fam$FID, fam$IID) %in%
                    paste(test_ids$FID, test_ids$IID))
if (length(ind_row) == 0L)
  stop("No test individuals found in BED file.")
cat(sprintf("[score_pgs_variance] %d test individuals\n", length(ind_row)))

# ── Load per-SNP posterior variance ─────────────────────────────────────────
betas <- fread(opt$betas)
if (!"BETA_VAR" %in% names(betas))
  stop("betas.tsv does not contain a BETA_VAR column.")

n_na <- sum(is.na(betas$BETA_VAR))
if (n_na > 0) {
  cat(sprintf("[score_pgs_variance] Dropping %d SNPs with NA BETA_VAR\n", n_na))
  betas <- betas[!is.na(BETA_VAR)]
}
cat(sprintf("[score_pgs_variance] %d SNPs with valid BETA_VAR\n", nrow(betas)))

if (nrow(betas) == 0L) {
  cat("[score_pgs_variance] No SNPs with valid BETA_VAR; writing NA variance.\n")
  out_df <- data.table(
    FID     = fam$FID[ind_row],
    IID     = fam$IID[ind_row],
    PGS_VAR = NA_real_
  )
  fwrite(out_df, opt$out, sep = "\t")
  cat(sprintf("[score_pgs_variance] Wrote %d rows to: %s\n", nrow(out_df), opt$out))
  quit(save = "no", status = 0L)
}

# ── Match SNPs to BED map ────────────────────────────────────────────────────
ind_col <- match(betas$SNP, map$rsid)
n_miss  <- sum(is.na(ind_col))
if (n_miss > 0) {
  cat(sprintf("[score_pgs_variance] Dropping %d SNPs not found in BED\n", n_miss))
  betas   <- betas[!is.na(ind_col)]
  ind_col <- ind_col[!is.na(ind_col)]
}
cat(sprintf("[score_pgs_variance] %d SNPs matched to BED\n", length(ind_col)))

if (length(ind_col) == 0L) {
  cat("[score_pgs_variance] No SNPs matched to BED; writing NA variance.\n")
  out_df <- data.table(
    FID     = fam$FID[ind_row],
    IID     = fam$IID[ind_row],
    PGS_VAR = NA_real_
  )
  fwrite(out_df, opt$out, sep = "\t")
  cat(sprintf("[score_pgs_variance] Wrote %d rows to: %s\n", nrow(out_df), opt$out))
  quit(save = "no", status = 0L)
}

# Build a full-length vector so big_apply can index by FBM column index
beta_var_by_col           <- rep(0.0, ncol(G))
beta_var_by_col[ind_col]  <- betas$BETA_VAR

# ── Compute Var(PGS_i) = sum_j  BETA_VAR_j * G_ij^2 ────────────────────────
cat("[score_pgs_variance] Computing per-individual PGS variance ...\n")

pgs_var <- big_apply(
  X       = G,
  a.FUN   = function(X, ind, ind_row, bv_all) {
    g_blk <- X[ind_row, ind]          # n_test × block_size  (integer 0/1/2)
    bv    <- bv_all[ind]              # block_size
    as.matrix(rowSums(sweep(g_blk^2, 2L, bv, `*`)))
  },
  a.combine = `+`,
  ind       = ind_col,
  ind_row   = ind_row,
  bv_all    = beta_var_by_col,
  ncores    = opt$ncores
)

# ── Write output ─────────────────────────────────────────────────────────────
out_df <- data.table(
  FID     = fam$FID[ind_row],
  IID     = fam$IID[ind_row],
  PGS_VAR = as.vector(pgs_var)
)
fwrite(out_df, opt$out, sep = "\t")
cat(sprintf("[score_pgs_variance] Wrote %d rows to: %s\n",
            nrow(out_df), opt$out))
