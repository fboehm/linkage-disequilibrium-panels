#!/usr/bin/env Rscript
# score_pgs_variance.R
#
# Compute per-individual posterior variance of the PGS using the LD-aware
# formula:
#
#   Var(PGS_i) = w_i^T R w_i,   w_i = sqrt(BETA_VAR) * G_i  (element-wise)
#
# where R is the pre-computed LD correlation matrix (SFBM from
# prepare_ldpred2_ref.R), BETA_VAR_j is the posterior variance of the effect
# at SNP j, and G_ij is the allele dosage (0/1/2) for individual i at SNP j.
#
# When R = I this reduces to the independence approximation
#   Var(PGS_i) = Σ_j BETA_VAR_j * G_ij²
#
# SNPs are matched between the betas file and the SFBM by CHR + BP, so the
# formula is robust to different SNP-ID conventions across methods.
#
# If no SNPs can be matched to the SFBM (or all BETA_VAR are NA/zero), the
# output contains NA PGS_VAR for every test individual and a warning is emitted.
#
# Usage:
#   Rscript scripts/score_pgs_variance.R \
#       --bed      results/plink/<sim_method>/gwas/rep<r>/merged.bed \
#       --betas    results/pgs_weights/<method>/.../betas.tsv \
#       --test-ids results/splits/<sim_method>/rep<r>/test.txt \
#       --ld-sfbm  results/ldpred2_work/<sim_method>/rep<r>/<ancestry>/n<N>/ld_sfbm.rds \
#       --ld-snps  results/ldpred2_work/<sim_method>/rep<r>/<ancestry>/n<N>/matched_snps.tsv \
#       --out      results/pgs/<method>/.../scores_var.tsv

suppressPackageStartupMessages({
  library(bigsnpr)
  library(bigsparser)
  library(data.table)
  library(optparse)
  library(parallel)
})

opt_list <- list(
  make_option("--bed",      type = "character",
              help = "Path to merged PLINK .bed file"),
  make_option("--betas",    type = "character",
              help = "betas.tsv with SNP, A1, BETA, BETA_VAR, CHR, BP columns"),
  make_option("--test-ids", type = "character",
              help = "Test-set keep-file (FID IID, no header)"),
  make_option("--ld-sfbm",  type = "character",
              help = "Path to ld_sfbm.rds (SFBM saved by prepare_ldpred2_ref.R)"),
  make_option("--ld-snps",  type = "character",
              help = "Path to matched_snps.tsv (SNP table saved by prepare_ldpred2_ref.R)"),
  make_option("--out",      type = "character",
              help = "Output TSV: FID IID PGS_VAR"),
  make_option("--ncores",   type = "integer", default = 1L,
              help = "Number of cores for parallel variance computation [default: %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))

missing_args <- Filter(is.null, list(
  bed       = opt$bed,
  betas     = opt$betas,
  `test-ids`= opt$`test-ids`,
  `ld-sfbm` = opt$`ld-sfbm`,
  `ld-snps` = opt$`ld-snps`,
  out       = opt$out
))
if (length(missing_args) > 0)
  stop("Required arguments missing: ", paste(names(missing_args), collapse = ", "))

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)

# ── Helper: write NA output and exit cleanly ─────────────────────────────────
write_na_output <- function(fid, iid, path, msg) {
  cat(sprintf("[score_pgs_variance] WARNING: %s; writing NA variance.\n", msg))
  fwrite(data.table(FID = fid, IID = iid, PGS_VAR = NA_real_), path, sep = "\t")
  cat(sprintf("[score_pgs_variance] Wrote %d rows (NA) to: %s\n",
              length(fid), path))
  quit(save = "no", status = 0L)
}

# ── Load GWAS BED into FBM ───────────────────────────────────────────────────
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
map[, chr := as.integer(sub("^chr", "", chr))]

fam <- setDT(obj$fam)
setnames(fam, 1:2, c("FID", "IID"))

cat(sprintf("[score_pgs_variance] BED: %d individuals × %d SNPs\n",
            nrow(G), ncol(G)))

# ── Identify test individuals ────────────────────────────────────────────────
test_ids <- fread(opt$`test-ids`, header = FALSE, col.names = c("FID", "IID"))
ind_row  <- which(paste(fam$FID, fam$IID) %in%
                    paste(test_ids$FID, test_ids$IID))
if (length(ind_row) == 0L)
  stop("No test individuals found in BED file.")
cat(sprintf("[score_pgs_variance] %d test individuals\n", length(ind_row)))

# ── Load LD reference ────────────────────────────────────────────────────────
cat(sprintf("[score_pgs_variance] Loading SFBM: %s\n", opt$`ld-sfbm`))
corr  <- readRDS(opt$`ld-sfbm`)
m_ld  <- nrow(corr)
cat(sprintf("[score_pgs_variance] SFBM: %d × %d\n", m_ld, m_ld))

msnps <- fread(opt$`ld-snps`)
rename_map <- c("_NUM_ID_" = "NUM_ID", "_FLIP_" = "FLIP")
old_nm <- intersect(names(rename_map), names(msnps))
if (length(old_nm) > 0) setnames(msnps, old_nm, rename_map[old_nm])

if (nrow(msnps) != m_ld)
  stop(sprintf("matched_snps has %d rows but SFBM has %d; files may be mismatched.",
               nrow(msnps), m_ld))
cat(sprintf("[score_pgs_variance] Matched-SNP table: %d SNPs\n", nrow(msnps)))

# ── Load per-SNP posterior variance ─────────────────────────────────────────
betas <- fread(opt$betas)
if (!"BETA_VAR" %in% names(betas))
  stop("betas.tsv does not contain a BETA_VAR column.")
if (!all(c("CHR", "BP") %in% names(betas)))
  stop("betas.tsv must contain CHR and BP columns for SNP matching.")

betas[, CHR := as.integer(CHR)]
betas[, BP  := as.integer(BP)]
betas[, BETA_VAR := as.numeric(BETA_VAR)]

cat(sprintf("[score_pgs_variance] Betas: %d SNPs\n", nrow(betas)))

# ── Build BETA_VAR vector in SFBM SNP order ──────────────────────────────────
# Match betas to matched_snps by CHR + BP; unmatched SFBM SNPs get BETA_VAR=0
# (zero posterior uncertainty = no contribution to the LD-aware sum).
msnps[, chr := as.integer(chr)]
msnps[, pos := as.integer(pos)]

beta_var_sfbm <- numeric(m_ld)   # defaults to 0

m_idx <- match(
  paste(msnps$chr, msnps$pos),
  paste(betas$CHR,  betas$BP)
)
in_betas       <- !is.na(m_idx)
bv_matched     <- betas$BETA_VAR[m_idx[in_betas]]
has_valid_bv   <- in_betas & !is.na(c(bv_matched, numeric(sum(!in_betas)))[seq_len(m_ld)])
# safer indexing:
valid_rows <- which(in_betas)
for (i in seq_along(valid_rows)) {
  bv <- betas$BETA_VAR[m_idx[valid_rows[i]]]
  if (!is.na(bv))
    beta_var_sfbm[valid_rows[i]] <- bv
}

n_bv_nonzero <- sum(beta_var_sfbm > 0)
cat(sprintf("[score_pgs_variance] %d / %d SFBM SNPs have non-zero BETA_VAR\n",
            n_bv_nonzero, m_ld))

if (n_bv_nonzero == 0L)
  write_na_output(fam$FID[ind_row], fam$IID[ind_row], opt$out,
                  "no SFBM SNPs have non-zero BETA_VAR")

# ── Match SFBM SNPs to GWAS BED columns ─────────────────────────────────────
# Try rsid first; fall back to chr+pos.
ind_col_sfbm <- match(msnps$rsid, map$rsid)

unmatched <- is.na(ind_col_sfbm)
if (any(unmatched)) {
  map_key   <- paste(map$chr, map$pos)
  msnps_key <- paste(msnps$chr, msnps$pos)
  ind_col_sfbm[unmatched] <- match(msnps_key[unmatched], map_key)
}

n_in_bed <- sum(!is.na(ind_col_sfbm))
cat(sprintf("[score_pgs_variance] %d / %d SFBM SNPs found in BED\n",
            n_in_bed, m_ld))

if (n_in_bed == 0L)
  write_na_output(fam$FID[ind_row], fam$IID[ind_row], opt$out,
                  "no SFBM SNPs found in BED file")

# ── Extract test-set genotypes for all SFBM SNPs ────────────────────────────
# G_sub[i, j] = dosage for test individual i at SFBM SNP j (0 where not in BED)
in_bed  <- !is.na(ind_col_sfbm)
G_sub   <- matrix(0.0, nrow = length(ind_row), ncol = m_ld)
G_sub[, in_bed] <- as.matrix(G[ind_row, ind_col_sfbm[in_bed]])

cat(sprintf("[score_pgs_variance] Genotype matrix: %d × %d\n",
            nrow(G_sub), ncol(G_sub)))

# ── Compute LD-aware variance: Var(PGS_i) = w_i^T R w_i ─────────────────────
# w_i = sqrt(BETA_VAR) * G_sub[i, ]  (element-wise scaling of each row)
v     <- sqrt(beta_var_sfbm)          # length m_ld
W     <- t(t(G_sub) * v)             # scale columns: W[i,j] = v[j] * G_sub[i,j]

cat(sprintf(
  "[score_pgs_variance] Computing LD-aware PGS variance (%d cores) ...\n",
  opt$ncores
))

pgs_var <- unlist(mclapply(seq_len(nrow(W)), function(k) {
  w_k   <- W[k, ]
  rw_k  <- bigsparser::sp_prodVec(corr, w_k)
  sum(w_k * rw_k)
}, mc.cores = opt$ncores))

# ── Write output ─────────────────────────────────────────────────────────────
out_df <- data.table(
  FID     = fam$FID[ind_row],
  IID     = fam$IID[ind_row],
  PGS_VAR = as.vector(pgs_var)
)
fwrite(out_df, opt$out, sep = "\t")
cat(sprintf("[score_pgs_variance] Wrote %d rows to: %s\n",
            nrow(out_df), opt$out))
