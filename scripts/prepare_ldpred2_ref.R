#!/usr/bin/env Rscript
# prepare_ldpred2_ref.R
#
# Precompute the shared LD reference objects used by all run_ldpred2.R jobs
# for a given replicate:
#   1. Load panel genotypes into a bigsnpr FBM backing file.
#   2. Match a reference set of GWAS summary statistics to the panel to
#      determine the shared matched-SNP set (same for every phenotype run).
#   3. Compute the sparse LD correlation matrix (SFBM) with snp_cor.
#   4. Save panel.rds / panel.bk, ld_sfbm.sbk, and matched_snps.tsv.
#
# Outputs are placed in --out-dir and reused by all 144 run_ldpred2 jobs,
# reducing the total LD computation from 144× to 1×.
#
# Usage:
#   Rscript scripts/prepare_ldpred2_ref.R \
#       --panel-bed    results/plink/panel/rep1/merged.bed \
#       --ref-sumstats results/gwas/rep1/quantitative/h2_0.1/pc_0.001/gaussian/sumstats.tsv \
#       --out-dir      results/ldpred2_work/rep1/shared \
#       --window       3000000 \
#       --ncores       4

suppressPackageStartupMessages({
  library(bigsnpr)
  library(bigsparser)
  library(data.table)
  library(optparse)
})

opt_list <- list(
  make_option("--panel-bed",    type = "character",
              help = "Path to panel PLINK .bed file"),
  make_option("--ref-sumstats", type = "character",
              help = "Reference linear-regression sumstats (BETA/SE columns)"),
  make_option("--out-dir",      type = "character",
              help = "Directory for all output files"),
  make_option("--window",       type = "double",   default = 3e6,
              help = "LD window size in bp [default: 3e6]"),
  make_option("--ncores",       type = "integer",  default = 1L,
              help = "Number of cores for snp_cor [default: 1]"),
  make_option("--n-panel",      type = "integer",  default = NULL,
              help = "Number of panel individuals to subsample (NULL = use all)"),
  make_option("--seed",         type = "integer",  default = 42L,
              help = "Random seed for individual subsampling [default: 42]"),
  make_option("--hm3-positions", type = "character", default = NULL,
              help = "TSV with HM3 SNP positions: chrom pos (and optionally rsid a0 a1). When provided, the matched SNP set is restricted to HM3 positions before computing the LD matrix.")
)

opt <- parse_args(OptionParser(option_list = opt_list))

missing_args <- Filter(is.null, list(
  `panel-bed`    = opt$`panel-bed`,
  `ref-sumstats` = opt$`ref-sumstats`,
  `out-dir`      = opt$`out-dir`
))
if (length(missing_args) > 0)
  stop("Required arguments missing: ", paste(names(missing_args), collapse = ", "))

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
# Use a local temp directory for the FBM backing file; OneDrive/network paths
# on WSL2 do not support the memory-mapped I/O that bigsnpr requires.
bk_path   <- file.path(tempdir(), "ldpred2_panel")
sfbm_path <- file.path(opt$`out-dir`, "ld_sfbm")

# ── 1. Load panel into FBM ────────────────────────────────────────────────────

rds_file <- paste0(bk_path, ".rds")
if (!file.exists(rds_file)) {
  cat("[prepare_ldpred2_ref] Reading panel BED into FBM ...\n", flush = TRUE)
  snp_readBed(opt$`panel-bed`, backingfile = bk_path)
}
obj      <- snp_attach(rds_file)
G_panel  <- obj$genotypes
map_panel <- setDT(obj$map)
setnames(map_panel,
         c("chromosome", "marker.ID", "genetic.dist",
           "physical.pos", "allele1", "allele2"),
         c("chr", "rsid", "cm", "pos", "a1", "a0"))
map_panel[, chr := as.integer(sub("^chr", "", chr))]

cat(sprintf("[prepare_ldpred2_ref] Panel: %d individuals × %d SNPs\n",
            nrow(G_panel), ncol(G_panel)), flush = TRUE)

# ── Optional individual subsampling ───────────────────────────────────────────
# When --n-panel is set and smaller than the full panel, randomly subsample
# that many individuals.

n_avail <- nrow(G_panel)
if (!is.null(opt$`n-panel`) && as.integer(opt$`n-panel`) < n_avail) {
  set.seed(opt$seed)
  ind_row <- sort(sample(n_avail, as.integer(opt$`n-panel`)))
  cat(sprintf("[prepare_ldpred2_ref] Subsampling %d of %d panel individuals (seed=%d)\n",
              length(ind_row), n_avail, opt$seed), flush = TRUE)
} else {
  ind_row <- seq_len(n_avail)
}

# When the subsample is smaller than the full panel, copy only those rows into
# a new FBM before running snp_cor.  snp_cor reads every row of the FBM even
# when ind.row is small, so keeping the full 10K-individual matrix in the
# backing file wastes memory proportional to n_avail, not n_panel.
if (length(ind_row) < n_avail) {
  bk_sub  <- paste0(bk_path, "_sub")
  G_panel <- bigstatsr::big_copy(G_panel, ind.row = ind_row,
                                 backingfile = bk_sub)
  ind_row <- bigstatsr::rows_along(G_panel)
  cat(sprintf("[prepare_ldpred2_ref] Copied subsetted FBM: %d × %d\n",
              nrow(G_panel), ncol(G_panel)), flush = TRUE)
}

# ── 2. Match reference sumstats to panel ──────────────────────────────────────

cat(sprintf("[prepare_ldpred2_ref] Reading reference sumstats: %s\n",
            opt$`ref-sumstats`), flush = TRUE)

ss <- fread(opt$`ref-sumstats`, data.table = TRUE)
setnames(ss, 1L, "CHROM")
ss <- ss[TEST == "ADD"]
if ("BETA" %in% names(ss)) {
  ss <- ss[!is.na(BETA) & !is.na(SE) & SE > 0]
} else if ("OR" %in% names(ss)) {
  setnames(ss, "LOG(OR)_SE", "SE")
  ss[, BETA := log(OR)]
  ss <- ss[!is.na(BETA) & !is.na(SE) & SE > 0 & is.finite(BETA)]
} else {
  stop("Reference sumstats must have BETA or OR column.")
}
ss[, CHROM := as.integer(sub("^chr", "", CHROM))]

cat(sprintf("[prepare_ldpred2_ref] Reference sumstats: %d SNPs\n", nrow(ss)),
    flush = TRUE)

map_ref <- data.frame(
  chr  = map_panel$chr,
  pos  = map_panel$pos,
  a0   = map_panel$a0,
  a1   = map_panel$a1,
  rsid = map_panel$rsid,
  stringsAsFactors = FALSE
)

df_dummy <- data.frame(
  chr     = ss$CHROM,
  pos     = as.integer(ss$POS),
  a0      = ss$REF,
  a1      = ss$A1,
  beta    = ss$BETA,
  beta_se = ss$SE,
  n_eff   = 1L,
  rsid    = ss$ID,
  stringsAsFactors = FALSE
)

info_ref <- tryCatch(
  snp_match(df_dummy, map_ref, join_by_pos = TRUE),
  error = function(e) {
    message("[prepare_ldpred2_ref] snp_match error: ", conditionMessage(e))
    NULL
  }
)

if (is.null(info_ref) || nrow(info_ref) < 50L)
  stop(sprintf("[prepare_ldpred2_ref] Only %d SNPs matched.",
               if (is.null(info_ref)) 0L else nrow(info_ref)))

cat(sprintf("[prepare_ldpred2_ref] Matched %d SNPs\n", nrow(info_ref)),
    flush = TRUE)

# ── 2b. Restrict to HM3 SNPs (optional) ──────────────────────────────────────
if (!is.null(opt$`hm3-positions`)) {
  hm3 <- fread(opt$`hm3-positions`, data.table = TRUE)
  hm3[, chrom := as.integer(sub("^chr", "", chrom))]
  hm3[, pos   := as.integer(pos)]
  hm3_key <- hm3[, .(chrom, pos)]
  info_ref <- as.data.table(info_ref)
  info_ref <- info_ref[hm3_key, on = c("chr" = "chrom", "pos" = "pos"), nomatch = NULL]
  cat(sprintf("[prepare_ldpred2_ref] After HM3 restriction: %d SNPs\n", nrow(info_ref)),
      flush = TRUE)
  if (nrow(info_ref) < 50L)
    stop("[prepare_ldpred2_ref] Fewer than 50 SNPs after HM3 restriction. ",
         "Check that --hm3-positions covers the correct genome build and chromosomes.")
}

# ── 3. Compute sparse LD matrix (SFBM) ────────────────────────────────────────

cat(sprintf("[prepare_ldpred2_ref] Computing LD matrix (window = %.0f bp, ncores = %d) ...\n",
            opt$window, opt$ncores), flush = TRUE)

chromosomes <- sort(unique(info_ref$chr))
corr        <- NULL

for (chr_i in chromosomes) {
  ind_chr   <- which(info_ref$chr == chr_i)
  panel_idx <- info_ref[["_NUM_ID_"]][ind_chr]

  cat(sprintf("  chr%d: %d SNPs\n", chr_i, length(ind_chr)), flush = TRUE)

  corr_chr <- snp_cor(
    Gna       = G_panel,
    ind.row   = ind_row,
    ind.col   = panel_idx,
    infos.pos = info_ref$pos[ind_chr],
    size      = opt$window,
    ncores    = opt$ncores
  )

  if (is.null(corr)) {
    corr <- as_SFBM(corr_chr, sfbm_path, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

cat(sprintf("[prepare_ldpred2_ref] LD matrix: %d × %d (sparse)\n",
            nrow(corr), ncol(corr)), flush = TRUE)

# Save the SFBM R6 object so run_ldpred2.R can reload it
corr$save()   # writes ld_sfbm.rds alongside ld_sfbm.sbk

# ── 4. Save matched-SNP table ─────────────────────────────────────────────────

# Add SFBM row index (1-based, same order as the SFBM)
info_ref[["sfbm_row"]] <- seq_len(nrow(info_ref))

# Save subset of columns needed for matching in run_ldpred2.R
cols_to_save <- c("chr", "pos", "a0", "a1", "rsid",
                  "_NUM_ID_", "_FLIP_", "sfbm_row")
cols_to_save <- intersect(cols_to_save, names(info_ref))
if (!is.data.table(info_ref)) setDT(info_ref)
fwrite(info_ref[, cols_to_save, with = FALSE],
       file.path(opt$`out-dir`, "matched_snps.tsv"), sep = "\t")

cat(sprintf("[prepare_ldpred2_ref] Saved matched_snps.tsv (%d rows)\n",
            nrow(info_ref)), flush = TRUE)
cat("[prepare_ldpred2_ref] Done.\n", flush = TRUE)
