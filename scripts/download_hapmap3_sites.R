#!/usr/bin/env Rscript
# download_hapmap3_sites.R
#
# Downloads the HapMap3 SNP list used by the LDPred2 / bigsnpr tutorial
# and writes it as a tab-delimited file suitable for --hapmap3-sites.
#
# Output columns: chrom  pos  rsid  a0  a1
#   chrom – integer chromosome number (no 'chr' prefix)
#   pos   – 1-based GRCh37 position
#   rsid  – rs identifier
#   a0    – reference allele
#   a1    – effect/alternate allele
#
# Usage:
#   Rscript scripts/download_hapmap3_sites.R \
#       --out resources/hapmap3_sites.tsv
#
# Requires: optparse
# The HapMap3 RDS is fetched from the bigsnpr figshare page on first run and
# cached locally so subsequent calls do not re-download it.

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
  make_option("--out", type = "character",
              default = "resources/hapmap3_sites.tsv",
              help = "Output path for the sites TSV [default: %default]"),
  make_option("--cache-dir", type = "character",
              default = "resources",
              help = "Directory in which to cache the downloaded RDS [default: %default]"),
  make_option("--chroms", type = "character", default = NULL,
              help = "Comma-separated list of chromosomes to keep, e.g. '21,22'. Default: keep all autosomes (1-22).")
)

opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$`cache-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$out),  recursive = TRUE, showWarnings = FALSE)

# ── Download (or load from cache) ────────────────────────────────────────────
# The HapMap3 variant info RDS is published by Privé et al. alongside the
# LDPred2 paper.  Search "LDpred2 HapMap3 figshare" for the current URL;
# the file is typically named something like "map_hm3_ldpred2.rds".

cache_file <- file.path(opt$`cache-dir`, "map_hm3_ldpred2.rds")

if (file.exists(cache_file)) {
  cat(sprintf("[download_hapmap3_sites] Using cached file: %s\n", cache_file))
  info <- readRDS(cache_file)
} else {
  # map_hm3_ldpred2.rds published by Privé et al. alongside the LDPred2 paper.
  # Figshare item 13034123, file 25503788.
  HM3_URL <- "https://ndownloader.figshare.com/files/25503788"
  cat(sprintf("[download_hapmap3_sites] Downloading HapMap3 map from:\n  %s\n", HM3_URL))
  tmp <- tempfile(fileext = ".rds")
  tryCatch(
    download.file(HM3_URL, destfile = tmp, mode = "wb", quiet = FALSE),
    error = function(e) stop(
      "Download failed: ", conditionMessage(e), "\n",
      "Download the file manually from the URL above and save it to:\n  ",
      cache_file
    )
  )
  info <- readRDS(tmp)
  file.copy(tmp, cache_file, overwrite = TRUE)
  unlink(tmp)
  cat(sprintf("[download_hapmap3_sites] Cached to: %s\n", cache_file))
}

cat(sprintf("[download_hapmap3_sites] HapMap3 table: %d SNPs\n", nrow(info)))
cat("Columns: ", paste(names(info), collapse = ", "), "\n")

# ── Normalise column names ────────────────────────────────────────────────────
# The RDS may use different column names across bigsnpr versions.
# Common variants: chr/chrom, pos/bp, a0/ref, a1/alt, rsid/id
col_map <- list(
  chrom = c("chr", "chrom", "CHR", "chromosome"),
  pos   = c("pos", "bp", "BP", "POS", "position"),
  rsid  = c("rsid", "id", "ID", "snp", "SNP"),
  a0    = c("a0", "ref", "REF", "A2"),
  a1    = c("a1", "alt", "ALT", "A1")
)

find_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L)
    stop("Cannot find column (tried: ", paste(candidates, collapse = ", "), ")")
  hit[1L]
}

info2 <- data.frame(
  chrom = as.integer(info[[ find_col(info, col_map$chrom) ]]),
  pos   = as.integer(info[[ find_col(info, col_map$pos)   ]]),
  rsid  = info[[ find_col(info, col_map$rsid) ]],
  a0    = info[[ find_col(info, col_map$a0)   ]],
  a1    = info[[ find_col(info, col_map$a1)   ]],
  stringsAsFactors = FALSE
)

# ── Optionally subset chromosomes ─────────────────────────────────────────────
if (!is.null(opt$chroms)) {
  keep_chroms <- as.integer(strsplit(opt$chroms, ",")[[1]])
  info2 <- info2[info2$chrom %in% keep_chroms, ]
  cat(sprintf("[download_hapmap3_sites] Kept %d SNPs on chrom(s): %s\n",
              nrow(info2), paste(keep_chroms, collapse = ", ")))
}

info2 <- info2[order(info2$chrom, info2$pos), ]

# ── Write ─────────────────────────────────────────────────────────────────────
write.table(info2,
            file      = opt$out,
            sep       = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote     = FALSE)

cat(sprintf("[download_hapmap3_sites] Wrote %d SNPs to: %s\n",
            nrow(info2), opt$out))
