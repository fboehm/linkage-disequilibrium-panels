#!/usr/bin/env Rscript
# split_subjects.R
#
# Randomly split GWAS individuals into training (discovery) and test
# (evaluation) sets.  Outputs two plain-text PLINK keep-files.
#
# Usage:
#   Rscript scripts/split_subjects.R \
#       --fam        results/plink/gwas/rep1/merged.fam \
#       --train-frac 0.8 \
#       --seed       1001 \
#       --train      results/splits/rep1/train.txt \
#       --test       results/splits/rep1/test.txt

suppressPackageStartupMessages(library(optparse))

opt_list <- list(
  make_option("--fam",        type = "character",
              help = "Path to PLINK .fam file"),
  make_option("--train-frac", type = "double",    default = 0.8,
              help = "Proportion of individuals for training [default: %default]"),
  make_option("--seed",       type = "integer",
              help = "Random seed"),
  make_option("--train",      type = "character",
              help = "Output path for training-set keep-file (FID IID, no header)"),
  make_option("--test",       type = "character",
              help = "Output path for test-set keep-file (FID IID, no header)")
)

opt <- parse_args(OptionParser(option_list = opt_list))

missing_args <- Filter(is.null, list(fam   = opt$fam,
                                     seed  = opt$seed,
                                     train = opt$train,
                                     test  = opt$test))
if (length(missing_args) > 0)
  stop("Required arguments missing: ", paste(names(missing_args), collapse = ", "))

fam <- read.table(opt$fam, header = FALSE,
                  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
n <- nrow(fam)

set.seed(opt$seed)
train_idx <- sort(sample.int(n, size = round(opt$`train-frac` * n)))
test_idx  <- setdiff(seq_len(n), train_idx)

for (path in c(opt$train, opt$test))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

write.table(fam[train_idx, c("FID", "IID")],
            file = opt$train, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(fam[test_idx, c("FID", "IID")],
            file = opt$test,  sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat(sprintf("[split_subjects] %d train + %d test = %d total  (seed = %d)\n",
            length(train_idx), length(test_idx), n, opt$seed))
