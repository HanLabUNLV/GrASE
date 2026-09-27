#!/usr/bin/env Rscript
## Add length-normalized pi to an exontest results table.
##
## The test's pi is a RAW COUNT ratio -- correct for the beta-binomial
## likelihood, which is defined on counts -- but not a molecular proportion when
## the distinct set and reference differ in length. This appends the per-base
## reading alongside. It changes no call: p-values, padj, pi and delta_pi are
## untouched.
##
## Junction-substituted sides are left NA: a junction is a point feature with no
## length, so a per-base pi is undefined for it.
##
## Usage:
##   Rscript scripts/annotate_pi.R --results=<annotated.txt> \
##     --gff_dir=<per-gene dexseq gff dir> --outfile=<out.txt>

suppressMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-r", "--results"), type = "character", default = NULL),
  make_option(c("-g", "--gff_dir"), type = "character", default = NULL),
  make_option(c("-o", "--outfile"), type = "character", default = NULL))))
for (r in c("results", "gff_dir", "outfile"))
  if (is.null(opt[[r]])) stop("--", r, " must be specified", call. = FALSE)

rdir <- file.path(dirname(dirname(normalizePath(
          sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "R")
if (!dir.exists(rdir)) rdir <- "R"
for (f in list.files(rdir, pattern = "\\.R$", full.names = TRUE))
  try(source(f), silent = TRUE)
stopifnot(exists("annotate_pi_perbase"))

res <- utils::read.delim(opt$results, stringsAsFactors = FALSE)
message("rows: ", nrow(res))
res <- annotate_pi_perbase(res, opt$gff_dir)
utils::write.table(res, opt$outfile, sep = "\t", quote = FALSE, row.names = FALSE)

ok <- !is.na(res$delta_pi_perbase)
message(sprintf("per-base pi computed for %d of %d rows (%.1f%%)",
                sum(ok), nrow(res), 100 * mean(ok)))
if (any(ok)) {
  message(sprintf("median |delta_pi| raw %.4f  ->  per-base %.4f",
                  stats::median(abs(res$delta_pi[ok])),
                  stats::median(abs(res$delta_pi_perbase[ok]))))
  for (thr in c(0.05, 0.1, 0.2)) {
    a <- sum(abs(res$delta_pi[ok]) >= thr)
    b <- sum(abs(res$delta_pi_perbase[ok]) >= thr)
    message(sprintf("  |delta_pi| >= %.2f : raw %6d   per-base %6d   (%+.1f%%)",
                    thr, a, b, 100 * (b - a) / max(a, 1)))
  }
}
message("wrote ", opt$outfile)
