#!/usr/bin/env Rscript
## Re-derive the `significant` flag using the LENGTH-NORMALIZED delta_pi.
##
## pi is a raw count ratio, so a single absolute |delta_pi| threshold is not
## scale-fair: the same threshold demands a different real effect depending on
## how the distinct set and reference compare in length. This regates on the
## per-base value wherever it is defined and on the raw value elsewhere.
##
## Junction-substituted sides keep the raw gate by necessity -- a junction is a
## point feature with no length, so there is no per-base pi for it. min_dpi_sj
## already carries a partial scale correction for those.
##
## Nothing else is recomputed: p-values and padj are untouched, so this is a
## re-thresholding, not a re-test.
##
## Usage:
##   Rscript scripts/resignify_perbase.R --results=<annotated.txt> \
##     --gff_dir=<per-gene dexseq gff> --outfile=<out.txt> \
##     [--padj_thr=0.01] [--delta=0] [--min_dpi=0.1] [--min_dpi_sj=0.05] \
##     [--min_reads=10]

suppressMessages({library(optparse); library(grase)})
opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-r", "--results"),  type = "character", default = NULL),
  make_option(c("-g", "--gff_dir"),  type = "character", default = NULL),
  make_option(c("-o", "--outfile"),  type = "character", default = NULL),
  make_option("--padj_thr",  type = "double", default = 0.01),
  make_option("--delta",     type = "double", default = 0),
  make_option("--min_dpi",   type = "double", default = 0.1),
  make_option("--min_dpi_sj", type = "double", default = 0.05),
  make_option("--min_reads", type = "double", default = 10))))
for (r in c("results", "gff_dir", "outfile"))
  if (is.null(opt[[r]])) stop("--", r, " must be specified", call. = FALSE)

res <- utils::read.delim(opt$results, stringsAsFactors = FALSE)
message("rows: ", nrow(res))

if (!"delta_pi_perbase" %in% names(res)) {
  message("adding per-base pi from ", opt$gff_dir, " ...")
  res <- annotate_pi_perbase(res, opt$gff_dir)
}
ok <- !is.na(res$delta_pi_perbase)
message(sprintf("per-base pi available for %d of %d rows (%.1f%%)",
                sum(ok), nrow(res), 100 * mean(ok)))

old <- add_significant(res, opt$padj_thr, opt$delta, opt$min_dpi,
                       opt$min_dpi_sj, opt$min_reads)$significant
new <- add_significant(res, opt$padj_thr, opt$delta, opt$min_dpi,
                       opt$min_dpi_sj, opt$min_reads,
                       use_perbase = TRUE)$significant
res$significant_raw_gate <- old
res$significant <- new
utils::write.table(res, opt$outfile, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("\nsignificant: raw gate %d  ->  per-base gate %d  (%+d, %+.1f%%)",
                sum(old, na.rm = TRUE), sum(new, na.rm = TRUE),
                sum(new, na.rm = TRUE) - sum(old, na.rm = TRUE),
                100 * (sum(new, na.rm = TRUE) - sum(old, na.rm = TRUE)) /
                  max(sum(old, na.rm = TRUE), 1)))
message(sprintf("  gained (per-base clears, raw did not): %d",
                sum(new & !old, na.rm = TRUE)))
message(sprintf("  lost   (raw cleared, per-base does not): %d",
                sum(old & !new, na.rm = TRUE)))
message("wrote ", opt$outfile)
