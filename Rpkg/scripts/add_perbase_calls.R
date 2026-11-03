#!/usr/bin/env Rscript
## Add length-normalized pi to existing exontest.R results and recompute `significant`, without refitting.
## For runs made without --gff_dir: annotate_pi_perbase() adds len_D, len_S, pi_ctrl_perbase, pi_trt_perbase
## and delta_pi_perbase, then add_significant() redoes the call with --use_perbase (default both: raw and
## per-base delta_pi must both pass, raw alone where a side has no length) and --min_effect, and
## dpi_raw_pass / dpi_perbase_pass are written as exontest.R writes them. p.value, padj and every other column are kept.
## A pi_ref column from before the rename to pi_ctrl is renamed on the way in.
## Each file is rewritten in place; the original is kept as <file>.bak_rawpi (an existing backup is not
## overwritten).
## Usage: Rscript add_perbase_calls.R --results <a.annotated.txt>,<b.annotated.txt> --gff_dir <dir> --read_length <r>
##          [--padj_threshold 0.01] [--delta 0] [--min_dpi 0.05] [--min_dpi_sj 0.05] [--min_effect 0.5]
##          [--min_reads 10]
## Use the padj threshold, delta and read floor of the exontest.R run that wrote the files.
suppressPackageStartupMessages({library(optparse); library(data.table); library(grase)})
o <- parse_args(OptionParser(option_list = list(
  make_option("--results",        type = "character", help = "annotated exontest.R results, comma separated"),
  make_option("--gff_dir",        type = "character", help = "per-gene DEXSeq GFFs (<gene>.dexseq.gff or <gene>/<gene>.dexseq.gff)"),
  make_option("--read_length",    type = "integer", help = "read length of the counted reads (mate length if paired)"),
  make_option("--padj_threshold", type = "double", default = 0.01),
  make_option("--delta",          type = "double", default = 0),
  make_option("--min_dpi",        type = "double", default = 0.05),
  make_option("--min_dpi_sj",     type = "double", default = 0.05),
  make_option("--min_effect",     type = "double", default = 0.5),
  make_option("--min_reads",      type = "double", default = 10),
  make_option("--use_perbase",    type = "character", default = "both",
              help = "both, TRUE (per-base) or FALSE (raw) [%default]"))))
for (r in c("results", "gff_dir", "read_length")) if (is.null(o[[r]])) stop("--", r, " is required")
o$use_perbase <- switch(o$use_perbase, both = "both", "TRUE" = TRUE, "FALSE" = FALSE,
                        stop("--use_perbase must be both, TRUE or FALSE"))

for (f in strsplit(o$results, ",")[[1]]) {
  d <- fread(f, sep = "\t", na.strings = "NA")
  if ("pi_ref" %in% names(d)) setnames(d, "pi_ref", "pi_ctrl")
  ## pi_ref_perbase is the pre-rename name of pi_ctrl_perbase: drop it too, or a file written
  ## before the rename keeps a stale copy beside the recomputed column.
  drop <- intersect(c("len_D", "len_S", "pi_ctrl_perbase", "pi_trt_perbase", "delta_pi_perbase",
                      "pi_ref_perbase", "dpi_raw_pass", "dpi_perbase_pass"), names(d))
  if (length(drop)) d[, (drop) := NULL]
  old <- d$significant %in% c(TRUE, "TRUE")
  d <- as.data.table(annotate_pi_perbase(as.data.frame(d), path.expand(o$gff_dir), o$read_length))
  d <- as.data.table(add_significant(as.data.frame(d), o$padj_threshold, o$delta, o$min_dpi, o$min_dpi_sj,
                                     o$min_reads, use_perbase = o$use_perbase, min_effect = o$min_effect))
  ## as exontest.R: the cutoff per row is min_dpi, or min_dpi_sj where the tested distinct set is empty
  sd <- fifelse(grepl("diff1", d$comparison), d$setdiff1, d$setdiff2)
  thr <- fifelse(is.na(sd) | sd %in% c("NA", ""), o$min_dpi_sj, o$min_dpi)
  d[, dpi_raw_pass := abs(delta_pi) >= thr]
  d[, dpi_perbase_pass := abs(delta_pi_perbase) >= thr]
  cat(sprintf("%s\n  per-base pi for %d of %d rows; significant %d -> %d (gained %d, lost %d)\n", f,
              sum(!is.na(d$delta_pi_perbase)), nrow(d), sum(old), sum(d$significant),
              sum(d$significant & !old), sum(old & !d$significant)))
  bak <- paste0(f, ".bak_rawpi")
  if (!file.exists(bak)) file.copy(f, bak)
  fwrite(d, f, sep = "\t", na = "NA", quote = FALSE)
}
