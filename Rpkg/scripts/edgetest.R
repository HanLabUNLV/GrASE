#!/usr/bin/env Rscript
## Annotate exontest results with the TSS/TTS edge-coverage diagnostic.
##
## Runs AFTER exontest.R. Adds columns; changes no call. A call's p-value and
## delta_pi are untouched -- what this adds is whether the annotated terminal
## boundary is one where coverage actually stops.
##
##   step = adj_n / (adj_n + d_n)    0 = coverage stops dead
##                                   0.5 = continuation
##
## READ THE FLAG ASYMMETRICALLY. Thresholds are fitted per dataset from the
## calibration bins, and the reported misread rates say how far to trust each
## label. On the DICE activation panel negatives were misread 4.1% but
## positives 38.9%, because many single-promoter genes carry real upstream
## transcription. "boundary" is specific; "continuation" is evidence, not
## proof.
##
## Usage:
##   Rscript scripts/edgetest.R --results=<exontest annotated.txt> \
##     --manifest=<flankbins_manifest.tsv> --flankcounts=<featureCounts out> \
##     --countdir=<per-side counts dir> --outfile=<annotated.edge.txt>

suppressMessages({library(optparse)})

option_list <- list(
  make_option(c("-r", "--results"), type = "character", default = NULL,
              help = "exontest annotated results"),
  make_option(c("-m", "--manifest"), type = "character", default = NULL,
              help = "flankbins_manifest.tsv from flankbins.R"),
  make_option(c("-f", "--flankcounts"), type = "character", default = NULL,
              help = "featureCounts output over flankbins.saf"),
  make_option(c("-c", "--countdir"), type = "character", default = NULL,
              help = "per-side count directory (for D and S counts)"),
  make_option(c("-o", "--outfile"), type = "character", default = NULL,
              help = "output path"),
  make_option(c("--min_reads"), type = "double", default = 2,
              help = "read floor for D and the flank [%default]"),
  make_option(c("--read_length"), type = "double", default = 100,
              help = "bases of coverage one read contributes [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (r in c("results", "manifest", "flankcounts", "outfile"))
  if (is.null(opt[[r]])) stop("--", r, " must be specified", call. = FALSE)

rdir <- file.path(dirname(dirname(normalizePath(
          sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "R")
if (!dir.exists(rdir)) rdir <- "R"
for (f in list.files(rdir, pattern = "\\.R$", full.names = TRUE))
  try(source(f), silent = TRUE)
stopifnot(exists("edge_step"), exists("calibrate_boundary"))

man <- utils::read.delim(opt$manifest, stringsAsFactors = FALSE)
fc  <- utils::read.delim(opt$flankcounts, comment.char = "#",
                         stringsAsFactors = FALSE)
## featureCounts: Geneid Chr Start End Strand Length <sample columns>
scol <- setdiff(names(fc), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
message("flank bins counted: ", nrow(fc), "   samples: ", length(scol))
adj <- stats::setNames(rowMeans(fc[, scol, drop = FALSE], na.rm = TRUE), fc$Geneid)

man$adj <- adj[man$bin_id]

## --- calibration -------------------------------------------------------
## flankbins.R emits each control as a PAIR sharing a pair_id: role "ADJ" is
## the window outside the boundary, role "D" the matched interior window it is
## compared against. Positives should give step near 0 (a real boundary),
## negatives near 0.5 (both windows inside one continuous exon).
stopifnot("role" %in% names(man), "pair_id" %in% names(man))
cal_pair <- function(rt) {
  m <- man[man$row_type == rt, , drop = FALSE]
  if (!nrow(m)) return(numeric(0))
  a <- stats::setNames(m$adj[m$role == "ADJ"],     m$pair_id[m$role == "ADJ"])
  d <- stats::setNames(m$adj[m$role == "D"],       m$pair_id[m$role == "D"])
  la <- stats::setNames(m$len_ADJ[m$role == "ADJ"], m$pair_id[m$role == "ADJ"])
  ld <- stats::setNames(m$len_D[m$role == "D"],     m$pair_id[m$role == "D"])
  ids <- intersect(names(a), names(d))
  vapply(ids, function(i) edge_step(a[[i]], d[[i]], NA, la[[i]], ld[[i]])$step,
         numeric(1))
}
cal <- calibrate_boundary(cal_pair("calib_pos"), cal_pair("calib_neg"))
message(sprintf("calibration: threshold %.3f   positives misread %.1f%%   negatives misread %.1f%%   (n=%d/%d)",
                cal$threshold, 100 * cal$pos_misread, 100 * cal$neg_misread,
                cal$n_pos, cal$n_neg))

## --- score the real sides ---------------------------------------------
res <- utils::read.delim(opt$results, stringsAsFactors = FALSE)
res$.side <- ifelse(grepl("^diff1", res$comparison), 1L, 2L)
key <- paste(res$gene, res$event, res$.side)
## sides carry a single ADJ bin; D comes from the results table
side_bins <- man[man$row_type == "side" & man$role == "ADJ", , drop = FALSE]
mkey <- paste(side_bins$gene, side_bins$event, side_bins$side)
idx <- match(key, mkey)
man_s <- side_bins

res$edge_bin      <- man_s$bin_id[idx]
res$edge_tier     <- man_s$tier[idx]
res$edge_flank_bp <- man_s$flank_width[idx]
## D comes from the side's distinct-set count. NOTE this is only on the same
## scale as the flanking bin for EXON-sourced sides: both are interval read
## counts. For a junction-substituted side diff*_mean is a spliced-read count,
## which is not comparable to interval coverage and would push step toward
## "continuation" spuriously. flankbins.R already skips sides with an empty
## exonic distinct set so they carry no bin, but label them explicitly rather
## than letting them fall through as a coverage problem.
d_raw <- if (!is.null(res$diff1_mean) && !is.null(res$diff2_mean))
           ifelse(res$.side == 1L, res$diff1_mean, res$diff2_mean) else NA_real_
sd_col <- ifelse(res$.side == 1L, res$setdiff1, res$setdiff2)
no_exonic_D <- is.na(sd_col) | sd_col %in% c("", "NA")
st <- mapply(function(a, d, la, ld) {
        if (is.na(a) || is.na(d)) return(c(NA_real_, NA_real_))
        r <- edge_step(a, d, NA, la, ld, 100)
        c(r$step, r$d_n)
      }, man_s$adj[idx], d_raw, man_s$len_ADJ[idx], man_s$len_D[idx])
res$edge_step <- st[1, ]
res$edge_d_n  <- st[2, ]
res$edge_adj_n <- man_s$adj[idx] * 100 / man_s$len_ADJ[idx]
res$boundary_support <- boundary_support(res$edge_step, res$edge_d_n, cal,
                                         res$edge_adj_n, opt$min_reads,
                                         opt$read_length)
## the diagnostic does not apply without an exonic distinct set: there is no
## alternative terminal exon whose free edge could be probed
res$boundary_support[no_exonic_D] <- "not_applicable"
res$.side <- NULL

utils::write.table(res, opt$outfile, sep = "\t", quote = FALSE, row.names = FALSE)
tb <- table(res$boundary_support, useNA = "ifany")
message("boundary_support:")
for (n in names(tb)) message(sprintf("   %-14s %6d", n, tb[[n]]))
message("wrote ", opt$outfile)
