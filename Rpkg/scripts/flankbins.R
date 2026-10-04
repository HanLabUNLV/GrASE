#!/usr/bin/env Rscript
## Emit flanking bins for the TSS/TTS localization test.
##
## One bin per bipartition SIDE: a fixed-width window immediately beyond that
## side's free end. The bin is then counted like any other feature and tested
## against the same reference as the side's distinct set, under the same call
## condition. A real boundary means the condition effect stops at the edge --
## the distinct set is called and the flank is not.
##
## The window is taken REGARDLESS of what lies beyond the boundary -- intron,
## another isoform's exon, or intergenic space. It is not clipped at the next
## exonic part and there is no minimum-width gate. That is deliberate: the test
## is a between-condition contrast against a shared reference, so anything in
## the flank that is not changing contributes equally to both conditions and
## cancels. Restricting to "empty" flanks is what the superseded step statistic
## did, and it is why that statistic was confounded -- see below.
##
## SUPERSEDES the absolute-step reading. `step = adj/(adj + D)` asked whether
## coverage falls to ZERO beyond the boundary, which is only the right question
## when nothing else is transcribed through the flank. With other isoforms
## expressed a genuine boundary gives a step CHANGE, not a drop to zero, so that
## statistic lost power in proportion to how little of the local coverage the
## terminating isoform contributed: on the DICE activation panel its call rate
## rose monotonically with pi (20% at pi<0.2 to 40% at pi>0.8) and its confirmed
## set was 7.4x enriched for each gene's OUTERMOST boundary. The calibration
## bins this script used to emit existed only to fit a threshold for that
## statistic, and are gone with it.
##
## Count the GFF with the SAME settings that produced the exonic part counts --
## for DICE that is DICE/scripts/dexseq_count_dice.sh, which records and verifies
## them. Do NOT fold these bins into the DEXSeq flattening: adding features
## changes read assignment for existing exonic parts, so every current count
## shifts. A separate annotation counted separately cannot disturb them.
##
## Geometry comes from the SPLICE GRAPH, not the flattened GFF. R/L edges are
## the distinct TSS/TTS positions and ex_part edges are the exonic parts
## (coordinates from their endpoint vertices). Transcript membership is not
## needed -- an R edge points at its route's first node, so anything 5' of it is
## by definition off that path.
##
## EVENT IDS. `event` is the ROW INDEX of the per-gene split file, which is what
## exoncnt assigns (`bipartitions$event = rownames(bipartitions)`,
## R/exoncnt_functions.R:141). It is not a column in the split tables --
## filterTSSTTS.R does not write one -- and an earlier version of this script
## required it to be, which meant the guard failed on every file and the script
## silently emitted nothing. Rows are numbered before any filtering so the ids
## line up with the exonic arm.
##
## Usage:
##   Rscript scripts/flankbins.R --split_dir=<filtered TSS/TTS splits> \
##     --graphml_dir=<per-gene graphml> --outdir=<dir> [--width=150] \
##     [--chr_gff=<aggregate dexseq gff>]

suppressMessages({library(optparse)})

option_list <- list(
  make_option(c("-i", "--split_dir"), type = "character", default = NULL,
              help = "directory of filtered bipartition split tables"),
  make_option(c("-g", "--graphml_dir"), type = "character", default = NULL,
              help = "directory of per-gene <gene>.graphml files"),
  make_option(c("-c", "--chr_gff"), type = "character", default = NULL,
              help = "aggregate dexseq gff; needed ONLY for graphs built before chrom was stored on the graph"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "output directory"),
  make_option(c("-w", "--width"), type = "integer", default = 150L,
              help = "flank bin width [default %default]"),
  ## Restricting to ONE arm is not cosmetic: `event` is the row index WITHIN a
  ## per-gene split file, so bins must be numbered off the same file the exonic
  ## arm was counted from. exoncnt.R maps -a TSS to alt "TSSTTS"
  ## (scripts/exoncnt.R:62-63) and globs "bipartition.<alt>.txt$".
  make_option(c("-p", "--pattern"), type = "character",
              default = "\\.bipartition\\.TSSTTS\\.txt$",
              help = "split-file pattern [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (r in c("split_dir", "graphml_dir", "outdir"))
  if (is.null(opt[[r]])) stop("--", r, " must be specified", call. = FALSE)

rdir <- file.path(dirname(dirname(normalizePath(
          sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "R")
if (!dir.exists(rdir)) rdir <- "R"
for (f in list.files(rdir, pattern = "\\.R$", full.names = TRUE))
  try(source(f), silent = TRUE)
suppressMessages(library(igraph))
stopifnot(exists("graph_exonic_parts"), exists("graph_side_boundary"),
          exists("flank_bin"), exists("free_end_direction"))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## Newer graphs carry chrom, strand and gene as graph attributes, so nothing is
## read from the annotation at all. --chr_gff is a fallback for graphs built
## before chrom was recorded (the DICE v34 graphs are such a set).
chr_map <- if (!is.null(opt$chr_gff) && file.exists(opt$chr_gff)) {
  m <- parse_gff_chr_map(opt$chr_gff)
  message("chromosome fallback map: ", length(m), " genes")
  m
} else list()

gcache <- new.env(parent = emptyenv())
read_gene <- function(gene) {
  if (!is.null(gcache[[gene]])) return(gcache[[gene]])
  fn <- file.path(opt$graphml_dir, paste0(gene, ".graphml"))
  v <- NULL
  if (file.exists(fn)) {
    g <- tryCatch(igraph::read_graph(fn, format = "graphml"), error = function(e) NULL)
    if (!is.null(g)) {
      ch <- graph_chrom(g)
      if (is.na(ch)) ch <- chr_map[[gene]]
      v <- list(g = g, strand = graph_strand(g),
                chrom = if (is.null(ch)) NA_character_ else ch)
    }
  }
  gcache[[gene]] <- if (is.null(v)) NA else v
  gcache[[gene]]
}

split_files <- list.files(opt$split_dir, pattern = opt$pattern, full.names = TRUE)
message("split tables: ", length(split_files))

rows <- list(); k <- 0L
n_nograph <- 0L; n_noboundary <- 0L; n_notRL <- 0L; n_empty <- 0L
for (sf in split_files) {
  d <- try(utils::read.delim(sf, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) next
  if (!all(c("gene", "source", "sink") %in% names(d))) next
  for (i in seq_len(nrow(d))) {
    ## event == row index of THIS file, matching exoncnt's numbering
    event <- i
    kind <- if (identical(as.character(d$source[i]), "R")) "TSS" else
            if (identical(as.character(d$sink[i]), "L"))   "TTS" else NA
    if (is.na(kind)) { n_notRL <- n_notRL + 1L; next }
    gg <- read_gene(d$gene[i])
    if (length(gg) == 1L && is.na(gg)) { n_nograph <- n_nograph + 1L; next }
    strand <- gg$strand; chrom <- gg$chrom
    if (is.na(strand) || !strand %in% c("+", "-") || is.na(chrom)) {
      n_nograph <- n_nograph + 1L; next
    }
    dir <- free_end_direction(kind, strand)
    for (side in c(1L, 2L)) {
      sd <- if (side == 1L) d$setdiff1[i] else d$setdiff2[i]
      ## a side with no exonic distinct set has no D to compare the flank
      ## against, so it is out of scope for this test
      if (is.na(sd) || sd %in% c("", "NA")) { n_empty <- n_empty + 1L; next }
      pf <- if (side == 1L) d$path1[i] else d$path2[i]
      b <- graph_side_boundary(gg$g, pf, kind, strand)
      if (is.na(b)) { n_noboundary <- n_noboundary + 1L; next }
      ## avail = width and min_width = 1: never clipped, never unscorable
      fb <- flank_bin(b, dir, opt$width, opt$width, 1L)
      k <- k + 1L
      rows[[k]] <- data.frame(
        bin_id = sprintf("F%d_%d", event, side),
        gene = d$gene[i], event = event, side = side, kind = kind,
        chrom = chrom, strand = strand, boundary = b,
        start = fb$start, end = fb$end, width = fb$width,
        distinct_parts = sd, stringsAsFactors = FALSE)
    }
  }
}
if (!k) stop("no flank bins emitted -- check --split_dir and --graphml_dir", call. = FALSE)
man <- do.call(rbind, rows)
message("flank bins: ", nrow(man), "  over ", length(unique(man$gene)), " genes")
message("  skipped: ", n_notRL, " not R/L anchored, ", n_nograph, " no graph/strand/chrom, ",
        n_empty, " empty distinct set, ", n_noboundary, " no route terminus")

## ---- DEXSeq-format GFF ------------------------------------------------------
## Only `exonic_part` lines are read by dexseq_count.py, which names each
## feature gene_id + ":" + exonic_part_number (dexseq_count.py:96-97) by plain
## string concatenation -- no numeric assumption. So the part number carries the
## event and side directly: the count rowname `<gene>:F<event>_<side>` maps back
## to the bipartition path with no lookup table. `aggregate_gene` lines are
## ignored by the counter but written for format completeness.
##
## Bins are grouped under their SOURCE GENE so that two overlapping bins of one
## gene are two features of the same gene rather than an `_ambiguous` read.
ord <- order(man$chrom, man$start, man$end)
mo <- man[ord, , drop = FALSE]
src <- "flankbins.R"
gff <- character(0)
for (g in unique(mo$gene)) {
  s <- mo[mo$gene == g, , drop = FALSE]
  gff <- c(gff, paste(s$chrom[1], src, "aggregate_gene", min(s$start), max(s$end),
                      ".", s$strand[1], ".",
                      sprintf('gene_id "%s"', g), sep = "\t"))
  gff <- c(gff, paste(s$chrom, src, "exonic_part", s$start, s$end, ".", s$strand, ".",
                      sprintf('gene_id "%s"; transcripts "NA"; exonic_part_number "%s"',
                              g, s$bin_id), sep = "\t"))
}
gff_path <- file.path(opt$outdir, "flankbins.dexseq.gff")
writeLines(gff, gff_path)

man_path <- file.path(opt$outdir, "flankbins_manifest.tsv")
utils::write.table(man, man_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("wrote ", gff_path)
message("wrote ", man_path)
message("next: count it with the SAME settings as the exonic parts, e.g.")
message("  bash DICE/scripts/dexseq_count_dice.sh ", gff_path, " <bam_list> <outdir>")
