#!/usr/bin/env Rscript
## Emit INTRONIC PARTS as quantification units, the intronic analogue of DEXSeq
## exonic parts.
##
## An intronic part is a maximal region of a gene that is exonic in NO isoform:
## the gene span minus every one of its exonic parts. Because the gene's exonic
## parts are disjoint (verified: 0 overlapping pairs in 519,701 adjacent pairs
## over 18,920 genes), exonic and intronic parts together form a COMPLETE,
## DISJOINT TILING of the gene -- every base belongs to exactly one unit. That
## property is the point of this construction. It is what makes a walk outward
## from a bipartition terminus well defined at every step, with no gaps to fall
## through and no overlapping units to choose between.
##
## SUPERSEDES the fixed-width flanking bin (flankbins.R). That bin was a
## synthetic 150bp window placed beyond each side's free end, which produced
## three pathologies this construction does not have:
##   - it could land on sequence belonging to no annotated unit, giving verdicts
##     (`flank_empty`, `no_flank_bin`, `flank_untested`) that described the
##     window rather than the biology;
##   - its width was unrelated to unit boundaries, so "is the neighbour within
##     the window" became a free parameter -- and agreement with an independent
##     neighbour test collapsed to chance (kappa 0.06) beyond that width;
##   - it was keyed per SIDE, so the same interval was emitted and counted once
##     per side sharing a terminus: 35,011 bins for 22,034 distinct intervals,
##     37.1% redundant.
## Keying on the gene's annotation instead removes the redundancy by
## construction -- a part is a part, and every side that abuts it looks it up.
##
## QUANTIFICATION IS UNIVERSAL, TESTING IS NOT. Every intronic part of every
## gene is counted, so the count table is a reusable asset: a later panel or
## contrast draws test rows from it without recounting. Restricting the
## annotation to the current call set, as the flank arm did, forced a recount
## whenever the call set moved. Test rows are emitted downstream only for units
## abutting a significant partition, and they inherit that partition's shared
## reference so the contrast stays comparable to D's.
##
## NOT IN THE GRAPH. No transcript traverses an intronic part, so it has no
## legitimate edge, and adding one would change bubble enumeration and therefore
## every existing call. Note the graph's `in` edges are a DIFFERENT object: they
## are per-splice-event and transcript-specific, so one `in` edge can span
## straight across an exonic part of another isoform (PDCD4 `in 10,13`,
## 110,876,754-110,881,233, crosses E008 entirely). An `in` edge's span
## decomposes into alternating exonic and intronic parts; that decomposition is
## where the graph and this tiling meet.
##
## BINNING. Intronic parts span 37 bp to 702 kb (median 1,416, p90 10,557), so a
## single part-level mean is right for most and meaningless for the tail -- 4.8%
## of parts hold 45.5% of intronic sequence and average 44 kb. Parts up to T are
## one bin. Longer parts get a bin of width W at EACH END and the interior is not
## quantified at all: the two end bins already give the contrast (both elevated
## and roughly equal means flat pervasive transcription, one elevated means a
## decaying tail), and skipping the interior drops the counted sequence from
## 1,344 Mb to 152 Mb.
##
## T >= 2W IS A HARD CONSTRAINT, not a preference: 30% of intronic parts are
## shorter than 600 bp, so with W = 300 and no threshold the two end bins would
## overlap and double-count the same reads.
##
## Coordinates are written with scientific notation explicitly disabled. R's
## default `scipen` emits a round coordinate as "2.3e+07", which is exact but
## crashes a strict integer parse downstream -- flankbins.R has this bug in its
## manifest `boundary` field (1 row) and geom_parts.tsv has it in start/end (2
## rows).
##
## Count the GFF with the SAME settings that produced the exonic part counts.
## Do NOT fold these units into the DEXSeq flattening: adding features changes
## read assignment for existing exonic parts, so every current count shifts. A
## separate annotation counted separately cannot disturb them. Ambiguity is not
## a concern here in any case -- an intronic part is exonic in no isoform, so
## nothing can be assigned `_ambiguous` against it.
##
## Usage:
##   Rscript scripts/intronparts.R --gff_dir=<per-gene dexseq gff dir> \
##     --outdir=<dir> [--T=2000] [--W=300] [--genes=<file, one id per line>]

suppressMessages({library(optparse); library(grase)})

option_list <- list(
  make_option(c("-g", "--gff_dir"), type = "character", default = NULL,
              help = "directory of per-gene <gene>.dexseq.gff files"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "output directory"),
  make_option("--T", type = "integer", default = 2000L,
              help = "parts up to this length are one bin [default %default]"),
  make_option("--W", type = "integer", default = 300L,
              help = "end-bin width for parts longer than T [default %default]"),
  make_option("--genes", type = "character", default = NULL,
              help = "optional file of gene ids to restrict to (default: all)")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (r in c("gff_dir", "outdir")) {
  if (is.null(opt[[r]])) stop("--", r, " is required", call. = FALSE)
}
## Checked here so a bad setting fails before reading tens of thousands of
## GFFs; intronic_parts() enforces it again at the point of use.
if (opt$T < 2L * opt$W) {
  stop("T (", opt$T, ") must be >= 2*W (", 2L * opt$W,
       "): otherwise the two end bins overlap and double-count reads",
       call. = FALSE)
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## integers as plain digits, never "2.3e+07"
fi <- function(x) format(as.numeric(x), scientific = FALSE, trim = TRUE)

gff_files <- list.files(opt$gff_dir, pattern = "\\.dexseq\\.gff$", full.names = TRUE)
if (!length(gff_files)) stop("no *.dexseq.gff files in ", opt$gff_dir, call. = FALSE)
if (!is.null(opt$genes)) {
  want <- readLines(opt$genes)
  want <- want[nzchar(want)]
  keep <- sub("\\.dexseq\\.gff$", "", basename(gff_files)) %in% want
  gff_files <- gff_files[keep]
  message("restricted to ", length(gff_files), " of ", length(want), " requested genes")
}
message("genes: ", length(gff_files))

rows <- vector("list", length(gff_files))
n_overlap <- 0L; n_single <- 0L; n_nointron <- 0L; k <- 0L
for (fp in gff_files) {
  gene <- sub("\\.dexseq\\.gff$", "", basename(fp))
  ln <- readLines(fp, warn = FALSE)
  ln <- ln[grepl("\texonic_part\t", ln, fixed = TRUE)]
  if (!length(ln)) { n_nointron <- n_nointron + 1L; next }
  f <- strsplit(ln, "\t", fixed = TRUE)
  chrom  <- vapply(f, `[`, "", 1L)
  start  <- as.numeric(vapply(f, `[`, "", 4L))
  end    <- as.numeric(vapply(f, `[`, "", 5L))
  strand <- vapply(f, `[`, "", 7L)
  attr9  <- vapply(f, `[`, "", 9L)
  part   <- sub('.*exonic_part_number "([^"]+)".*', "\\1", attr9)

  ## The tiling itself is grase::intronic_parts() -- the disjointness
  ## precondition, the no-gap filter and the long-intron end bins all live
  ## there and are unit tested. A gene it cannot tile comes back with zero rows
  ## and a `skip` reason, which is tallied rather than treated as an error:
  ## overlapping exonic parts break the premise of the tiling, so that gene is
  ## excluded instead of silently yielding negative-width bins.
  b <- intronic_parts(start, end, part, max_whole = opt$T, end_width = opt$W)
  if (!nrow(b)) {
    sk <- attr(b, "skip")
    if (identical(sk, "overlap"))     n_overlap  <- n_overlap + 1L
    else if (identical(sk, "single")) n_single   <- n_single + 1L
    else                              n_nointron <- n_nointron + 1L
    next
  }

  k <- k + 1L
  rows[[k]] <- data.frame(
    bin_id = b$bin_id, gene = gene, chrom = chrom[1], strand = strand[1],
    intron_id = b$intron_id, touches = b$touches,
    start = b$start, end = b$end, width = b$width,
    intron_start = b$intron_start, intron_end = b$intron_end,
    intron_len = b$intron_len,
    abut_part_lo = b$abut_part_lo, abut_part_hi = b$abut_part_hi,
    stringsAsFactors = FALSE)
}
rows <- rows[seq_len(k)]
if (!k) stop("no intronic parts emitted -- check --gff_dir", call. = FALSE)
## rbindlist over tens of thousands of frames; do.call(rbind, .) is quadratic
man <- if (requireNamespace("data.table", quietly = TRUE)) {
  as.data.frame(data.table::rbindlist(rows), stringsAsFactors = FALSE)
} else do.call(rbind, rows)

if (anyDuplicated(paste(man$gene, man$bin_id))) {
  stop("duplicate bin_id within a gene -- dexseq_count names features ",
       "gene_id:exonic_part_number, so ids must be unique per gene", call. = FALSE)
}
if (any(man$width < 1)) stop("non-positive bin width", call. = FALSE)

message("intronic parts: ", length(unique(paste(man$gene, man$intron_id))),
        "  bins: ", nrow(man), "  genes: ", length(unique(man$gene)))
message("  whole-part bins: ", sum(man$touches == "both"),
        "   end bins: ", sum(man$touches != "both"))
u <- unique(man[, c("gene", "intron_id", "intron_len")])
message("  sequence counted: ", round(sum(man$width) / 1e6, 1), " Mb of ",
        round(sum(u$intron_len) / 1e6, 1), " Mb intronic")
message("  skipped: ", n_overlap, " genes with overlapping exonic parts, ",
        n_single, " single-part, ", n_nointron, " no intronic part")

## ---- DEXSeq-format GFF ------------------------------------------------------
## Only `exonic_part` lines are read by dexseq_count.py, which names each
## feature gene_id + ":" + exonic_part_number by plain string concatenation, so
## the count rowname is `<gene>:I<lo>_<hi>[.lo|.hi]` and maps back to the
## abutting exonic parts with no lookup table. `aggregate_gene` lines are
## ignored by the counter but written for format completeness.
ord <- order(man$chrom, man$start, man$end)
mo <- man[ord, , drop = FALSE]
src <- "intronparts.R"
## One pass with split(). The per-gene `mo[mo$gene == g, , drop = FALSE]` form
## inherited from flankbins.R is quadratic and does not survive this scale:
## 37,689 genes each rescanning 414,993 rows made no progress in 5 minutes.
idx <- split(seq_len(nrow(mo)), mo$gene)
idx <- idx[order(vapply(idx, `[`, integer(1), 1L))]   # keep chrom/start gene order
exon_lines <- paste(mo$chrom, src, "exonic_part", fi(mo$start), fi(mo$end), ".",
                    mo$strand, ".",
                    sprintf('gene_id "%s"; transcripts "NA"; exonic_part_number "%s"',
                            mo$gene, mo$bin_id), sep = "\t")
agg_lines <- vapply(idx, function(i)
  paste(mo$chrom[i[1L]], src, "aggregate_gene", fi(min(mo$start[i])),
        fi(max(mo$end[i])), ".", mo$strand[i[1L]], ".",
        sprintf('gene_id "%s"', mo$gene[i[1L]]), sep = "\t"), character(1))
out <- vector("list", length(idx))
for (j in seq_along(idx)) out[[j]] <- c(agg_lines[j], exon_lines[idx[[j]]])

gff_path <- file.path(opt$outdir, "intronparts.dexseq.gff")
writeLines(unlist(out, use.names = FALSE), gff_path)

man_path <- file.path(opt$outdir, "intronparts_manifest.tsv")
mo2 <- mo
for (cc in c("start", "end", "width", "intron_start", "intron_end", "intron_len")) {
  mo2[[cc]] <- fi(mo2[[cc]])
}
utils::write.table(mo2, man_path, sep = "\t", quote = FALSE, row.names = FALSE)

message("wrote ", gff_path)
message("wrote ", man_path)
message("next: count it with the SAME settings as the exonic parts")
