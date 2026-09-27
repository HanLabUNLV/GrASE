#!/usr/bin/env Rscript
## Emit flanking bins for the TSS/TTS edge-coverage diagnostic.
##
## Runs AFTER filterTSSTTS.R (bins are defined per bipartition side) and
## BEFORE counting. Emits a SAF for featureCounts plus a manifest tying each
## bin back to its side.
##
## Also emits CALIBRATION bins, without which the step statistic has no scale:
##   positive  the shared TSS of a single-promoter gene -- a real boundary
##   negative  the midpoint of a long exonic part -- pure continuation
## Thresholds are fitted per dataset from these, never hardcoded.
##
## Count the SAF with the SAME semantics as the DEXSeq counting that produced
## the exonic part counts -- read counts, matching strandedness:
##
##   featureCounts -F SAF -a flankbins.saf -s 2 -o flankcounts.txt <bams>
##
## Do NOT fold these bins into the DEXSeq flattening: adding features changes
## read assignment for existing exonic parts, so every current count shifts.
##
## Usage:
##   Rscript scripts/flankbins.R --split_dir=<filtered splits> \
##     --gff_dir=<per-gene dexseq gff> --gencode=<annotation.gff3> \
##     --outdir=<dir> [--width=100] [--min_width=50] [--n_calib=1500]

suppressMessages({library(optparse)})

option_list <- list(
  make_option(c("-i", "--split_dir"), type = "character", default = NULL,
              help = "directory of filtered bipartition split tables"),
  make_option(c("-g", "--gff_dir"), type = "character", default = NULL,
              help = "directory of per-gene <gene>.dexseq.gff files"),
  make_option(c("-a", "--gencode"), type = "character", default = NULL,
              help = "gencode annotation gff3, for single-promoter controls"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "output directory"),
  make_option(c("-w", "--width"), type = "integer", default = 100L,
              help = "nominal flank bin width [default %default]"),
  make_option(c("-m", "--min_width"), type = "integer", default = 50L,
              help = "below this a side is unscorable [default %default]"),
  make_option(c("-n", "--n_calib"), type = "integer", default = 1500L,
              help = "calibration bins per class [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (r in c("split_dir", "gff_dir", "outdir"))
  if (is.null(opt[[r]])) stop("--", r, " must be specified", call. = FALSE)

rdir <- file.path(dirname(dirname(normalizePath(
          sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "R")
if (!dir.exists(rdir)) rdir <- "R"
for (f in list.files(rdir, pattern = "\\.R$", full.names = TRUE))
  try(source(f), silent = TRUE)
stopifnot(exists("extend_contiguous_run"), exists("flank_bin"))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

read_gene_parts <- function(gene) {
  fn <- file.path(opt$gff_dir, paste0(gene, ".dexseq.gff"))
  if (!file.exists(fn)) return(NULL)
  ## quote = "" is REQUIRED: GFF attributes are double-quoted and read.delim's
  ## default quoting would mangle the attribute field.
  x <- utils::read.delim(fn, header = FALSE, comment.char = "#",
                         quote = "", stringsAsFactors = FALSE)
  x <- x[x$V3 == "exonic_part", , drop = FALSE]
  if (!nrow(x)) return(NULL)
  n <- as.integer(sub('.*exonic_part_number "([^"]+)".*', "\\1", x$V9))
  tx <- strsplit(sub('.*transcripts "([^"]+)".*', "\\1", x$V9), "\\+")
  list(parts  = stats::setNames(Map(c, x$V4, x$V5), as.character(n)),
       tx     = stats::setNames(tx, as.character(n)),
       chrom  = x$V1[1], strand = x$V7[1])
}

split_files <- list.files(opt$split_dir, pattern = "\\.txt$", full.names = TRUE)
message("split tables: ", length(split_files))

rows <- list(); saf <- list(); k <- 0L
n_uns <- 0L
for (sf in split_files) {
  d <- try(utils::read.delim(sf, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) next
  if (!all(c("gene", "event", "source", "sink") %in% names(d))) next
  for (i in seq_len(nrow(d))) {
    kind <- if (identical(as.character(d$source[i]), "R")) "TSS" else
            if (identical(as.character(d$sink[i]), "L"))   "TTS" else NA
    if (is.na(kind)) next
    g <- read_gene_parts(d$gene[i]); if (is.null(g)) next
    dir <- free_end_direction(kind, g$strand)
    for (side in c(1L, 2L)) {
      sd <- if (side == 1L) d$setdiff1[i] else d$setdiff2[i]
      if (is.na(sd) || sd %in% c("", "NA")) next
      pn <- suppressWarnings(as.integer(sub("^E", "",
              trimws(strsplit(sd, ",")[[1]]))))
      pn <- pn[!is.na(pn) & as.character(pn) %in% names(g$parts)]
      if (!length(pn)) next
      b <- extend_contiguous_run(g$parts, g$tx, pn, dir)
      if (is.na(b)) next
      fw <- flank_width(g$parts, b, dir)
      fb <- flank_bin(b, dir, opt$width, fw, opt$min_width)
      if (fb$tier == "unscorable") { n_uns <- n_uns + 1L; next }
      k <- k + 1L
      id <- sprintf("F%07d", k)
      len_d <- sum(vapply(g$parts[as.character(pn)],
                          function(p) p[2] - p[1] + 1, numeric(1)))
      rows[[k]] <- data.frame(bin_id = id, row_type = "side", role = "ADJ",
        pair_id = id, gene = d$gene[i],
        event = d$event[i], side = side, kind = kind, chrom = g$chrom,
        strand = g$strand, boundary = b, flank_width = fw,
        len_ADJ = fb$width, len_D = len_d,
        distinct_parts = sd, tier = fb$tier, stringsAsFactors = FALSE)
      saf[[k]] <- data.frame(GeneID = id, Chr = g$chrom, Start = fb$start,
        End = fb$end, Strand = g$strand, stringsAsFactors = FALSE)
    }
  }
}
message("side bins: ", k, "   unscorable (flank < ", opt$min_width, "): ", n_uns)

## ---- calibration bins -------------------------------------------------
## A control needs BOTH windows: the one outside the boundary and the matched
## interior one it is compared against. Without the paired interior window the
## control step has no valid denominator and the fitted threshold is wrong,
## which would mis-set the classification for every real side.
add_control <- function(rt, gene, chrom, strand, adj_s, adj_e, d_s, d_e) {
  pid <- sprintf("P%07d", k + 1L)
  for (role in c("ADJ", "D")) {
    k <<- k + 1L
    id <- sprintf("F%07d", k)
    ss <- if (role == "ADJ") adj_s else d_s
    ee <- if (role == "ADJ") adj_e else d_e
    rows[[k]] <<- data.frame(bin_id = id, row_type = rt, role = role,
      pair_id = pid, gene = gene,
      event = NA_character_, side = NA_integer_, kind = rt, chrom = chrom,
      strand = strand, boundary = NA_integer_, flank_width = opt$width,
      len_ADJ = adj_e - adj_s + 1L, len_D = d_e - d_s + 1L,
      distinct_parts = NA_character_, tier = "full", stringsAsFactors = FALSE)
    saf[[k]] <<- data.frame(GeneID = id, Chr = chrom, Start = ss, End = ee,
      Strand = strand, stringsAsFactors = FALSE)
  }
  pid
}

if (!is.null(opt$gencode) && file.exists(opt$gencode)) {
  message("reading gencode for single-promoter controls ...")
  con <- file(opt$gencode, "r"); tss <- new.env(parent = emptyenv())
  repeat {
    ln <- readLines(con, n = 50000L); if (!length(ln)) break
    ln <- ln[substr(ln, 1, 1) != "#"]
    f <- strsplit(ln, "\t", fixed = TRUE)
    for (x in f) {
      if (length(x) < 9 || x[3] != "transcript") next
      gid <- sub('.*gene_id=([^;]+).*', "\\1", x[9])
      p <- as.integer(if (x[7] == "-") x[5] else x[4])
      tss[[gid]] <- c(tss[[gid]], p)
    }
  }
  close(con)
  sp <- Filter(function(v) (max(v) - min(v)) <= 50L, as.list(tss))
  message("single-promoter genes: ", length(sp))
  npos <- 0L
  for (gene in names(sp)) {
    if (npos >= opt$n_calib) break
    g <- read_gene_parts(gene); if (is.null(g)) next
    if (length(g$parts) < 4L) next
    outer <- sp[[gene]][1]
    dir <- if (g$strand == "+") -1L else 1L
    fb <- flank_bin(outer, dir, opt$width, opt$width, opt$min_width)
    if (fb$tier == "unscorable") next
    ## D is the matched window INSIDE the terminal exon
    if (dir < 0) { ds <- outer; de <- outer + opt$width - 1L }
    else         { ds <- outer - opt$width + 1L; de <- outer }
    add_control("calib_pos", gene, g$chrom, g$strand, fb$start, fb$end, ds, de)
    npos <- npos + 1L
  }
  message("calibration positive bins: ", npos)
}

nneg <- 0L
for (sf in split_files) {
  if (nneg >= opt$n_calib) break
  d <- try(utils::read.delim(sf, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) next
  for (gene in unique(d$gene)) {
    if (nneg >= opt$n_calib) break
    g <- read_gene_parts(gene); if (is.null(g)) next
    lens <- vapply(g$parts, function(p) p[2] - p[1] + 1, numeric(1))
    big <- names(lens)[lens >= 2 * opt$width + 20]
    if (!length(big)) next
    p <- g$parts[[big[ceiling(length(big) / 2)]]]
    mid <- floor((p[1] + p[2]) / 2)
    ## both windows lie inside one continuous exon, so this is pure
    ## continuation and the step must come out near 0.5
    add_control("calib_neg", gene, g$chrom, g$strand,
                mid - opt$width + 1L, mid, mid + 1L, mid + opt$width)
    nneg <- nneg + 1L
  }
}
message("calibration negative bins: ", nneg)

man <- do.call(rbind, rows)
sf_ <- do.call(rbind, saf)
utils::write.table(man, file.path(opt$outdir, "flankbins_manifest.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(sf_, file.path(opt$outdir, "flankbins.saf"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", nrow(man), " bins to ", opt$outdir)
message("next: featureCounts -F SAF -a ", file.path(opt$outdir, "flankbins.saf"),
        " -s 2 -o flankcounts.txt <bams>")
