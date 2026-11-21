#!/usr/bin/env Rscript
## Junction sets of every bipartition, from the splice graph alone: for each row
## of <gene>.bipartition.<type>.txt, the intron edges distinct to each side (the
## introns every route of the side crosses and no route of the other side does)
## and the introns shared by both, as label_all_bipartition_introns() defines
## them. Written for every gene, with or without junction reads, as
## <gene>.bipartition.<type>.junctions.txt: gene, event, intron_distinct1,
## intron_distinct2, intron_shared. `event` is the row number, as exoncnt.R
## numbers events. scripts/sigcnt.py counts a side on these introns when its
## exonic distinct set is empty.
##
## Usage: Rscript bipartition_junction_sets.R --inputdir <filtered bipartitions> \
##          --graphmldir <graphml> --gff <aggregate DEXSeq GFF> --type internal|terminal \
##          --output <dir> [--cores N]
suppressPackageStartupMessages({ library(optparse); library(parallel); library(grase) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--inputdir",   type = "character", help = "filtered bipartition files"),
  make_option("--graphmldir", type = "character", help = "per-gene .graphml files"),
  make_option("--gff",        type = "character", help = "aggregate DEXSeq GFF (chromosome of each gene)"),
  make_option("--type",       type = "character", help = "internal or terminal"),
  make_option("--output",     type = "character", help = "output directory"),
  make_option("--cores",      type = "integer", default = 8L, help = "parallel workers [default %default]")
)))
for (k in c("inputdir", "graphmldir", "gff", "type", "output"))
  if (is.null(opt[[k]])) stop("--", k, " is required")
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

chr_map <- parse_gff_chr_map(opt$gff)
pattern <- paste0("\\.bipartition\\.", opt$type, "\\.txt$")
files <- list.files(opt$inputdir, pattern = pattern, full.names = TRUE)
cat("bipartition files:", length(files), "\n")
err_log <- file.path(opt$output, "bipartition_junction_sets.errors.log")

res <- mclapply(files, function(f) {
  gid <- sub(pattern, "", basename(f))
  out <- file.path(opt$output, paste0(gid, ".bipartition.", opt$type, ".junctions.txt"))
  if (file.exists(out) && file.size(out) > 0L) return("skipped")
  tryCatch({
    s <- utils::read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                           colClasses = "character", na.strings = c("NA", ""))
    if (!nrow(s)) return("empty")
    s$event <- seq_len(nrow(s))
    l <- suppressWarnings(label_all_bipartition_introns(s, opt$graphmldir, chr_map))
    utils::write.table(l[, c("gene", "event", "intron_distinct1", "intron_distinct2", "intron_shared")],
                       out, sep = "\t", quote = FALSE, row.names = FALSE)
    "ok"
  }, error = function(e) {
    cat(sprintf("[%s] %s: %s\n", format(Sys.time(), "%F %T"), gid, conditionMessage(e)),
        file = err_log, append = TRUE)
    "error"
  })
}, mc.cores = opt$cores, mc.preschedule = FALSE)
print(table(unlist(res)))
