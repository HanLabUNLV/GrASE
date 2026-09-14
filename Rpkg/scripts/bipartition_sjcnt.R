library(parallel)
library(tidyverse)
library(igraph)
library(grase)
library(optparse)

option_list <- list(
  make_option(c("-i", "--inputdir"),  type="character", metavar="character",
              help="directory of per-gene bipartition split files"),
  make_option(c("-g", "--graphmldir"), type="character", metavar="character",
              help="directory containing per-gene .graphml files"),
  make_option(c("--gff"),             type="character", metavar="character",
              help="DEXSeq aggregated GFF file (for gene->chromosome map)"),
  make_option(c("-j", "--sjdir"),     type="character", metavar="character",
              help="directory with condition subdirs containing *.SJ.out.tab files"),
  make_option(c("-1", "--cond1"),     type="character", metavar="character",
              help="condition 1 name (subdir under sjdir)"),
  make_option(c("-2", "--cond2"),     type="character", metavar="character",
              help="condition 2 name (subdir under sjdir)"),
  make_option(c("--conditions"),      type="character", default=NULL,
              metavar="character",
              help=paste("comma-separated list of ALL conditions, e.g.",
                         "B,CD4,CD8,NK,MONO.CLASSIC -- overrides --cond1/--cond2.",
                         "Mirrors the same option in exoncnt.R, so a multi-group",
                         "design (DICE has 13 cell types) can be counted in one pass")),
  make_option(c("-o", "--output"),    type="character", metavar="character",
              help="output directory"),
  make_option(c("-t", "--type"),      type="character", metavar="character",
              default="internal",
              help="bipartition type suffix to match: internal | TSSTTS [default: internal]"),
  make_option(c("-c", "--cores"),     type="integer",   default=32L,
              metavar="integer",
              help="number of parallel cores [default: 32]"),
  make_option(c("--multi"),           action="store_true", default=FALSE,
              help="add n_multi to n_uniq counts (default: n_uniq only)"),
  make_option(c("--combine"),         action="store_true", default=FALSE,
              help=paste("also write bipartition.sjcnt.combined.txt.",
                         "OFF by default: merge_exon_sj_counts.R reads the",
                         "PER-GENE files, and only the abandoned exontest.sj.R",
                         "path ever consumed the combined one. On DICE it was",
                         "74.6 GB of output that nothing downstream read"))
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

need <- list(opt$inputdir, opt$graphmldir, opt$gff, opt$sjdir, opt$output)
if (is.null(opt$conditions)) need <- c(need, list(opt$cond1, opt$cond2))
if (any(sapply(need, is.null))) {
  stop("--inputdir, --graphmldir, --gff, --sjdir, --output are required, plus ",
       "either --conditions or both --cond1 and --cond2")
}

input_dir   <- path.expand(opt$inputdir)
graphml_dir <- path.expand(opt$graphmldir)
gff_path    <- path.expand(opt$gff)
sj_dir      <- path.expand(opt$sjdir)
cond1       <- opt$cond1
cond2       <- opt$cond2
all_conditions <- if (!is.null(opt$conditions))
  trimws(strsplit(opt$conditions, ",")[[1]]) else c(cond1, cond2)
all_conditions <- all_conditions[nzchar(all_conditions)]
output_dir  <- path.expand(opt$output)
bp_type     <- opt$type
n_cores     <- opt$cores
use_multi   <- isTRUE(opt$multi)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# --- shared setup (done once) ---

cat("Parsing chromosome map from GFF:", gff_path, "\n")
chr_map <- parse_gff_chr_map(gff_path)

## One SJ matrix per condition, then a single union-keyed matrix over all of
## them. Sample columns are named "<condition>_<file basename>", which is how
## exoncnt.R names its samples -- the merge joins the two on that name, so the
## SJ files must be laid out as <sjdir>/<condition>/<sample>.SJ.out.tab.
sj_list <- lapply(all_conditions, function(cd) {
  d <- file.path(sj_dir, cd)
  cat("Building SJ matrix for", cd, "from:", d, "\n")
  m <- build_sj_matrix(d, use_multi=use_multi)
  colnames(m$mat) <- paste0(cd, "_", m$samples)
  m$mat
})
names(sj_list) <- all_conditions

all_keys <- Reduce(union, lapply(sj_list, rownames))
all_cols <- unlist(lapply(sj_list, colnames), use.names=FALSE)
sj_mat   <- matrix(0L, nrow=length(all_keys), ncol=length(all_cols),
                   dimnames=list(all_keys, all_cols))
for (m in sj_list) sj_mat[rownames(m), colnames(m)] <- m

sample_names <- colnames(sj_mat)
conditions   <- rep(all_conditions, vapply(sj_list, ncol, integer(1)))
sampleinfo   <- setNames(conditions, sample_names)
cat(sprintf("Conditions: %s\n", paste(all_conditions, collapse=", ")))
cat(sprintf("Samples: %d, junctions: %d\n", ncol(sj_mat), nrow(sj_mat)))

# --- per-gene input files ---

pattern   <- paste0("\\.bipartition\\.", bp_type, "\\.txt$")
in_files  <- list.files(input_dir, pattern=pattern, full.names=TRUE)
cat("Found", length(in_files), "bipartition files for type:", bp_type, "\n")
if (length(in_files) == 0L) stop("no matching input files in: ", input_dir)

error_log <- file.path(output_dir, "bipartition_sjcnt.errors.log")

# --- parallel per-gene processing ---

cat("Processing genes on", n_cores, "cores...\n")

res <- mclapply(in_files, function(fpath) {
  gid <- sub(paste0("\\.bipartition\\.", bp_type, "\\.txt$"), "", basename(fpath))
  out_file <- file.path(output_dir, paste0(gid, ".bipartition.sjcnt.txt"))
  if (file.exists(out_file) && file.size(out_file) > 0L) return(invisible(NULL))

  tryCatch({
    splits_df <- read.table(fpath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                            colClasses="character", na.strings=c("NA", ""))
    if (nrow(splits_df) == 0L) return(invisible(NULL))

    labeled <- label_all_bipartition_introns(splits_df, graphml_dir, chr_map)
    count_bipartition_sj(labeled, sj_mat, sampleinfo, output_dir, "bipartition")
  }, error = function(e) {
    msg <- sprintf("[%s] ERROR %s (PID %d): %s\n",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   gid, Sys.getpid(), conditionMessage(e))
    cat(msg, file=error_log, append=TRUE)
    NULL
  })
}, mc.cores=n_cores)

# --- combine per-gene output files ---

if (!isTRUE(opt$combine)) {
  cat("\nSkipping combine (--combine not set).\n")
  cat("Done.\n")
  quit(save="no", status=0)
}
cat("\nCombining output files...\n")
out_files <- list.files(output_dir, pattern="\\.bipartition\\.sjcnt\\.txt$",
                        full.names=TRUE)
out_files <- out_files[!grepl("combined", out_files)]

if (length(out_files) > 0L) {
  combined_file <- file.path(output_dir, "bipartition.sjcnt.combined.txt")
  file.copy(out_files[1L], combined_file, overwrite=TRUE)
  if (length(out_files) > 1L) {
    for (i in 2:length(out_files)) {
      lines <- readLines(out_files[i])
      if (length(lines) > 1L) write(lines[-1L], file=combined_file, append=TRUE)
    }
  }
  cat("Combined", length(out_files), "files into", combined_file, "\n")
} else {
  cat("No sjcnt output files produced -- check error log:", error_log, "\n")
}

cat("Done.\n")
