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
  make_option(c("-o", "--output"),    type="character", metavar="character",
              help="output directory"),
  make_option(c("-t", "--type"),      type="character", metavar="character",
              default="internal",
              help="bipartition type suffix to match: internal | TSSTTS [default: internal]"),
  make_option(c("-c", "--cores"),     type="integer",   default=32L,
              metavar="integer",
              help="number of parallel cores [default: 32]"),
  make_option(c("--multi"),           action="store_true", default=FALSE,
              help="add n_multi to n_uniq counts (default: n_uniq only)")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (any(sapply(list(opt$inputdir, opt$graphmldir, opt$gff,
                    opt$sjdir, opt$cond1, opt$cond2, opt$output), is.null))) {
  stop("--inputdir, --graphmldir, --gff, --sjdir, --cond1, --cond2, --output are all required")
}

input_dir   <- path.expand(opt$inputdir)
graphml_dir <- path.expand(opt$graphmldir)
gff_path    <- path.expand(opt$gff)
sj_dir      <- path.expand(opt$sjdir)
cond1       <- opt$cond1
cond2       <- opt$cond2
output_dir  <- path.expand(opt$output)
bp_type     <- opt$type
n_cores     <- opt$cores
use_multi   <- isTRUE(opt$multi)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# --- shared setup (done once) ---

cat("Parsing chromosome map from GFF:", gff_path, "\n")
chr_map <- parse_gff_chr_map(gff_path)

sj_dir1 <- file.path(sj_dir, cond1)
sj_dir2 <- file.path(sj_dir, cond2)
cat("Building SJ matrix for", cond1, "from:", sj_dir1, "\n")
sj1 <- build_sj_matrix(sj_dir1, use_multi=use_multi)
cat("Building SJ matrix for", cond2, "from:", sj_dir2, "\n")
sj2 <- build_sj_matrix(sj_dir2, use_multi=use_multi)

colnames(sj1$mat) <- paste0(cond1, "_", sj1$samples)
colnames(sj2$mat) <- paste0(cond2, "_", sj2$samples)

all_keys <- union(rownames(sj1$mat), rownames(sj2$mat))
sj_mat   <- matrix(0L, nrow=length(all_keys),
                   ncol=ncol(sj1$mat) + ncol(sj2$mat),
                   dimnames=list(all_keys,
                                 c(colnames(sj1$mat), colnames(sj2$mat))))
sj_mat[rownames(sj1$mat), colnames(sj1$mat)] <- sj1$mat
sj_mat[rownames(sj2$mat), colnames(sj2$mat)] <- sj2$mat

sample_names <- colnames(sj_mat)
conditions   <- c(rep(cond1, ncol(sj1$mat)), rep(cond2, ncol(sj2$mat)))
sampleinfo   <- setNames(conditions, sample_names)

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
