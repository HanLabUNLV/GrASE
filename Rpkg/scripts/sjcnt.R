library(parallel)
library(tidyverse)
library(grase)
library(optparse)

# default values (override with command-line args)
splits_file   <- NULL
sj_dir        <- NULL
cond1         <- "group1"
cond2         <- "group2"
output_prefix <- NULL
analysis_type <- "bipartition"
use_multi     <- FALSE

args <- commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("-s", "--splits"),  type="character", metavar="character",
              help="labeled splits file (output of intron_label.R)"),
  make_option(c("-j", "--sjdir"),   type="character", metavar="character",
              help="directory with condition subdirs containing *.SJ.out.tab files"),
  make_option(c("-1", "--cond1"),   type="character", metavar="character",
              help="condition 1 name (must match subdir name under sjdir)"),
  make_option(c("-2", "--cond2"),   type="character", metavar="character",
              help="condition 2 name (must match subdir name under sjdir)"),
  make_option(c("-o", "--output"),  type="character", metavar="character",
              help="output directory path prefix"),
  make_option(c("-t", "--type"),    type="character", metavar="character",
              default="bipartition",
              help="partition type: bipartition | n_choose_2 | multinomial [default: bipartition]"),
  make_option(c("--multi"),         action="store_true", default=FALSE,
              help="add n_multi to n_uniq counts (default: n_uniq only)")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)
print(opt)

if (!is.null(opt$splits))  splits_file   <- path.expand(opt$splits)
if (!is.null(opt$sjdir))   sj_dir        <- path.expand(opt$sjdir)
if (!is.null(opt$cond1))   cond1         <- opt$cond1
if (!is.null(opt$cond2))   cond2         <- opt$cond2
if (!is.null(opt$output))  output_prefix <- path.expand(opt$output)
if (!is.null(opt$type))    analysis_type <- opt$type
use_multi <- isTRUE(opt$multi)

if (is.null(splits_file) || is.null(sj_dir) || is.null(output_prefix)) {
  stop("--splits, --sjdir, and --output are required")
}

cat("Reading labeled splits file:", splits_file, "\n")
splits_df <- read.table(splits_file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                        colClasses="character", na.strings=c("NA", ""))

sj_dir1 <- file.path(sj_dir, cond1)
sj_dir2 <- file.path(sj_dir, cond2)

cat("Building SJ matrix for", cond1, "from:", sj_dir1, "\n")
sj1 <- build_sj_matrix(sj_dir1, use_multi=use_multi)

cat("Building SJ matrix for", cond2, "from:", sj_dir2, "\n")
sj2 <- build_sj_matrix(sj_dir2, use_multi=use_multi)

# prefix column names with condition to match exoncnt convention
colnames(sj1$mat) <- paste0(cond1, "_", sj1$samples)
colnames(sj2$mat) <- paste0(cond2, "_", sj2$samples)

# combined matrix: union of all junction keys, zeros for junctions absent in one condition
all_keys  <- union(rownames(sj1$mat), rownames(sj2$mat))
sj_mat    <- matrix(0L, nrow=length(all_keys), ncol=ncol(sj1$mat) + ncol(sj2$mat),
                    dimnames=list(all_keys, c(colnames(sj1$mat), colnames(sj2$mat))))
sj_mat[rownames(sj1$mat), colnames(sj1$mat)] <- sj1$mat
sj_mat[rownames(sj2$mat), colnames(sj2$mat)] <- sj2$mat

sample_names <- colnames(sj_mat)
conditions   <- c(rep(cond1, ncol(sj1$mat)), rep(cond2, ncol(sj2$mat)))
sampleinfo   <- setNames(conditions, sample_names)

if (!dir.exists(output_prefix)) dir.create(output_prefix, recursive=TRUE)

splits_list <- split(splits_df, splits_df$gene)

res <- mclapply(names(splits_list), function(gid) {
  tryCatch({
    gene_splits <- splits_list[[gid]]
    out_file <- file.path(output_prefix, paste0(gid, ".", analysis_type, ".sjcnt.txt"))
    if (!file.exists(out_file) || file.size(out_file) == 0L) {
      count_bipartition_sj(gene_splits, sj_mat, sampleinfo, output_prefix, analysis_type)
    }
  }, error = function(e) {
    msg <- sprintf("[%s] ERROR in gene %s (PID %d): %s\n",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   gid, Sys.getpid(), conditionMessage(e))
    cat(msg, file="sjcnt.errors.log", append=TRUE)
    NULL
  })
}, mc.cores=32)


# concatenate per-gene output files
cat("\nCombining sjcnt files...\n")
files <- list.files(output_prefix,
                    pattern=paste0("\\.", analysis_type, "\\.sjcnt\\.txt$"),
                    full.names=TRUE)

if (length(files) > 0L) {
  combined_file <- file.path(output_prefix,
                             paste0(analysis_type, ".sjcnt.combined.txt"))
  file.copy(files[1L], combined_file, overwrite=TRUE)
  if (length(files) > 1L) {
    for (i in 2:length(files)) {
      lines <- readLines(files[i])
      if (length(lines) > 1L) write(lines[-1L], file=combined_file, append=TRUE)
    }
  }
  cat("Combined", length(files), "files into", combined_file, "\n")
} else {
  cat("No sjcnt output files found -- check that intron_distinct1/2 columns are populated.\n")
}

cat("Done.\n")
