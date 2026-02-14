library(optparse)

split_dir = "~/GrASE_simulation/bipartition2"
outdir = "~/GrASE_simulation/bipartition.filtered"

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--split_dir"), type="character", default=NULL,
              help="graphml directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory path", metavar="character"),
  make_option(c("-s", "--split_type"), type="character", default=NULL,
              help="output directory path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$split_dir)) {
  print_help(opt_parser)
  stop("Input directory (--split_dir) must be specified.", call.=FALSE)
}
if (is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Output directory (--outdir) must be specified.", call.=FALSE)
}
if (is.null(opt$split_type)) {
  print_help(opt_parser)
  stop("Split method (--split_type) must be specified among 'bipartition', 'multinomial' or 'n_choose_2'.", call.=FALSE)
}

split_dir <- opt$split_dir
outdir <- opt$outdir
split_type <- opt$split_type
if (!dir.exists(outdir)) {
  dir.create(outdir)
}



split_files <- list.files(path = split_dir,
                                 pattern = paste0(split_type,".txt$"), full.names=TRUE)
for (f in split_files) {
   
  print(paste0("infile :", f))
  gene <- sub("\\.all\\.txt$", "", basename(f))
  gene <- sub("\\.txt$", "", basename(f))
  filename_internal = file.path(outdir, paste0(gene, ".internal.txt"))
  filename_TSS = file.path(outdir, paste0(gene, ".TSSTTS.txt"))

  split_df <- tryCatch({
    read.table(f, sep="\t", header=TRUE)
  }, error = function(e) {
    data.frame()
  })

 
  if (nrow(split_df) > 0) { 
  TSSidx = split_df$source == 'R' | split_df$sink == 'L'
    split_internal <- split_df[!TSSidx,]
    split_TSS <- split_df[TSSidx,]
    if (nrow(split_internal) > 0) { 
      write.table(split_internal, file = filename_internal, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    if (nrow(split_TSS) > 0) { 
      write.table(split_TSS, file = filename_TSS, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  }
}
