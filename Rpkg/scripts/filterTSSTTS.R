library(optparse)

bipartitions_dir = "~/GrASE_simulation/bipartition"
outdir = "~/GrASE_simulation/bipartition.filtered"

# Parse command line arguments
option_list <- list(
  make_option(c("-g", "--bipartitions_dir"), type="character", default=NULL,
              help="graphml directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory path", metavar="character"),
  make_option(c("-s", "--split"), type="character", default=NULL,
              help="output directory path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$bipartitions_dir)) {
  print_help(opt_parser)
  stop("Input directory (--bipartitions_dir) must be specified.", call.=FALSE)
}
if (is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Output directory (--outdir) must be specified.", call.=FALSE)
}
if (is.null(opt$split)) {
  print_help(opt_parser)
  stop("Split method (--split) must be specified among 'bipartition', 'multinomial' or 'n_choose_2'.", call.=FALSE)
}

bipartitions_dir <- opt$bipartitions_dir
outdir <- opt$outdir
split <- opt$split
if (!dir.exists(outdir)) {
  dir.create(outdir)
}



bipartitions_files <- list.files(path = bipartitions_dir,
                                 pattern = "bipartitions.txt$", full.names=TRUE)
for (f in bipartitions_files) {
   
  print(paste0("infile :", f))
  gene <- sub("\\.all\\.txt$", "", basename(f))
  gene <- sub("\\.txt$", "", basename(f))
  filename_internal = file.path(outdir, paste0(gene, ".internal.txt"))
  filename_TSS = file.path(outdir, paste0(gene, ".TSSTTS.txt"))

  bipartitions <- tryCatch({
    read.table(f, sep="\t", header=TRUE)
  }, error = function(e) {
    data.frame()
  })

 
  if (nrow(bipartitions) > 0) { 
  TSSidx = bipartitions$source == 'R' | bipartitions$sink == 'L'
    bipartitions_internal <- bipartitions[!TSSidx,]
    bipartitions_TSS <- bipartitions[TSSidx,]
    if (nrow(bipartitions_internal) > 0) { 
      write.table(bipartitions_internal, file = filename_internal, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    if (nrow(bipartitions_TSS) > 0) { 
      write.table(bipartitions_TSS, file = filename_TSS, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  }
}
