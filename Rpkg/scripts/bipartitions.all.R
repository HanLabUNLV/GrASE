library(dplyr)
library(parallel)
library(grase)
library(SplicingGraphs)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("-g", "--graphdir"), type="character", default=NULL,
              help="graphml directory path", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$graphdir)) {
  print_help(opt_parser)
  stop("Input directory (--graphdir) must be specified.", call.=FALSE)
}

if (is.null(opt$outdir)) {
  print_help(opt_parser)
  stop("Output directory (--outdir) must be specified.", call.=FALSE)
}

graphdir <- opt$graphdir
outdir <- opt$outdir
if (!dir.exists(outdir)) {
  dir.create(outdir)
}
# Get gene list from graphml files
graphml_files <- list.files(graphdir, pattern = "\\.graphml$", full.names = FALSE)
gene_names <- sub("\\.graphml$", "", graphml_files)

num_cores <- 20

process_gene <- function(gene) {
  tryCatch({
    graph_path <- file.path(graphdir, paste0(gene, ".graphml"))

    filename = file.path(outdir, "/", paste0(gene, ".bipartitions.txt"))
    runninglog = file.path(outdir, "/", paste0(gene, ".running"))
    if (file.exists(filename) | file.exists(runninglog)) {
      message(paste("skipping existing ", filename))
      flush.console()
      return(gene)
    }

    print(paste("running", filename))
    file.create(runninglog)

    g <- igraph::read_graph(graph_path, format = "graphml")

    cat("  calling grase::bipartition_paths()\n"); flush.console()
    grase::bipartition_paths(
        gene     = gene,
        g        = g,
        outdir   = outdir, 
        max_path = 20,
        collapse_bubbles = FALSE 
    )
    print(paste0("FINISH ", gene))
    on.exit(unlink(runninglog))
    return(gene)

    }, error = function(e) {
      message(paste("error in ", gene, ": ", e))
      return(paste("ERROR", gene))
    })
  }

results <- mclapply(gene_names, process_gene, mc.cores = num_cores)
