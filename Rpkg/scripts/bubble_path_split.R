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
              help="output directory path", metavar="character"),
  make_option(c("-s", "--split"), type="character", default=NULL,
              help="output directory path", metavar="character"),
  make_option(c("-n", "--genename"), type="character", default=NULL,
              help="single gene name to process (optional)", metavar="character")
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
if (is.null(opt$split)) {
  print_help(opt_parser)
  stop("Split method (--split) must be specified among 'bipartition', 'multinomial' or 'n_choose_2'.", call.=FALSE)
}

graphdir <- opt$graphdir
outdir <- opt$outdir
split <- opt$split
if (!dir.exists(outdir)) {
  dir.create(outdir)
}
# Get gene list from graphml files or use single gene if specified
if (!is.null(opt$genename)) {
  gene_names <- opt$genename
  message(paste("Processing single gene:", gene_names))
} else {
  graphml_files <- list.files(graphdir, pattern = "\\.graphml$", full.names = FALSE)
  gene_names <- sub("\\.graphml$", "", graphml_files)
  message(paste("Processing", length(gene_names), "genes from directory"))
}

num_cores <- 20

split_bipartition <- function(gene) {
  tryCatch({
    graph_path <- file.path(graphdir, paste0(gene, ".graphml"))

    filename = file.path(outdir, "/", paste0(gene, ".bipartition.txt"))
    runninglog = file.path(outdir, "/", paste0(gene, ".running"))
    if (file.exists(filename) | file.exists(runninglog)) {
      message(paste("skipping existing ", filename))
      flush.console()
      return(gene)
    }

    print(paste("running", filename))
    file.create(runninglog)

    g <- igraph::read_graph(graph_path, format = "graphml")

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
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "split_bipartition.errors.log", append = TRUE)
      NULL
    }
  )
}

split_multinomial <- function(gene) {
  tryCatch({
    graph_path <- file.path(graphdir, paste0(gene, ".graphml"))

    filename = file.path(outdir, "/", paste0(gene, ".multinomial.txt"))
    runninglog = file.path(outdir, "/", paste0(gene, ".running"))
    if (file.exists(filename) | file.exists(runninglog)) {
      message(paste("skipping existing ", filename))
      flush.console()
      return(gene)
    }

    print(paste("running", filename))
    file.create(runninglog)

    g <- igraph::read_graph(graph_path, format = "graphml")

    grase::multinomial_paths(
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
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "split_multinomial.errors.log", append = TRUE)
      NULL
    }
  )
}

split_n_choose_2 <- function(gene) {
  tryCatch({
    graph_path <- file.path(graphdir, paste0(gene, ".graphml"))

    filename = file.path(outdir, "/", paste0(gene, ".n_choose_2.txt"))
    runninglog = file.path(outdir, "/", paste0(gene, ".running"))
    if (file.exists(filename) | file.exists(runninglog)) {
      message(paste("skipping existing ", filename))
      flush.console()
      return(gene)
    }

    print(paste("running", filename))
    file.create(runninglog)

    g <- igraph::read_graph(graph_path, format = "graphml")

    grase::n_choose_2_paths(
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
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "split_n_choose_2.errors.log", append = TRUE)
      NULL
    }
  )
}


if ( split == 'bipartition') {
  if (length(gene_names) == 1) {
    results <- split_bipartition(gene_names)
  } else {
    results <- mclapply(gene_names, split_bipartition, mc.cores = num_cores)
  }
} else if (split == 'multinomial') {
  if (length(gene_names) == 1) {
    results <- split_multinomial(gene_names)
  } else {
    results <- mclapply(gene_names, split_multinomial, mc.cores = num_cores)
  }
} else if (split == 'n_choose_2') {
  if (length(gene_names) == 1) {
    results <- split_n_choose_2(gene_names)
  } else {
    results <- mclapply(gene_names, split_n_choose_2, mc.cores = num_cores)
  }
} else {
  stop("Split method (--split) must be specified among 'bipartition', 'multinomial' or 'n_choose_2'.", call.=FALSE)
}


