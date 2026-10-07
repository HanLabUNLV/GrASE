library(dplyr)
library(parallel)
library(grase)
library(SplicingGraphs)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="Input directory path", metavar="character"),
  make_option(c("--mc_cores"), type="integer", default=NULL,
              help="parallel workers [default: mc.cores option, else min(detectCores(), 8)]",
              metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## Cores: explicit --mc_cores, else the mc.cores option, else a portable cap.
## Never detectCores() unbounded -- this box reports 188, a laptop reports 4,
## and the old hard-coded value was wrong for both.
mc_cores <- if (!is.null(opt$mc_cores)) as.integer(opt$mc_cores) else
  as.integer(getOption("mc.cores", max(1L, min(parallel::detectCores(), 8L))))
message("Using ", mc_cores, " cores.")

if (is.null(opt$indir)) {
  print_help(opt_parser)
  stop("Input directory (--indir) must be specified.", call.=FALSE)
}

indir <- opt$indir
graphdir = paste0(indir, "/graphml/")
if (!dir.exists(graphdir)) {
  dir.create(graphdir)
}
genes <- read.table(paste0(indir,"/ref/genelist"), header=FALSE)
print(head(genes))


process_gene <- function(gene) {
  tryCatch({
    gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
    gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
    graph_path <- file.path(graphdir, paste0(gene, ".graphml"))

    gr <- rtracklayer::import(gtf_path)
    if (length(unique(gr$transcript_id[!is.na(gr$transcript_id)])) < 2) {
      return(NULL)
    }
    gr <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
    
    # Create TxDb in same process, use immediately, then let it be cleaned up normally
    sg <- SplicingGraphs::SplicingGraphs(txdbmaker::makeTxDbFromGRanges(gr), min.ntx = 1)
    rm(gr)  # Free memory
    gc(verbose = FALSE)  # Force cleanup before finalize issues
    
    edges_by_gene <- SplicingGraphs::sgedgesByGene(sg)
    gene_sg = sg[gene]
    gene_graph = edges_by_gene[[gene]]
    sgigraph = grase::SG2igraph(gene, gene_sg, gene_graph)
    
    # now add dexseq edges 
    gff <- readLines(gff_path)
    sgigraph = grase::map_DEXSeq_from_gff(sgigraph, gff)
    igraph::write_graph(sgigraph, graph_path, "graphml")
    
    return(gene)

  }, error = function(e) {
    message(paste("error in ", gene, ": ", e))
    return(paste("ERROR", gene))
  })
}

results <- mclapply(genes$V1, process_gene, mc.cores = mc_cores)
