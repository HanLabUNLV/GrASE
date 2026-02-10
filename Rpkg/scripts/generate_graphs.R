library(dplyr)
library(parallel)
library(grase)
library(SplicingGraphs)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="Input directory path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

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

num_cores <- 20

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

results <- mclapply(genes$V1, process_gene, mc.cores = num_cores)
