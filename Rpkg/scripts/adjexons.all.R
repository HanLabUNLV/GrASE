#devtools::document()
#devtools::install()

.libPaths('/home/mhan/R/singularity-library/4.4.2/Bioc')
library(dplyr)
library(doParallel)
library(foreach)
library(grase)
library(SplicingGraphs)

indir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'
outdir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'

#sim_deds <- read.table("simulation_deds.cell150.txt", header=TRUE, sep="\t")
sim_deds <- read.table(paste0(indir,"genelist"), header=TRUE, sep="\t")
gene_summary = sim_deds %>% group_by(geneID)
genes =  unique(gene_summary$geneID)
#gene = 'ENSG00000006744.19'
#gene = 'ENSG00000000003.15'
#gene = 'ENSG00000016490.16'
#gene = 'ENSG00000005302.19'

num_cores <- 20
cl <- makeCluster(num_cores, outfile = "worker5.log")
registerDoParallel(cl)

results <- foreach(
  gene        = genes,                  # iterates over your gene vector
  .packages = c("grase", "rtracklayer", "SplicingGraphs", "txdbmaker", "igraph"),      # any packages you call inside
  .combine  = 'c'                     # or 'rbind', 'list', etc.
#  .verbose  = TRUE
) %dopar% {

  indir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'
  outdir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'

  filename = file.path(outdir, "focalexons", paste0(gene, ".focalexons.txt"))
  if (file.exists(filename)) {
    print (paste("skipping", filename))
    return(gene)
  }

  gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
  gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
  graph_path <- file.path(indir, "graphml", paste0(gene, ".dexseq.graphml"))

  gr <- rtracklayer::import(gtf_path)
  if (length(unique(gr$transcript_id[!is.na(gr$transcript_id)])) < 2) {
    print (paste("skipping", filename))
    return(gene)
  }
  print (paste("running", filename))
  gr <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
  txdb <- txdbmaker::makeTxDbFromGRanges(gr)
  sg <- SplicingGraphs::SplicingGraphs(txdb, min.ntx = 1)
  pdf(paste0(gene, ".sg.pdf"))
  plot(SplicingGraphs::sgraph(sg))
  dev.off()


  edges_by_gene <- SplicingGraphs::sgedgesByGene(sg)
  gene_sg = sg[gene]
  gene_graph = edges_by_gene[[gene]]
  sgigraph = grase::SG2igraph(gene, gene_sg, gene_graph)
  
  # now add dexseq edges 
  gff <- readLines(gff_path)
  sgigraph = grase::map_DEXSeq_from_gff(sgigraph, gff)
  igraph::write_graph(sgigraph, graph_path, "graphml")


  g <- igraph::read_graph(graph_path, format = "graphml")

  # call your downstream function, fully namespaced
  cat("  calling grase::focal_exons_gene_powerset()\n"); flush.console()
  res <- tryCatch({
    grase::focal_exons_gene_powerset(
      gene     = gene,
      g        = g,
      sg       = sg,
      outdir   = outdir
    )
  }, error = function(e) {
    cat("  !! focal_exons failed: ", e$message, "\n"); flush.console()
    stop(e)
  })

  rm(gr)
  rm(txdb)
  rm(sg)
  rm(g)
  gc()
  return (gene)

}
for (r in results) {
  print (r)
}
stopCluster(cl)

#run_grase("find_focal_exons_between_tx", gene, indir = indir, outdir= outdir, tx_ids )
#run_grase("find_focal_exons_gene_ovr", gene, indir = indir, outdir= outdir)
#run_grase("plot_graphs", gene, indir = indir, outdir= outdir, tx_ids )

# create TxDb from gencode
#source("gencode.TxDb.R")
## Locate file
#gtf_file <- gencode_source_url(
#  version = "34",
#  genome = "hg38"
#}
#gencode.v34 = gencode_txdb(
#  gtf_file,
#  chrs = paste0("chr", c(seq_len(22), "X", "Y", "M"))
#)
#saveDb(gencode.v34, 'txdb.gencode34.sqlite')

# load previously saved TxDb 
#txdb = loadDb(file = '../txdb.gencode34.sqlite')


