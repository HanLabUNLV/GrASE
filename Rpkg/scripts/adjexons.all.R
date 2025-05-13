#devtools::document()
#devtools::install()

library(dplyr)
library(doParallel)
library(foreach)

indir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
outdir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'

#sim_deds <- read.table("simulation_deds.cell150.txt", header=TRUE, sep="\t")
sim_deds <- read.table(paste0(indir,"genelist"), header=TRUE, sep="\t")
gene_summary = sim_deds %>% group_by(geneID)
genes =  unique(gene_summary$geneID)
gene = 'ENSG00000006744.19'
gene = 'ENSG00000000003.15'
gene = 'ENSG00000016490.16'

num_cores <- 10
cl <- makeCluster(num_cores, outfile = "worker4.log")
registerDoParallel(cl)

results <- foreach(
  gene        = genes,                  # iterates over your gene vector
  .packages = c("grase", "rtracklayer", "SplicingGraphs", "txdbmaker", "igraph"),      # any packages you call inside
  .combine  = 'c'                     # or 'rbind', 'list', etc.
#  .verbose  = TRUE
) %dopar% {

  indir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
  outdir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'

  filename = file.path(outdir, "focalexons", paste0(gene, ".focalexons.txt"))
  if (file.exists(filename)) {
    print (paste("skipping", filename))
    return(gene)
  }

  gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
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


  edges_by_gene <- SplicingGraphs::sgedgesByGene(sg)
  gene_sg = sg[gene]
  gene_graph = edges_by_gene[[gene]]
  sgigraph = grase::SG2igraph(gene, gene_sg, gene_graph)
  
  # now add dexseq edges 
  gff_filename = paste0(gene,".dexseq.gff") 
  gff <- readLines(gff_filename)
  sgigraph = grase::map_DEXSeq_from_gff(sgigraph, gff)
  igraph::write_graph(sgigraph, graph_path, "graphml")


  g <- igraph::read_graph(graph_path, format = "graphml")

  # call your downstream function, fully namespaced
  cat("  calling grase::focal_exons_gene_nested()\n"); flush.console()
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


