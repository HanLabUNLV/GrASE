#devtools::document()
#devtools::install()

#.libPaths('/home/mhan/R/singularity-library/4.4.2/Bioc')
library(dplyr)
library(doParallel)
library(foreach)
library(grase)
library(SplicingGraphs)

DEBUG_MODE = TRUE

indir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
outdir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
#indir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'
#outdir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'


#sim_deds <- read.table("simulation_deds.cell150.txt", header=TRUE, sep="\t")
sim_deds <- read.table(paste0(indir,"genelist"), header=TRUE, sep="\t")
gene_summary = sim_deds %>% group_by(geneID)
genes =  unique(gene_summary$geneID)

num_cores <- 20
cl <- makeCluster(num_cores, outfile = "bipartitions.log")
registerDoParallel(cl)

results <- foreach(
  gene        = genes,                  # iterates over your gene vector
  .packages = c("grase", "rtracklayer", "SplicingGraphs", "txdbmaker", "igraph"),      # any packages you call inside
  .combine  = 'c',                     # or 'rbind', 'list', etc.
  .errorhandling = "pass"
#  .verbose  = TRUE
) %dopar% {

#for (gene in genes) {
  tryCatch({
  if(DEBUG_MODE) print(paste0("START ", gene))
  flush.console()
  indir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
  outdir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
#  indir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'
#  outdir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'


  filename = file.path(outdir, "bipartitions.nocollapse", paste0(gene, ".bipartitions.txt"))
  runninglog = file.path(outdir, "bipartitions.nocollapse", paste0(gene, ".running"))
  if (file.exists(filename) | file.exists(runninglog)) {
    message(paste("skipping existing ", filename))
    flush.console()
    return(gene)
    #next
  }

  gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
  gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
  graph_path <- file.path(indir, "graphml", paste0(gene, ".graphml"))

  gr <- rtracklayer::import(gtf_path)
  if (length(unique(gr$transcript_id[!is.na(gr$transcript_id)])) < 2) {
    message(paste("skipping single isoform ", filename))
    flush.console()
    return(gene)
    #next
  }
  if(DEBUG_MODE) print(paste("running", filename))
  file.create(runninglog)
  gr <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
  txdb <- txdbmaker::makeTxDbFromGRanges(gr)

  sg <- SplicingGraphs::SplicingGraphs(txdb, min.ntx = 1)
  #pdf(file.path(indir, 'sgplot', paste0(gene, ".sg.pdf")))
  #plot(SplicingGraphs::sgraph(sg))
  #dev.off()


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
  cat("  calling grase::bipartition_paths()\n"); flush.console()
   
  bipartitiondir = file.path(outdir, "bipartitions.nocollapse") 
  if(DEBUG_MODE) print("bipartitiondir")
  if(DEBUG_MODE) print(bipartitiondir)
  grase::bipartition_paths(
      gene     = gene,
      g        = g,
      sg       = sg,
      outdir   = bipartitiondir, 
      max_path = 20,
      collapse_bubbles = FALSE 
  )
  if(DEBUG_MODE) print(paste0("FINISH ", gene))
  on.exit(unlink(runninglog))
  flush.console()
  }, error = function(e) {
    message("ERROR in gene ", gene, ": ", e$message)
    return(NULL)
  })

#  return (gene)

}
stopCluster(cl)


