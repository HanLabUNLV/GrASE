#devtools::document()
#devtools::install()

#.libPaths('/home/mhan/R/singularity-library/4.4.2/Bioc')
library(dplyr)
#library(grase)
library(SplicingGraphs)
library(txdbmaker)
library(GenomicFeatures)
library(devtools)

DEBUG_MODE = TRUE

grase_package_path <- "/home/mirahan/adjexon_bubbles/GrASE/Rpkg/"
#grase_package_path <- "/Users/mirahan/Work/GrASE/Rpkg/"
#detach("package:grase", unload = TRUE, force = TRUE)
devtools::load_all(grase_package_path)


indir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
outdir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
#indir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'
#outdir = '/data2/han_lab/stevepark/SplicingGraphs/indir/'

genes = c('ENSG00000023318.8', 'ENSG00000023608.5', 'ENSG00000288095.1', 'ENSG00000288600.1')
#gene = 'ENSG00000006744.19'
#gene = 'ENSG00000000003.15'
#gene = 'ENSG00000016490.16'
#gene = 'ENSG00000005302.19'
#gene = 'ENSG00000000460.17'
#gene = 'ENSG00000001084.13'
#gene = 'ENSG00000248933.2'
#gene = 'ENSG00000248672.5'
#gene = 'ENSG00000103254.10'
#gene = 'ENSG00000169488.6'
#gene = 'ENSG00000023191.17'
#gene = 'ENSG00000183878.15'
#gene = 'ENSG00000130830.15'

for (gene in genes) {
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
    next
  }

  gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
  gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
  graph_path <- file.path(indir, "graphml", paste0(gene, ".graphml"))

  gr <- rtracklayer::import(gtf_path)
  if (length(unique(gr$transcript_id[!is.na(gr$transcript_id)])) < 2) {
    message(paste("skipping single isoform ", filename))
    flush.console()
    next
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
  grase::style_and_plot(g, gene, outdir) 

  # call your downstream function, fully namespaced
  cat("  calling grase::bipartition_paths()\n"); flush.console()
   
  bipartitiondir = file.path(outdir, "bipartitions.nocollapse") 
  if(DEBUG_MODE) print("bipartitiondir")
  if(DEBUG_MODE) print(bipartitiondir)
  bipartitions <- grase::bipartition_paths(
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

#  return (gene)

}


