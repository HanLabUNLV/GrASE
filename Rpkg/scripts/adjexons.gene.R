#devtools::document()
#devtools::install()

#.libPaths('/Users/mhan/R/singularity-library/4.4.2/Bioc')
library(dplyr)
library(SplicingGraphs)
library(txdbmaker)
library(GenomicFeatures)
library(devtools)

#grase_package_path <- "/home/mirahan/adjexon_bubbles/GrASE/Rpkg/"
grase_package_path <- "/Users/mirahan/Work/GrASE/Rpkg/"
#detach("package:grase", unload = TRUE, force = TRUE)
devtools::load_all(grase_package_path)


#indir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
#outdir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
indir = '/Users/mirahan/Work/GrASE/Rpkg/'
outdir = '/Users/mirahan/Work/GrASE/Rpkg/'

genes = c('ENSG00000196628.19', 'ENSG00000197912.16', 'ENSG00000249859.11', 'ENSG00000156113.23', 
'ENSG00000253314.7', 'ENSG00000241469.9', 'ENSG00000242086.8', 
'ENSG00000145362.20', 
'ENSG00000127990.19',
'ENSG00000226674.11',
'ENSG00000224078.15', 
'ENSG00000179818.14', 
'ENSG00000109339.24')
genes = c('ENSG00000188227.13')
#genes = c('ENSG00000007372.23')

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
gene = 'ENSG00000183878.15'


for (gene in genes) {
  print(paste0("START ", gene))
  flush.console()

  filename = file.path(outdir, "focalexons", paste0(gene, ".focalexons.txt"))
  if (file.exists(filename)) {
    print (paste("skipping", filename))
    print(paste0("FINISH ", gene))
    flush.console()
    next
  }


  gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
  gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
  graph_path <- file.path(indir, "graphml", paste0(gene, ".dexseq.graphml"))

  gr <- rtracklayer::import(gtf_path)
  if (length(unique(gr$transcript_id[!is.na(gr$transcript_id)])) < 2) {
    print (paste("skipping", filename))
    print(paste0("FINISH ", gene))
    flush.console()
    next
  }
  print (paste("running", filename))
  gr <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
  rtracklayer::genome(gr) <- 'hg38'
  txdb <- txdbmaker::makeTxDbFromGRanges(gr)

  sg <- SplicingGraphs::SplicingGraphs(txdb, min.ntx = 1)
  pdf(file.path(indir, 'sgplot', paste0(gene, ".sg.pdf")))
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
  grase::style_and_plot(g, gene, outdir) 

  # call your downstream function, fully namespaced
  cat("  calling grase::focal_exons_gene_powerset()\n"); flush.console()
  
  
  #focalexons_read <- read.table(filename, sep="\t", header=TRUE)

  focalexons <- grase::focal_exons_gene_powerset(
      gene     = gene,
      g        = g,
      sg       = sg,
      outdir   = outdir,
      max_path = 30,
      collapse_bubbles = TRUE
  )

  print(paste0("FINISH ", gene))
  flush.console()


}

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


