library(igraph)
library(dplyr)
library(grase)
library(doParallel)
library(foreach)

#/RAID10/mirahan/graphml.dexseq.v34/
genelist_file = '~/graphml.dexseq.v34/grase_results.bipart/all_genes.txt'
genes <- read.table(genelist_file, header=FALSE, sep="\t")


indir = '~/graphml.dexseq.v34/'
outdir = '~/graphml.dexseq.v34/'
grase_output_dir = file.path(outdir, 'grase_results.bipart/')
eventTypes =  c('A3SS', 'A5SS', 'SE', 'RI')

rmats_df = data.frame(matrix(nrow = 0, ncol = 14)) 
colnames(rmats_df) = c('ID', 'GeneID', 'geneSymbol', 'chr', 'strand', 'longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE',  'DexseqFragment', 'DexseqRefFrag', 'bipartID')
for (eventType in eventTypes) {
  out2 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.fromGTF.", eventType, ".txt"))
  write.table(rmats_df, file = out2, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

mapped_df = data.frame(matrix(nrow = 0, ncol = 4)) 
for (eventType in eventTypes) {
  colnames(mapped_df) = c('GeneID', 'DexseqFragment', paste0('rMATS_ID_',eventType), paste0('rMATS_ID_',eventType,'_ref'))
  out4 <- file.path(grase_output_dir, "results", "tmp", paste0("combined.dexseq.", eventType, ".mapped.txt"))
  write.table(mapped_df, file = out4, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

num_cores <- 20
cl <- makeCluster(num_cores, outfile = "rmats.log")
registerDoParallel(cl)

results <- foreach(
  gene        = genes[,1],                  # iterates over your gene vector
  .packages = c("grase", "dplyr", "igraph"),      # any packages you call inside
  .combine  = 'c',                     # or 'rbind', 'list', etc.
  .errorhandling = "pass"
#  .verbose  = TRUE
) %dopar% {

#for (gene in genes[,1]) {
  tryCatch({

    print(gene)
    gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
    gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
    graph_path <- file.path(indir, "graphml", paste0(gene, ".graphml"))
    bipartitions_path = file.path(outdir, "bipartitions.filtered", paste0(gene, ".bipartitions.internal.txt"))
    if (!file.exists(bipartitions_path)) { return(NULL) }
    if (file.info(bipartitions_path)$size <= 1) { return(NULL) }

    g <- igraph::read_graph(graph_path, format = "graphml")

    rmats_outdir = paste0(indir, '/grase_results.bipart/gene_files/', gene)
    fromGTF_A3SS = paste0(rmats_outdir, '/', 'fromGTF.A3SS.txt')
    fromGTF_A5SS = paste0(rmats_outdir, '/', 'fromGTF.A5SS.txt')
    fromGTF_SE = paste0(rmats_outdir, '/', 'fromGTF.SE.txt')
    fromGTF_RI = paste0(rmats_outdir, '/', 'fromGTF.RI.txt')

    bipartitions <- read.table(bipartitions_path, sep="\t", header=TRUE)

    if (nrow(bipartitions)) {
      map_rMATS(g, gene, gff_path, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, bipartitions, grase_output_dir)
    }
  }, error = function(e) {
    message("ERROR in gene ", gene, ": ", e$message)
    return(NULL)
  })
}

stopCluster(cl)
