library(igraph)
library(dplyr)
library(grase)

#/RAID10/mirahan/graphml.dexseq.v34/
genelist_file = '/RAID10/mirahan/graphml.dexseq.v34/grase_results/all_genes.txt'
genes <- read.table(genelist_file, header=TRUE, sep="\t")


indir = '/RAID10/mirahan/graphml.dexseq.v34/'
outdir = '/RAID10/mirahan/graphml.dexseq.v34/'
grase_output_dir = file.path(outdir, 'grase_results/')
eventTypes =  c('A3SS', 'A5SS', 'SE', 'RI')

rmats_df = data.frame(matrix(nrow = 0, ncol = 4)) 
colnames(rmats_df) = c('GeneID', 'ID', 'DexseqFragment', 'DexseqRefFrag')
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

for (gene in genes[,1]) {

  print(gene)
  gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
  gff_path <- file.path(indir, "dexseq.gff", paste0(gene, ".dexseq.gff"))
  graph_path <- file.path(indir, "graphml", paste0(gene, ".graphml"))
  focalexons_path = file.path(outdir, "focalexons.collapse", paste0(gene, ".focalexons.txt"))
  if (!file.exists(focalexons_path)) { next }
  if (file.info(focalexons_path)$size <= 1) { next }

  g <- igraph::read_graph(graph_path, format = "graphml")

  rmats_outdir = paste0(indir, '/grase_results/gene_files/', gene)
  fromGTF_A3SS = paste0(rmats_outdir, '/', 'fromGTF.A3SS.txt')
  fromGTF_A5SS = paste0(rmats_outdir, '/', 'fromGTF.A5SS.txt')
  fromGTF_SE = paste0(rmats_outdir, '/', 'fromGTF.SE.txt')
  fromGTF_RI = paste0(rmats_outdir, '/', 'fromGTF.RI.txt')
  map_rMATS(g, gene, gff_path, fromGTF_A3SS, fromGTF_A5SS, fromGTF_SE, fromGTF_RI, focalexons_path, grase_output_dir)

}
