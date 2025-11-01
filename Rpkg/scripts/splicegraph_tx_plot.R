## script to make transcript blocks above splicing graph
library(tidyverse)
library(GenomicFeatures)
library(SplicingGraphs)


indir = '/mnt/storage/jaquino/scRNAseq_sim_pt2/grase/graphml.dexseq.v34/'
gene = 'ENSG00000004809.14'
gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
gr <- rtracklayer::import(gtf_path)
gr <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
txdb <- txdbmaker::makeTxDbFromGRanges(gr)

sg <- SplicingGraphs::SplicingGraphs(txdb, min.ntx = 1)
gene_edges <- as.data.frame(SplicingGraphs::sgedges(sg[gene]))
gene_nodes <- SplicingGraphs::sgnodes(sg[gene])

#pdf(file.path(indir, 'sgplot', paste0(gene, ".tx.pdf")), width=24, height=6)
grase::plottx(gene, paste0(indir,'/sgplot'), gene_edges, gene_nodes)
#dev.off()

