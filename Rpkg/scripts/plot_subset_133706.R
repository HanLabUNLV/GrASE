# E026-E038 (the changed transcript region) spans sg_id 51-67.
# Adjust n_min / n_max to change the window.
library(igraph)
library(tidyverse)
library(GenomicFeatures)
library(SplicingGraphs)
library(grase)

gene    <- "ENSG00000133706.17"
indir   <- "/mnt/data1/home/mirahan/GrASE_simulation"
outdir  <- "/mnt/data1/home/mirahan/GrASE"
n_min   <- 35
n_max   <- 80

node_range  <- c(n_min, n_max)
sub_label   <- paste0(gene, ".sub", n_min, "_", n_max)

# --- splicing graph ---------------------------------------------------------
g <- igraph::read_graph(
  file.path(indir, "graphml", paste0(gene, ".graphml")), format = "graphml")
style_and_plot(g, sub_label, outdir, node_range = node_range)
cat("Graph written:", file.path(outdir, paste0(sub_label, ".pdf")), "\n")


# --- tx plot ----------------------------------------------------------------
gr   <- rtracklayer::import(file.path(indir, "gtf", paste0(gene, ".gtf")))
gr   <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
txdb <- txdbmaker::makeTxDbFromGRanges(gr)
sg   <- SplicingGraphs::SplicingGraphs(txdb, min.ntx = 1)

gene_edges <- as.data.frame(SplicingGraphs::sgedges(sg[gene]))
gene_nodes <- SplicingGraphs::sgnodes(sg[gene])

plottx(sub_label, outdir, gene_edges, gene_nodes, node_range = node_range, g = g)
cat("Tx plot written:", file.path(outdir, paste0(sub_label, ".tx.pdf")), "\n")
