library(AnnotationDbi)
library(GenomicFeatures)

# create TxDb from gencode
#source("gencode.TxDb.R")
#gencode.v34 = gencode_txdb(
#  version = "34",
#  genome = "hg38",
#  chrs = paste0("chr", c(seq_len(22), "X", "Y", "M"))
#)
#saveDb(gencode.v34, 'txdb.gencode34.sqlite')

# load previously saved TxDb 
#txdb = loadDb(file = '../txdb.gencode34.sqlite')

library(SplicingGraphs)
library(igraph)
source("gencode.TxDb.R")

txdb = gencode_txdb()

#Alternatively, you can create a txdb from any organism as long as you have the gtf file
#Example for Drosophila melanogaster:
#txdb = makeTxDbFromGFF("Drosophila_melanogaster.BDGP6.46.111.gtf", format="gtf", organism="Drosophila melanogaster")

SG2igraph <- function(geneID, sg, edges_by_gene) {

  nodes = sgnodes(sg[geneID])
  g1 = edges_by_gene[[geneID]]

  g1.df <- cbind.data.frame( from=g1$from, to=g1$to, g1[,4:5])
  g1.ncol = ncol(g1.df)
  tx_ids = unique(unlist(g1.df$tx_id))
  for (tx_id in tx_ids) {
    g1.df = cbind.data.frame(g1.df, unlist(lapply(lapply(g1.df$tx_id, '%in%', tx_id) , any)))
  }
  colnames(g1.df)[(g1.ncol+1):(g1.ncol+length(tx_ids))] = tx_ids

  if (g1.df$strand[1] == '+') {
    node1 <- data.frame(oldid = g1.df$from, newcoord = g1.df$start)
    node2 <- data.frame(oldid = g1.df$to, newcoord = g1.df$end+1)
    node_coord = unique(rbind.data.frame(node1, node2))
    node_coord = node_coord[order(node_coord$newcoord),]
  }
  else if (g1.df$strand[1] == '-') {
    start = g1.df$start
    end = g1.df$end
    g1.df$start = end
    g1.df$end = start
    node1 <- data.frame(oldid = g1.df$from, newcoord = g1.df$start+1)
    node2 <- data.frame(oldid = g1.df$to, newcoord = g1.df$end)
    node_coord = unique(rbind.data.frame(node1, node2))
    node_coord = node_coord[order(-node_coord$newcoord),]
  }

  rownames(node_coord) = node_coord$oldid
  g1.df$from = node_coord[g1.df$from,"newcoord"]
  g1.df$to = node_coord[g1.df$to,"newcoord"]
  
  g1.df$tx_id <- vapply(g1.df$tx_id, paste, collapse = ",", character(1L))

  write.table(g1.df, paste0(geneID, ".edgetable.txt"),  row.names=FALSE, sep="\t", quote=FALSE, col.names = TRUE)

  node_coord = node_coord[!duplicated(node_coord$newcoord),] 
  nodes.df = data.frame( ID=c('R', node_coord$newcoord, 'L'))
  write.table(nodes.df, paste0(geneID, ".vertices.txt"),  row.names=FALSE, sep="\t", quote=FALSE, col.names = TRUE)

  drops <- c("seqnames","strand", "tx_id")
  g1.df = g1.df[ , !(names(g1.df) %in% drops)]
  g <- graph.data.frame(g1.df, directed=TRUE, vertices=nodes.df)
  write_graph(g, paste0(geneID, ".graphml"), "graphml")

  roots = nodes.df$ID[!(nodes.df$ID %in% c(g1.df$to, "R", "L"))] 
  leaves = nodes.df$ID[!(nodes.df$ID %in% c(g1.df$from, "R", "L"))]
  for (root in roots) {
    g <- g %>% add_edges(c("R", root))
  }
  for (leaf in leaves) {
    g <- g %>% add_edges(c(leaf, "L"))
  }
  return (g)
}





#specify the chromosomes you are interested in
chrs = paste0("chr", c(1:22, "X", "Y", "M"))
#example chromosomes for Drosophila melanogaster:
#chrs = c("2L", "2R", "3L", "3R", "X", "Y")

#generate graphMLs for each gene
for (chr in chrs) {
  seqlevels(txdb) <- c(chr)
  sg <- SplicingGraphs(txdb,  min.ntx=1)
  edges_by_gene <- sgedgesByGene(sg)

  geneIDs <- names(edges_by_gene)
  for (geneID in geneIDs) {
    sgigraph = SG2igraph(geneID, sg, edges_by_gene)
    write_graph(sgigraph, paste0(geneID, ".graphml"), "graphml")

  }
  seqlevels(txdb) <- seqlevels0(txdb)
}

