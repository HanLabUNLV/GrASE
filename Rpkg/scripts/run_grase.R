#' Run grase API
#'
#' @param task Task to perform ("plot_graphs" or "find_focal_exons")
#' @param gene Gene ID
#' @param indir Input directory containing GTF and GraphML
#' @param outdir Output directory
#' @param tx_ids Vector of transcript IDs (required for "find_focal_exons")
#' @export
run_grase <- function(task, gene, indir, outdir, tx_ids = NULL) {
  gtf_path <- file.path(indir, "gtf", paste0(gene, ".gtf"))
  graph_path <- file.path(indir, "graphml", paste0(gene, ".dexseq.graphml"))
  print(gtf_path)
  print(graph_path)

  g <- igraph::read_graph(graph_path, format = "graphml")
  igraph::vertex_attr(g)$sg_id[1] <- 0
  igraph::vertex_attr(g)$sg_id[length(igraph::vertex_attr(g)$sg_id)] <- as.numeric(igraph::vertex_attr(g)$sg_id[length(igraph::vertex_attr(g)$sg_id)-1])+1
  igraph::vertex_attr(g)$sg_id <- as.numeric(igraph::vertex_attr(g)$sg_id)

  gr <- rtracklayer::import(gtf_path)
  gr <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
  txdb <- txdbmaker::makeTxDbFromGRanges(gr)
  sg <- SplicingGraphs::SplicingGraphs(txdb, min.ntx = 1)

  if (task == "plot_graphs") {
    style_and_plot(g, gene, file.path(outdir, "graseplot"))
    png(filename = file.path(outdir, "sgplot", paste0(gene, ".sg.png")), width = 5000, height = 5000)
    plot(SplicingGraphs::sgraph(sg))
    dev.off()
    return(invisible(TRUE))
  }

  if (task == "find_focal_exons_between_tx") {
    stopifnot(!is.null(tx_ids) && length(tx_ids) == 2)
    return(focal_exons_between_tx(gene, g, sg, tx_ids, outdir))
  }

  if (task == "find_focal_exons_gene") {
    return(focal_exons_gene(gene, g, sg, outdir))
  }
  stop("Unsupported task: ", task)
}

#' Run sg2igraph
#' @export

run_sg2igraph <- function() {

chrs = paste0("chr", c(1:22, "X", "Y", "M"))
for (chr in chrs) {
  #chr = 'chrY'
  seqlevels(txdb) <- c(chr)
  sg <- SplicingGraphs(txdb,  min.ntx=1)
  edges_by_gene <- sgedgesByGene(sg)

  geneIDs <- names(edges_by_gene)
  for (geneID in geneIDs) {
    #geneID = 'ENSG00000129824.16'
    print (geneID)
    gene_sg = sg[geneID]

    gene_graph = edges_by_gene[[geneID]]
    sgigraph = SG2igraph(geneID, gene_sg, gene_graph)
    write_graph(sgigraph, paste0(geneID, ".graphml"), "graphml")
    
    # now add dexseq edges 
    gff_filename = paste0(geneID,".dexseq.gff") 
    gff <- readLines(gff_filename)
    sgigraph = map_DEXSeq_from_gff(sgigraph, gff)
    write_graph(sgigraph, paste0(geneID, ".dexseq.graphml"), "graphml")

  }
  seqlevels(txdb) <- seqlevels0(txdb)
}
}
