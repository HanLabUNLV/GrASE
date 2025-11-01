#' Style and plot splicing graph
#' @export
style_and_plot <- function(g, gene, outdir) {
  color_dict <- c("ex" = "purple", "in" = "grey", "R" = "black", "L" = "black", "ex_part" = "dark green")
  width_dict <- c("ex" = 1, "in" = 0.4, "R" = 0.2, "L" = 0.2, "ex_part" = 1)

  layout.kamada.kawai.deterministic <- function(...)
  {   
    set.seed(33L)
    igraph::layout.kamada.kawai(...)
  }
   
  # place vertices in numeric order on a line: R at far left, L at far right
  layout_linear_by_sgid <- function(g) {
  sid <- suppressWarnings(as.integer(igraph::V(g)$sg_id))
  if (length(sid) != igraph::vcount(g) || any(is.na(sid))) {
    stop("Vertex attribute sg_id must be numeric-like and length vcount(g).")
  }
  ord <- order(sid)
  x <- integer(length(sid)); x[ord] <- seq_along(sid)
  cbind(x, y = 0)
  }

  grDevices::pdf(file = file.path(outdir, paste0(gene, ".pdf")), width = 32, height = 8) 

  lin <- layout_linear_by_sgid(g) 
  plot(
    g,  
    layout = lin,
    main   = paste(gene, "(linear by sg_id)"),
    vertex.shape = "none",
    vertex.label = igraph::V(g)$name,     # shows R..L even though we order by sg_id
    vertex.label.cex = 1.5,
    edge.label  = igraph::E(g)$dexseq_fragment,
    edge.label.cex = 1.5,
    edge.arrow.size = 0.5,
    edge.color = sapply(igraph::E(g)$ex_or_in, function(x) color_dict[x]),
    edge.width = sapply(igraph::E(g)$ex_or_in, function(x) width_dict[x]),
    edge.curved = ifelse(igraph::E(g)$ex_or_in == "in", 0.35, 0)  # curve introns a bit 
  )


  layers_num <- as.integer(igraph::V(g)$sg_id)
  sug <- igraph::layout_with_sugiyama(g, layers = layers_num)$layout

  # rotate 90° so flow is LEFT→RIGHT instead of vertical
  # (equivalent to swapping axes and flipping the old x)
  sug_lr <- cbind(sug[,2], -sug[,1])

  # ensure 'R' is to the LEFT of 'L' (flip horizontally if needed)
  iR <- which(igraph::V(g)$name == "R")[1]
  iL <- which(igraph::V(g)$name == "L")[1]
  if (!is.na(iR) && !is.na(iL) && sug_lr[iR,1] > sug_lr[iL,1]) {
    sug_lr[,1] <- -sug_lr[,1]
  }

  # circles sized to fit labels
  lbl   <- as.character(igraph::V(g)$name)
  vsize <- 3 # tweak if labels feel tight

  plot(
    g,  
    layout = sug_lr,
    main   = paste(gene, "(sugiyama by sg_id)"),
    asp = 0,

    # nodes: circles
    vertex.shape       = "circle",
    vertex.size        = vsize,
    vertex.frame.color = "black",
    vertex.frame.width = 1.4,
    vertex.color       = grDevices::adjustcolor("white", alpha.f = 0.00), # hollow
    vertex.label       = lbl,
    vertex.label.cex   = 1,
    vertex.label.color = "black",
    vertex.label.dist  = 0,

    # edges: all curved + green labels
    edge.curved        = 0.30,              # all edges curved
    edge.label         = igraph::E(g)$dexseq_fragment,
    edge.label.cex     = 1,
    edge.label.color   = "darkgreen",       # green edge labels
    edge.arrow.size    = 0.4,
    edge.color         = sapply(igraph::E(g)$ex_or_in, function(x) color_dict[x]),
    edge.width         = sapply(igraph::E(g)$ex_or_in, function(x) width_dict[x])

  )
  grDevices::dev.off()
  return(0)
}

#' Style and tx structure
#' @export
plottx <- function (gene, outdir, gene_edges, gene_nodes) {

  exon_edges <- gene_edges %>% dplyr::filter(ex_or_in %in% "ex")
  tx_list <- c()
  for(i in 1:length(gene_edges$tx_id)){
    tx_list <- c(tx_list, gene_edges$tx_id[[i]])
  }
  tx_list <- unique(tx_list)

  grDevices::pdf(file = file.path(outdir, paste0(gene, ".tx.pdf")), width = 24, height = 6)
  #initiate splicing graph plot
  par(mar=c(7,10,7,1), bty="n")
  plot(0,0, type="n",xaxt="n", yaxt="n", xlab=gene, ylab="",
      xlim=c(1,length(gene_nodes)), ylim=c(-(nrow(gene_edges)/10),length(tx_list)+2))

  x <- 1:length(gene_nodes)
  y <- rep(1, length(gene_nodes))

  #splice junctions as points
  points(x, y+0*2, cex=3, col="orange")
  text(x,y,labels=c("R",1:(length(gene_nodes)-2),"L"))

  iArrows <- igraph:::igraph.Arrows

  for(i in 1:nrow(gene_edges)) {
    #exons
    if(gene_edges$ex_or_in[i]=="ex"){
      iArrows(as.numeric(gene_edges$from[i])+1, 1, as.numeric(gene_edges$to[i])+1, 1,
              h.lwd=0.25, sh.lwd=1, sh.col="darkgreen",
              width=2, size=0.01)
    }
    #introns
    if(gene_edges$ex_or_in[i]=="in"){
      iArrows(as.numeric(gene_edges$from[i])+1, 1, as.numeric(gene_edges$to[i])+1, 1,
              h.lwd=0.25, sh.lwd=1, sh.col="dimgrey",
              width=2, size=0.01, curve=0.3 - (i %% 2), sh.lty=2)
    }
    #introns connecting artificial root node (R)
    if(gene_edges$ex_or_in[i]=="" & gene_edges$from[i]=="R"){
      iArrows(1, 1, as.numeric(gene_edges$to[i])+1, 1,
              h.lwd=0.25, sh.lwd=1, sh.col="dimgrey",
              width=2, size=0.01, curve=0.3 - (i %% 2), sh.lty=2,)
    }
    #introns connecting artificial leaf node (L)
    if(gene_edges$ex_or_in[i]=="" & gene_edges$to[i]=="L"){
      iArrows(as.numeric(gene_edges$from[i])+1, 1, length(gene_nodes), 1,
              h.lwd=0.25, sh.lwd=1, sh.col="dimgrey",
              width=2, size=0.01, curve=0.3 - (i %% 2), sh.lty=2)
    }
  }

  #transcript segments
  count = 0
  for(i in 1:length(tx_list)){
    x0 <- c()
    x1 <- c()
    y <- rep(2+count, length(grep(tx_list[i], exon_edges$tx_id)))
    for(j in grep(tx_list[i], exon_edges$tx_id)){
      x0 <- c(x0, as.numeric(exon_edges$from[j])+1)
      x1 <- c(x1, as.numeric(exon_edges$to[j])+1)
    }
    text(x=-1, y=y, tx_list[i], cex=0.8, xpd = TRUE)
    segments(x0 = x0, y0 = y, x1 = x1, y1 = y, col="darkorchid", lwd=2)
    count = count + 1
  }

  grDevices::dev.off()
  return(0)

}
