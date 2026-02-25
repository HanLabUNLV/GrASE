#' Find the shortest edges in an igraph layout for curvature adjustment
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' coords <- igraph::layout_with_sugiyama(g)$layout
#' pairs <- grase::shortest_edges(g, coords)
#' }
shortest_edges <- function(g, coords) {
  nm <- igraph::V(g)$name
  D <- as.matrix(dist(coords))
  diag(D) <- Inf
  dmin <- min(D, na.rm = TRUE)
  tol  <- max(1e-9, dmin * 1e-6)       # small tolerance for float ties
  pairs <- which(abs(D - dmin) <= tol, arr.ind = TRUE)
  pairs <- pairs[pairs[,1] < pairs[,2], , drop = FALSE]  # drop duplicates (j,i)
  pairs <- cbind.data.frame(nm[pairs[,1]], nm[pairs[,2]])
  pairs
}

# Pull a node toward the average y of its sg_id neighbors (sg-1 and sg+1).
# alpha controls how strongly it moves (0 = no move, 1 = snap to neighbor mean).
pull_vertical_outlier_toward_neighbors <- function(g, coords, target = NULL, alpha = 0.6) {
  sid <- as.integer(igraph::V(g)$sg_id)
  nm  <- igraph::V(g)$name

  k <- if (is.null(target)) {
    # pick the biggest |z| outlier in y
    y <- coords[,2]; z <- (y - stats::median(y)) / stats::mad(y, constant = 1.4826)
    which.max(abs(z))
  } else {
    match(target, nm, nomatch = NA_integer_)
  }
  if (is.na(k)) return(coords)

  y_prev <- coords[which(sid == sid[k] - 1)[1], 2]
  y_next <- coords[which(sid == sid[k] + 1)[1], 2]
  y_ref  <- mean(c(y_prev, y_next), na.rm = TRUE)

  if (is.finite(y_ref)) {
    coords[k,2] <- (1 - alpha) * coords[k,2] + alpha * y_ref
  }
  coords
}

# increase curvature for shortest edges
increase_curve_shortest_edges <- function(g, coords, ec) {
  nm <- igraph::V(g)$name
  el <- igraph::as_edgelist(g)

  pairs = shortest_edges(g, coords)

  # --- select edges that connect ANY closest pair ---
  if (igraph::is_directed(g)) {
    # allow both directions
    edge_keys  <- paste(el[,1], el[,2], sep = "->")
    pair_keys  <- c(paste(pairs[,1], pairs[,2], sep = "->"),
                    paste(pairs[,2], pairs[,1], sep = "->"))
    sel <- edge_keys %in% pair_keys
  } 
  hit <- sum(sel)

  # --- double curvature only for those edges (keep sign) ---
  if (hit > 0) {
    ec[sel] <- pmax(pmin(ec[sel] * 2, 0.95), -0.95)   # clamp to safe range
  }
  ec
}

style_and_plot <- function(g, gene, outdir) {
  color_dict <- c("ex" = "darkorange", "in" = "black", "R" = "black", "L" = "black", "ex_part" = "dark green")
  width_dict <- c("ex" = 1, "in" = 1.2, "R" = 1.2, "L" = 1.2, "ex_part" = 1)
  lty_dict <- c("ex" = "solid", "in" = "dotted", "R" = "dotted", "L" = "dotted", "ex_part" ="solid")

  layout.kamada.kawai.deterministic <- function(...)
  {   
    set.seed(33L)
    igraph::layout.kamada.kawai(...)
  }
   
  grDevices::pdf(file = file.path(outdir, paste0(gene, ".pdf")), width = 32, height = 8) 

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
  vsize <- 2 # tweak if labels feel tight

  # Inspect
  sug_lr <- pull_vertical_outlier_toward_neighbors(g, sug_lr, alpha = 0.7)
  base_curve <- 0.30
  ec <- igraph::curve_multiple(g, base_curve)                    
  ec [ec == 0] = base_curve
  ec <- increase_curve_shortest_edges(g, sug_lr, ec)

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
    vertex.color       = grDevices::adjustcolor("white", alpha.f = 1), # filled
    vertex.label       = lbl,
    vertex.label.cex   = 1,
    vertex.label.color = "black",
    vertex.label.dist  = 0,

    # edges: all curved + green labels
    edge.curved        = ec,              # all edges curved
    edge.label         = igraph::E(g)$dexseq_fragment,
    edge.label.cex     = 1,
    edge.label.color   = "darkgreen",       # green edge labels
    edge.arrow.size    = 0.4,
    edge.color         = sapply(igraph::E(g)$ex_or_in, function(x) color_dict[x]),
    edge.width         = sapply(igraph::E(g)$ex_or_in, function(x) width_dict[x]),
    edge.lty = sapply(igraph::E(g)$ex_or_in, function(x) lty_dict[x])

  )
  grDevices::dev.off()
  return(0)
}

#' Style and tx structure
#' @export
#' @examples
#' \dontrun{
#' library(SplicingGraphs)
#' gr <- rtracklayer::import("ENSG00000000003.gtf")
#' gr <- gr[!(rtracklayer::mcols(gr)$type %in% c("start_codon", "stop_codon"))]
#' sg <- SplicingGraphs::SplicingGraphs(txdbmaker::makeTxDbFromGRanges(gr), min.ntx = 1)
#' gene_edges <- as.data.frame(SplicingGraphs::sgedges(sg["ENSG00000000003.14"]))
#' gene_nodes <- SplicingGraphs::sgnodes(sg["ENSG00000000003.14"])
#' grase::plottx("ENSG00000000003.14", outdir = tempdir(),
#'               gene_edges = gene_edges, gene_nodes = gene_nodes)
#' }
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
  points(x, y+0*2, cex=3, col="black")
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
    segments(x0 = x0, y0 = y, x1 = x1, y1 = y, col="darkorange", lwd=2)
    count = count + 1
  }

  grDevices::dev.off()
  return(0)

}
