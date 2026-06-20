#' Find the shortest edges in an igraph layout for curvature adjustment
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @param coords A numeric matrix of layout coordinates with one row per vertex (as returned by
#'   igraph layout functions such as \code{igraph::layout_with_sugiyama}).
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

#' Style and plot
#' @param gene Character string. Gene identifier.
#' @param outdir Character string. Path to output directory where the PDF plot
#'   is written.
#' @param node_range Integer vector of length 2 giving the lower and upper sg_id bounds. If provided,
#'   only nodes within that range are shown; R and L are excluded.
#' @export
style_and_plot <- function(g, gene, outdir, node_range = NULL) {
  if (!is.null(node_range)) {
    sg_ids <- as.integer(igraph::V(g)$sg_id)
    keep   <- sg_ids >= node_range[1] & sg_ids <= node_range[2]
    g <- igraph::induced_subgraph(g, igraph::V(g)[keep])
  }

  color_dict <- c("ex" = "darkorange", "in" = "black", "R" = "black", "L" = "black", "ex_part" = "dark green")
  width_dict <- c("ex" = 1, "in" = 1.2, "R" = 1.2, "L" = 1.2, "ex_part" = 1)
  lty_dict <- c("ex" = "solid", "in" = "dotted", "R" = "dotted", "L" = "dotted", "ex_part" ="solid")

  layout.kamada.kawai.deterministic <- function(...)
  {   
    set.seed(33L)
    igraph::layout.kamada.kawai(...)
  }
   
  n_v   <- igraph::vcount(g)
  pdf_w <- max(24, n_v * 0.48)
  pdf_h <- max(8, n_v * 0.11)
  grDevices::pdf(file = file.path(outdir, paste0(gene, ".pdf")), width = pdf_w, height = pdf_h)

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

  # circles sized to fit labels; shrinks for denser graphs, caps at 2
  lbl   <- as.character(igraph::V(g)$name)
  vsize <- min(2, max(1.2, 100 / n_v))

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
#' @param gene Character string. Gene identifier.
#' @param outdir Character string. Path to output directory where the PDF plot
#'   is written.
#' @param gene_edges A data frame of splice graph edges for the gene, as
#'   returned by \code{as.data.frame(SplicingGraphs::sgedges(...))}.
#' @param gene_nodes A character vector of splice graph node names for the
#'   gene, as returned by \code{SplicingGraphs::sgnodes(...)}.
#' @param node_range Integer vector of length 2 giving the lower and upper sg_id bounds. If provided,
#'   only nodes within that range are shown. R and L edges
#'   are dropped and original node numbers are used as axis labels.
#' @param g Optional igraph splicing graph (from \code{igraph::read_graph}).
#'   If provided, exonic part labels (DEXSeq fragment numbers) are drawn below
#'   the node row.
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
plottx <- function (gene, outdir, gene_edges, gene_nodes, node_range = NULL, g = NULL) {

  if (!is.null(node_range)) {
    n_min    <- node_range[1]
    n_max    <- node_range[2]
    from_num <- suppressWarnings(as.integer(gene_edges$from))
    to_num   <- suppressWarnings(as.integer(gene_edges$to))
    # drop R/L edges; keep only edges fully within the numeric range
    keep       <- !is.na(from_num) & !is.na(to_num) &
                  from_num >= n_min & to_num <= n_max
    gene_edges <- gene_edges[keep, ]
    gene_nodes <- as.character(n_min:n_max)
    # x positions use original sg_id: node "k" plots at x = k + 1
    x_node  <- as.numeric(gene_nodes) + 1
    x_lim   <- c(n_min + 1, n_max + 1)
    x_label <- gene_nodes          # show original node numbers on axis
    x_label_offset <- n_min        # for transcript label placement
  } else {
    x_node  <- seq_along(gene_nodes)
    x_lim   <- c(1, length(gene_nodes))
    x_label <- c("R", seq_len(length(gene_nodes) - 2), "L")
    x_label_offset <- 0
  }

  exon_edges <- gene_edges %>% dplyr::filter(ex_or_in %in% "ex")
  tx_list <- c()
  for(i in seq_along(gene_edges$tx_id)){
    tx_list <- c(tx_list, gene_edges$tx_id[[i]])
  }
  tx_list <- unique(tx_list)

  n_nodes  <- length(gene_nodes)
  n_tx     <- length(tx_list)
  tx_step  <- 1.0
  pdf_w    <- max(24, n_nodes * 0.27)
  pdf_h    <- min(10, max(6, n_tx * 0.35 + 3))
  grDevices::pdf(file = file.path(outdir, paste0(gene, ".tx.pdf")), width = pdf_w, height = pdf_h)
  par(mar=c(10, 10, 7, 1), bty="n")
  plot(0, 0, type="n", xaxt="n", yaxt="n", xlab=gene, ylab="",
       xlim = x_lim,
       ylim = c(-(nrow(gene_edges)/10), n_tx * tx_step + 2))

  y_nodes <- rep(1, length(x_node))
  points(x_node, y_nodes, cex=2.5, col="black")
  text(x_node, y_nodes, labels=x_label, cex=0.7)

  # exonic part labels from graphml ex_part edges
  if (!is.null(g)) {
    gel <- igraph::as_data_frame(g, what="edges")
    ep  <- gel[gel$ex_or_in == "ex_part" & gel$dexseq_fragment != "", ]
    ep$from_n <- suppressWarnings(as.integer(ep$from))
    ep$to_n   <- suppressWarnings(as.integer(ep$to))
    if (!is.null(node_range)) {
      ep <- ep[!is.na(ep$from_n) & !is.na(ep$to_n) &
               ep$from_n >= node_range[1] & ep$to_n <= node_range[2], ]
    }
    ep$mid_x <- (ep$from_n + ep$to_n) / 2 + 1
    text(ep$mid_x, y = -1.5, labels = paste0("E", ep$dexseq_fragment),
         srt = 90, cex = 0.6, adj = c(1, 0.5))
  }

  iArrows <- igraph_arrows

  for(i in 1:nrow(gene_edges)) {
    if(gene_edges$ex_or_in[i]=="ex"){
      iArrows(as.numeric(gene_edges$from[i])+1, 1, as.numeric(gene_edges$to[i])+1, 1,
              h.lwd=0.25, sh.lwd=1, sh.col="darkgreen",
              width=2, size=0.01)
    }
    if(gene_edges$ex_or_in[i]=="in"){
      iArrows(as.numeric(gene_edges$from[i])+1, 1, as.numeric(gene_edges$to[i])+1, 1,
              h.lwd=0.25, sh.lwd=1, sh.col="dimgrey",
              width=2, size=0.01, curve=0.3 - (i %% 2), sh.lty=2)
    }
    if(is.null(node_range)){
      if(gene_edges$ex_or_in[i]=="" & gene_edges$from[i]=="R"){
        iArrows(1, 1, as.numeric(gene_edges$to[i])+1, 1,
                h.lwd=0.25, sh.lwd=1, sh.col="dimgrey",
                width=2, size=0.01, curve=0.3 - (i %% 2), sh.lty=2)
      }
      if(gene_edges$ex_or_in[i]=="" & gene_edges$to[i]=="L"){
        iArrows(as.numeric(gene_edges$from[i])+1, 1, length(gene_nodes), 1,
                h.lwd=0.25, sh.lwd=1, sh.col="dimgrey",
                width=2, size=0.01, curve=0.3 - (i %% 2), sh.lty=2)
      }
    }
  }

  count <- 0
  for(i in seq_along(tx_list)){
    x0 <- c()
    x1 <- c()
    y_pos <- 2 + count
    y <- rep(y_pos, length(grep(tx_list[i], exon_edges$tx_id)))
    for(j in grep(tx_list[i], exon_edges$tx_id)){
      x0 <- c(x0, as.numeric(exon_edges$from[j])+1)
      x1 <- c(x1, as.numeric(exon_edges$to[j])+1)
    }
    if (length(x0) > 0) {
      text(x = x_lim[1] - 1, y = y_pos, tx_list[i], cex=0.8, xpd=TRUE)
      segments(x0=x0, y0=y, x1=x1, y1=y, col="darkorange", lwd=2)
    }
    count <- count + 1
  }

  grDevices::dev.off()
  return(0)

}

# Adapted from igraph:::igraph.Arrows (igraph package, GPL >= 2)
# https://github.com/igraph/rigraph
igraph_arrows <- function(x1, y1, x2, y2, code = 2, size = 1,
    width = 1.2/4/par("cin")[2], open = TRUE, sh.adj = 0.1,
    sh.lwd = 1, sh.col = par("fg"), sh.lty = 1, h.col = sh.col,
    h.col.bo = sh.col, h.lwd = sh.lwd, h.lty = sh.lty, curved = FALSE) {
  n <- length(x1)
  recycle <- function(x) rep(x, length.out = n)
  x1 <- recycle(x1); y1 <- recycle(y1)
  x2 <- recycle(x2); y2 <- recycle(y2)
  size <- recycle(size); width <- recycle(width); curved <- recycle(curved)
  sh.lwd <- recycle(sh.lwd); sh.col <- recycle(sh.col); sh.lty <- recycle(sh.lty)
  h.col <- recycle(h.col); h.col.bo <- recycle(h.col.bo)
  h.lwd <- recycle(h.lwd); h.lty <- recycle(h.lty)
  uin <- 1/xyinch()
  label_x <- numeric(n)
  label_y <- numeric(n)
  for (i in seq_len(n)) {
    cin <- size[i] * par("cin")[2]
    w <- width[i] * (1.2/4/cin)
    delta <- sqrt(h.lwd[i]) * par("cin")[2] * 0.005
    x <- sqrt(seq(0, cin^2, length.out = floor(35 * cin) + 2))
    x.arr <- c(-rev(x), -x)
    wx2 <- w * x^2
    y.arr <- c(-rev(wx2 + delta), wx2 + delta)
    deg.arr <- c(atan2(y.arr, x.arr), NA)
    r.arr <- c(sqrt(x.arr^2 + y.arr^2), NA)
    theta1 <- atan2((y1[i] - y2[i]) * uin[2], (x1[i] - x2[i]) * uin[1])
    theta2 <- atan2((y2[i] - y1[i]) * uin[2], (x2[i] - x1[i]) * uin[1])
    r.seg <- cin * sh.adj
    x1d <- y1d <- x2d <- y2d <- 0
    if (code %in% c(1, 3)) { x2d <- r.seg * cos(theta2)/uin[1]; y2d <- r.seg * sin(theta2)/uin[2] }
    if (code %in% c(2, 3)) { x1d <- r.seg * cos(theta1)/uin[1]; y1d <- r.seg * sin(theta1)/uin[2] }
    sx1 <- x1[i] + x1d; sy1 <- y1[i] + y1d
    sx2 <- x2[i] + x2d; sy2 <- y2[i] + y2d
    if (!curved[i]) {
      segments(sx1, sy1, sx2, sy2, lwd = sh.lwd[i], col = sh.col[i], lty = sh.lty[i])
      phi <- atan2(y1[i] - y2[i], x1[i] - x2[i])
      r <- sqrt((x1[i] - x2[i])^2 + (y1[i] - y2[i])^2)
      label_x[i] <- x2[i] + 2/3 * r * cos(phi)
      label_y[i] <- y2[i] + 2/3 * r * sin(phi)
    } else {
      lambda <- if (is.numeric(curved)) curved[i] else 0.5
      midx <- (x1[i] + x2[i])/2; midy <- (y1[i] + y2[i])/2
      spx <- midx - lambda * 1/2 * (sy2 - sy1)
      spy <- midy + lambda * 1/2 * (sx2 - sx1)
      spl <- xspline(x = c(sx1, spx, sx2), y = c(sy1, spy, sy2), shape = 1, draw = FALSE)
      lines(spl, lwd = sh.lwd[i], col = sh.col[i], lty = sh.lty[i])
      label_x[i] <- spl$x[round(2/3 * length(spl$x))]
      label_y[i] <- spl$y[round(2/3 * length(spl$y))]
      if (code %in% c(2, 3)) { x1[i] <- spl$x[round(3/4 * length(spl$x))]; y1[i] <- spl$y[round(3/4 * length(spl$y))] }
      if (code %in% c(1, 3)) { x2[i] <- spl$x[round(1/4 * length(spl$x))]; y2[i] <- spl$y[round(1/4 * length(spl$y))] }
    }
    draw_arrowhead <- function(px, py, theta) {
      px2 <- rep(px, length(deg.arr)); py2 <- rep(py, length(deg.arr))
      ttheta <- rep(theta, length(deg.arr)) + deg.arr
      xhead <- px2 + r.arr * cos(ttheta)/uin[1]
      yhead <- py2 + r.arr * sin(ttheta)/uin[2]
      if (open) { lines(xhead, yhead, lwd = h.lwd[i], col = h.col.bo[i], lty = h.lty[i]) }
      else { polygon(xhead, yhead, col = h.col[i], lwd = h.lwd[i], border = h.col.bo[i], lty = h.lty[i]) }
    }
    if (code %in% c(2, 3)) draw_arrowhead(x2[i], y2[i], atan2((y2[i] - y1[i]) * uin[2], (x2[i] - x1[i]) * uin[1]))
    if (code %in% c(1, 3)) draw_arrowhead(x1[i], y1[i], atan2((y1[i] - y2[i]) * uin[2], (x1[i] - x2[i]) * uin[1]))
  }
  list(lab.x = label_x, lab.y = label_y)
}
