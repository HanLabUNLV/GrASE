#' Style and plot splicing graph
#' @export
style_and_plot <- function(g, gene, outdir) {
  color_dict <- c("ex" = "purple", "in" = "grey", "R" = "black", "L" = "black", "ex_part" = "dark green")
  width_dict <- c("ex" = 5, "in" = 4, "R" = 2, "L" = 2, "ex_part" = 5)
  #layout <- igraph::layout_with_sugiyama(g)$layout
  layout <- igraph::layout_as_tree(g)$layout

  grDevices::png(filename = file.path(outdir, paste0(gene, ".png")), width = 3500, height = 1000)
  plot(
    g, layout = layout, main = gene,
    edge.color = sapply(igraph::E(g)$ex_or_in, function(x) color_dict[x]),
    edge.width = sapply(igraph::E(g)$ex_or_in, function(x) width_dict[x]),
    vertex.label = igraph::V(g)$sg_id,
    vertex.label.cex = 1.5,
    vertex.shape = "none",
    edge.label = igraph::E(g)$dexseq_fragment,
    edge.label.cex = 1.5,
    edge.arrow.size = 0.001,
    margin = 0.1
  )
  grDevices::dev.off()
  return(0)
}


