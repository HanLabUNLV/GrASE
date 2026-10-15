#devtools::document()
#devtools::install()

#.libPaths('/home/mhan/R/singularity-library/4.4.2/Bioc')
library(grase)

## Paths come from the caller. These were hardcoded to
## /mnt/data1/home/mirahan/GrASE_simulation, which tied the package to one
## dataset on one machine. Genes stay as the list below: this is a figure
## script, not a pipeline stage.
## usage: Rscript plotgraphml.R <indir> [outdir]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("usage: Rscript plotgraphml.R <indir> [outdir]", call. = FALSE)
indir  <- args[1]
outdir <- if (length(args) >= 2) args[2] else indir

#genes = c('ENSG00000000457.14', 'ENSG00000001631.15')
#genes = c('ENSG00000004809.14')
genes = c('ENSG00000000419.12')

for (gene in genes) {

  graph_path <- file.path(indir, "graphml", paste0(gene, ".graphml"))
  g <- igraph::read_graph(graph_path, format = "graphml")
  style_and_plot(g, gene, paste0(outdir, '/graseplot/')) 

}

