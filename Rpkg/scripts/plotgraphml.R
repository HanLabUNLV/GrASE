#devtools::document()
#devtools::install()

#.libPaths('/home/mhan/R/singularity-library/4.4.2/Bioc')
library(grase)

indir = '/mnt/data1/home/mirahan/GrASE_simulation/'
outdir = '/mnt/data1/home/mirahan/GrASE_simulation/'

#genes = c('ENSG00000000457.14', 'ENSG00000001631.15')
#genes = c('ENSG00000004809.14')
genes = c('ENSG00000000419.12')

for (gene in genes) {

  graph_path <- file.path(indir, "graphml", paste0(gene, ".graphml"))
  g <- igraph::read_graph(graph_path, format = "graphml")
  style_and_plot(g, gene, paste0(outdir, '/graseplot/')) 

}

