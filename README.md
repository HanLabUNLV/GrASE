# GrASE
Graph of Alternative Splice junctions and Exonic parts

## Creating graphMLs
GraphML is an XML-based file format for graphs. This step ensures that the coordinates of each gene are annotated and will be essential when creating the splicing graphs. It uses the [SplicingGraphs](https://bioconductor.org/packages/release/bioc/html/SplicingGraphs.html) R package to obtain the splice junctions and exons.

To generate the GraphML objects for each gene, run `SplicingGraphs.igraph.r`. Please, ensure that the `gencode.TxDb.R` has been sourced. 



* need dexseq results using non-aggregated gff (add -r no to prepare_annotation.py command)
* if you rerun grase.py, empty the files in the results directory `ls grase_results/results/* | while read line; do > $line; done` 
