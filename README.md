# GrASE
Graph of Alternative Splice junctions and Exonic parts

## Dependencies
* igraph
* python
* dexseq
* rmats
* star
* splicingGraphs


## DEXSeq and rMATS
DEXSeq 
* need dexseq results using non-aggregated gff (add -r no to prepare_annotation.py command)

## Creating graphMLs
GraphML is an XML-based file format for graphs. This step ensures that the coordinates of each gene are annotated and will be essential when creating the splicing graphs. It uses the [SplicingGraphs](https://bioconductor.org/packages/release/bioc/html/SplicingGraphs.html) R package to obtain the splice junctions and exons.

To generate the GraphML objects for each gene, run [SplicingGraphs.igraph.r](SplicingGraphs.igraph.r).

## Preparing to run GrASE
GrASE will process every gene in your dataset that produces results in DEXSeq and rMATS. In order to properly run GrASE, some setup needs to be done. Run [creatingFilesByGene.sh](creatingFilesByGene.sh) to set up your `grase_results` directory, which will hold everything you need to run GrASE. 
```
bash creatingFilesByGene.sh -r /path/to/rmats/results -d /path/to/dexseq_prepare_annotation.py -a /path/to/annotation/file.gtf -g /path/to/graphml/directory -p number_of_threads
```

## Running GrASE

* if you run grase.py again, empty the files in the results directory: ```ls grase_results/results/* | while read line; do > $line; done``` and then rerun

## Heatmap "(?)"
