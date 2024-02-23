# GrASE
Graph of Alternative Splice junctions and Exonic parts

## Dependencies  
Tested on Ubuntu (20.04.6 LTS)  

R packages:  
* SplicingGraphs (1.40.0)
* igrpah (1.5.1)
* GenomicFeatures (1.52.0)
* AnnotatoinDbi (1.62.2)
* DEXSeq (1.46.0)

Python packages:  
* rMATS  (4.1.1)
* igraph (0.10.6)
* python (3.11.5)
* STAR   (optional - 2.7.10b)


## DEXSeq and rMATS
DEXSeq 
* need dexseq results using non-aggregated gff (add -r no to prepare_annotation.py command)

## Creating graphMLs
GraphML is an XML-based file format for graphs. This step ensures that the coordinates of each gene are annotated and will be essential when creating the splicing graphs. It uses the [SplicingGraphs](https://bioconductor.org/packages/release/bioc/html/SplicingGraphs.html) R package to obtain the splice junctions and exons. 

To generate the GraphML objects for each gene, run [SplicingGraphs.igraph.R](SplicingGraphs.igraph.R). Some command line arguments will be required when you run the R script, including the path to the gtf and the name of the organism associated with the gtf (include genus and species).
```
Rscript SplicingGraphs.igraph.R /path/to/annotation.gtf Genus species
```
## Preparing to run GrASE
GrASE will process every gene in your dataset that produces results in DEXSeq and rMATS. In order to properly run GrASE, some setup needs to be done. Run [creatingFilesByGene.sh](creatingFilesByGene.sh) to create and set up your `grase_results` directory (created in your current working directory), which will hold everything you need to run GrASE. 
```
bash creatingFilesByGene.sh -r /path/to/rmats/results -d /path/to/dexseq_prepare_annotation.py -a /path/to/annotation/file.gtf -g /path/to/graphml/directory -p number_of_threads
```

## Running GrASE

usage:
```
python grase.py -g <gene_files> --rmats <rmats_results_directory> --dexseq <dexseq_results.txt> --nthread <number_of_threads>
```



## Heatmap ?
