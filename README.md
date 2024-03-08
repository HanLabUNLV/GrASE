# GrASE
Graph of Alternative Splice junctions and Exonic parts

## About
GrASE is a tool that bridges the gap between the splice-junction approach of software like [rMATS](https://rnaseq-mats.sourceforge.io/index.html) and the exon-fragment approach of software like [DEXSeq](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html). Specifically, this tool uses both rMATS and DEXSeq in order to map splicing events to the exons they encompass and vice versa. Before using GrASE, you must run rMATS on your data in order to get an rmats_output directory. In addition, you must also run DEXSeq on your data, and then output the results into a txt file. GrASE will use the output of both software in order to do its job. 

* Container used for running GrASE can be found here: [GrASE Container](https://drive.google.com/drive/folders/10H6NxN0T1cP0O68VwhCh55KVVb08Iqzb)

## Dependencies  
Tested on Ubuntu (20.04.6 LTS)  

R packages:  
* SplicingGraphs (1.40.0)
* igraph (1.5.1)
* GenomicFeatures (1.52.0)
* AnnotationDbi (1.62.2)
* DEXSeq (1.46.0)

Python packages:  
* python3 (3.11.5)
* rMATS   (4.1.1)
* htseq   (0.13.5)
* igraph  (0.10.6)
* pycairo (1.23.0)
* pandas  (2.1.4)

Other packages:
* STAR   (optional - 2.7.10b)

## DEXSeq and rMATS
### DEXSeq 
DEXSeq is a software that is used for finding differential exon usage using RNA-seq exon counts between samples with different experimental designs. It conducts statistical tests based on a model that uses the negative binomial distribution to estimate the variance between biological replicates and generalized linear models for testing. The output of DEXSeq is a matrix that has many features including coordinates, pvalue, padj, read counts, and more for every exon described in the GFF file. The GFF is prepared by annotating the GTF of a given dataset. GrASE will use the output of DEXSeq in order to map each exon to a graphML object. This will be compared to rMATS results in order to map DEXSeq exons to rMATS events. 

When running DEXSeq, preprocessing should be done with python to ensure full compatibility. When running `dexseq_prepare_annotation.py`, make sure to include `-r no` as an argument. This will ensure that the dexseq.gff output will be in a non-aggregated form (also to ensure compatibility with GrASE). The rest of the DEXSeq pipeline can be completed as normal. Make sure to create a dexseq_output.txt file when the pipeline is finished. GrASE will use this output file in its analysis. 

### rMATS
rMATS is a computational tool to detect differential alternative splicing events from RNA-Seq data. The statistical model of rMATS calculates the P-value and false discovery rate between two conditions/samples exceeds a given threshold. From the RNA-Seq data, rMATS can automatically detect and analyze Alternative Splicing (AS) Events corresponding to all major types of alternative splicing patterns. For the purposes of GrASE, we will be looking at Alternative 3' and 5' (A3SS/A5SS) Splicing, Retatined Intron (RI) Events, and Skipped Exon (SE) Events. rMATS has many output files. GrASE will use JCEC.*.txt files, one for each of the previously listed AS events. Each of those contain a matrix that have features including coordinates, pvalue, FDR, inclusion and exclusion levels, and more for every AS event that meets rMATS' filtering criteria. GrASE will use the output of rMATS in order to map each event to a graphML object. This will be compared to DEXSeq results in order to map DEXSeq exons to rMATS events.

## Creating graphMLs
GraphML is an XML-based file format for graphs. This step ensures that the coordinates of each gene are annotated and will be essential when creating the splicing graphs. It uses the [SplicingGraphs](https://bioconductor.org/packages/release/bioc/html/SplicingGraphs.html) R package to obtain the splice junctions and exons. 

To generate the GraphML objects for each gene, run [SplicingGraphs.igraph.R](SplicingGraphs.igraph.R). Some command line arguments will be required when you run the R script, including the path to the gtf, the name of the organism associated with the gtf (include genus and species), and the path to the output directory.
```
Rscript SplicingGraphs.igraph.R /path/to/gtf Genus species /path/to/output_directory
```
## Preparing to run GrASE
GrASE will process every gene in your dataset that produces results in DEXSeq and rMATS. In order to properly run GrASE, some setup needs to be done. Run [creatingFilesByGene.sh](creatingFilesByGene.sh) to create and set up your `grase_results` directory (created in your current working directory), which will hold everything you need to run GrASE. 
```
bash creatingFilesByGene.sh -r /path/to/rmats/results -d /path/to/dexseq_prepare_annotation.py -a /path/to/annotation/file.gtf -g /path/to/graphml/directory -p number_of_threads
```

This script will create the `grase results` directory, which will contain: 
* `gene_files`, a directory for each relevant gene to be processed in your dataset (taken from rMATS and DEXSeq results)
* `results`, a directory that will hold the final results after running grase.py (next step)

Here is an example of a toy file structure:
```bash
grase_results
├── gene_files
│   ├── ENSG00000286751.1
│   │   ├── ENSG00000286751.1.dexseq.gff
│   │   ├── ENSG00000286751.1.graphml
│   │   ├── ENSG00000286751.1.gtf
│   │   ├── fromGTF.A3SS.txt
│   │   ├── fromGTF.A5SS.txt
│   │   ├── fromGTF.RI.txt
│   │   ├── fromGTF.SE.txt
│   │   └── output
│   ├── ENSG00000287688.1
│   │   ├── ENSG00000287688.1.dexseq.gff
│   │   ├── ENSG00000287688.1.graphml
│   │   ├── ENSG00000287688.1.gtf
│   │   ├── fromGTF.A3SS.txt
│   │   ├── fromGTF.A5SS.txt
│   │   ├── fromGTF.RI.txt
│   │   ├── fromGTF.SE.txt
│   │   └── output
│   └── ENSG00000288597.1
│       ├── ENSG00000288597.1.dexseq.gff
│       ├── ENSG00000288597.1.graphml
│       ├── ENSG00000288597.1.gtf
│       ├── fromGTF.A3SS.txt
│       ├── fromGTF.A5SS.txt
│       ├── fromGTF.RI.txt
│       ├── fromGTF.SE.txt
│       └── output
└── results
    ├── ExonParts
    ├── SplicingEvents
    └── tmp
```

## Running GrASE
Now that the file structure is prepared, grase.py is ready to be ran. This is where the DEXSeq exon and rMATS event mapping and processing will be done. Each gene in gene_files will be processed, and the output directory of each will be populated with results tables and a graph png of the gene's structure after mapping. 
  
usage:
```
python3 grase.py -g <gene_files> --rmats <rmats_results_directory> --dexseq <dexseq_results.txt> --nthread <number_of_threads>
```

Example of an output graph for a simple gene: 
![alt text](docs/graph.ENSG00000287356.1.png "Logo Title Text 1")
    
## Usage
### All Arguments
```
usage: Rscript SplicingGraphs.igraph.R [options]

options:
 GTF                               An annotation of genes and transcripts in GTF format
 Genus
 species
 Output Directory                  The directory where all of the created graphML files
                                   will be placed

usage: bash creatingFilesByGene [options]

options:
 -r rMATS_results
 -d dexseq_prepare_annotation.py
 -a GTF
 -g graphML Directory
 -p NPROCS

usage: python grase.py [options]

options:
 -h, --help                         Display a help message and exit
 -g Gene Files Directory            The gene_files directory inside grase_results that
                                    was created by creatingFilesByGene.sh
 --rmats rMATS Results Directory    The OD directory that holds the final output of the
                                    post step of rMATS
 --dexseq Dexseq Results File       The file that holds results from DEXSeq
 --nthread NTHREAD                  The number of threads. The optimal number of threads
                                    should be equal to the number of cpu cores. Default: 1
```

## Final Output
`grase_results/results` contains the final output files from GrASE
*  `summary.txt`:
*  `DEX_to_rMATS_Events.txt`:
*  `rMATS_to_DEX_Exons.txt`:
*  `Mapped.ExonsToEvents.txt`:
*  `Mapped.EventsToExons.txt`:

`grase_results/results/ExonParts` contains the output files that informed our Exon counts in `summary.txt`
*  `DexSigExons.txt`:
*  `rMATS_SigExons.txt`:
*  `rMATS_TestedExons.txt`:
*  `rMATS_Tested__DexSigExons.txt`:
*  `rMATS_Sig__DexSigExons.txt`:

`grase_results/results/SplicingEvents` contains the output files that informed our Events counts in `summary.txt`
*  `DexSigEvents.txt`:
*  `rMATS_SigEvents.txt`:
*  `rMATS_Sig__DexSigEvents.txt`:

