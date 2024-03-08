# Scripts for generating figures in publication

## Splicing Graph script for a given gene

`splicegraph_tx_plot.R` is a script to make individual splice graph of a gene. Required input for this script is a **txdb** object of your reference genome. To make a txdb object using a gtf/gff of your reference, you can utilize the *makeTxDbFromGFF* function from *GenomicFeatures* library. In this script, you will have to modify the geneID (line 15) to the gene of interest and the chromosome location (line 9) of that gene. An example splicing graph plot of OR9HP1 is shown below. 
![](https://github.com/HanLabUNLV/GrASE/blob/main/scripts_figures/OR9HP1_splicegraph.png)

The yellow circles represent the splice junctions, the green lines represent exon edges, the purple lines are the exons for each transcript of the gene, and the gray dotted lines are the introns.

## Simulations
In the [simulations](https://github.com/HanLabUNLV/GrASE/tree/main/scripts_figures/simulations) folder, we have our scripts used to simulate differential isoform expression. [SE_A3SS]() will have scripts for the **MPL** gene, [SE_A5SS]() has scripts for **RBPP5**, and [A3SS_A5SS]() has scripts for **TRERF1**. To run the simulations, run `bash protocol.sh`. Required inputs are
1. transcript fasta file of gene (we used RSEM to obtain this)
2. gtf for that gene (grep your gene from the gtf file of your reference genome)
3. a txt file with all simulated fastq files (i.e. simulated_reads_FC0_RD1000_1/sample_01)
4. a txt file with all simulated BAM files (i.e simulated_reads_FC0_RD1000_1/sample_01Aligned.sortedByCoord.out.bam)

The directory [eval_metrics_combined](https://github.com/HanLabUNLV/GrASE/tree/main/scripts_figures/simulations/eval_meterics_combined) contains scripts to obtain sensitivity, specificity, and accuracy of DEXSeq and rMATS.

## Mapping of rMATS splicing events to DEXSeq on B naive and CD8 naive cells
After  running GrASE, we used `Mapped.ExonsToEvents.txt` and `Mapped.EventsToExons.txt` to generate heatmaps using `fig9_heatmap_bonnal.R`. These heatmaps will show the inclusion levels and the log2 fold change values of the this/others counts for each sample of significant events. 

To see the overlap between DEXSeq exons and rMATS events, we created venn diagrams using the output of *summary.txt* from running GrASE. `fig6_venn_diagram_exonicparts.R` is a script to create the Exonic Parts venn diagram.

[fig8_discrep_rMATS_DEXSeq](https://github.com/HanLabUNLV/GrASE/tree/main/scripts_figures/fig8_discrep_rMATS_DEXSeq) will have scripts for exploring the discrepancies between the two software in AS analysis. `this_others_barplots.R` is a script for generating Figures 8A and 8B. The required input is the `Mapped.ExonsToEvents.txt`. `zerocounts_barplots.R` is a script for generating Figures 8C and 8D. The required input is `Mapped.EventsToExons.txt`. 

[figure7_filteringDEXSeqExons](https://github.com/HanLabUNLV/GrASE/tree/main/scripts_figures/figure7_filteringDEXSeqExons) will have the scripts for rerunning DEXSeq after filtering for exons involved in at least one splicing event.  

## Supplemental Figures 2 and 3
`figS2_S3_featureCounts_readDepth.R` is a script to see the distributions of the log2FC of gene expression and read depth between B naive and CD8 naive cells.







