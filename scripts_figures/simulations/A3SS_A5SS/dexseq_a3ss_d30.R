library(DEXSeq)
library(foreach)
library(doParallel)
library(tidyverse)
library(BiocParallel)

# cp dexseq annotation file to all simulation folders
# cat all_folders_50K.txt | xargs -n 1 -P 100 -I {} cp ENSG****.gff simulations/{}/ENSG****.gff

# create a table of one column and a row for each of the directories in the simulations folder(s) (names stored in folders_50K.txt) 
# edit files so that ENSG**** is replaced with FC#_RD#_#.ENSG****
# combine all gff files in each simulation folder into one large concatenated file: combined_ENSG****.gff
# combine all *_dexseq.cnts.txt files of sample[i] for each simulation folder into one large concatenated file: sample_[i]_dexseq.cnts.txt


# begin dexseq with combined/merged files
#######################################################################################################################################################################


flattenedFile = list.files("/scratch/han_lab/dwito/A3SS_A5SS_D30/dexseq/A3SS", pattern="combined_ENSG00000124496.12.gff", full.names=TRUE) #read in dexseq annotation file
sampleTable = data.frame(
  row.names = c( "sample_01", "sample_02", "sample_03", "sample_04", "sample_05",
                 "sample_06", "sample_07", "sample_08", "sample_09", "sample_10"),
  condition = c( "control", "control", "control", "control", "control", 
                 "case_a3", "case_a3", "case_a3", "case_a3", "case_a3"))
countFiles = list.files("/scratch/han_lab/dwito/A3SS_A5SS_D30/dexseq/A3SS/samples", pattern="dexseq.cnts.txt$", full.names=TRUE)
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
BPPARAM = MulticoreParam(20)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM = BPPARAM )
#dxd = estimateDispersions( dxd, fitType='mean' )        ((keep commented as it was originally commented out))
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition" )
dxr1 = DEXSeqResults( dxd )
dxr1 <- as.data.frame(dxr1)
rownames(dxr1) <- NULL
dxr1$transcripts <- vapply(dxr1$transcripts, paste, collapse = ", ", character(1L))
dxr1 <- unite(dxr1, col='countData.control', c('countData.sample_01','countData.sample_02','countData.sample_03','countData.sample_04','countData.sample_05'), sep=',')
dxr1 <- unite(dxr1, col='countData.case_a3', c('countData.sample_06','countData.sample_07','countData.sample_08','countData.sample_09','countData.sample_10'), sep=',')
write.table(dxr1,"A3SS_dexseq_results_50K_sim.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

dexseq_successful_pvalues <- dxr1 %>% 
  filter(!is.na(pvalue) & is.na(padj))
dexseq_unsuccessful_pvalues <- dxr1 %>% 
  filter(is.na(pvalue) & is.na(padj))
dexseq_successful_both <- dxr1 %>% 
  filter(!is.na(pvalue) & !is.na(padj))

dexseq_successful <- rbind(dexseq_successful_both, dexseq_successful_pvalues)
dexseq_successful$corrected.padj <- p.adjust(dexseq_successful$pvalue, method="BH")
write.table(dexseq_successful,"A3SS_dexseq_corrected_padj.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)





# commands after dexseq runs to test successful runs
#############################################################3####################################

length(dexseq_successful$corrected.padj[dexseq_successful$corrected.padj<0.05])
length(dexseq_successful$padj[dexseq_successful$padj<0.05 & !is.na(dexseq_successful$padj)])
length(dexseq_successful$pvalue[dexseq_successful$pvalue<0.05 & !is.na(dexseq_successful$padj)])
