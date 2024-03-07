library(DEXSeq)
library(tidyverse)

#filtered unique 112402 DEXSeq exons with mapped AS events that don't have any NA's
dexseq_with_AS_events <- read.table("/mnt/storage/jaquino/GrASE_results_bonnal/dexseq_res_all_filt_after.txt", header=TRUE, sep="\t")
#rMATS_tested exons with NA's and unique 138997 exon ID 
rMATS_TestedExons <- read.table("/mnt/storage/jaquino/GrASE_results_bonnal/ExonParts/rMATS_TestedExons.txt", header=TRUE, sep="\t")
rMATS_TestedExons$featureID <- gsub("E", "", rMATS_TestedExons$featureID)
rMATS_TestedExons <- rMATS_TestedExons %>% unite(col = "DEXSeqID", c(groupID, featureID), sep=":")

#non aggregated gff
flattenedFile = list.files("/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis",
                           pattern="gencode.v34.annotation.dexseq.gff",
                           full.names=TRUE) #read in dexseq annotation file
sampleTable = data.frame(
  row.names = c( "Bnaive1", "Bnaive2", "Bnaive3", "Bnaive4", "Bnaive5",
                "CD8naive1", "CD8naive2", "CD8naive3", "CD8naive4", "CD8naive5"),
  condition = c( "Bnaive", "Bnaive", "Bnaive", "Bnaive", "Bnaive",
                 "CD8naive", "CD8naive", "CD8naive", "CD8naive", "CD8naive"))
countFiles = list.files("/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis", pattern=".cnts.txt$", full.names=TRUE)

listOfFiles <- lapply(countFiles, function(x) read.table(x, header=FALSE, sep="\t"))
listOfFiles <- lapply(listOfFiles, function(z) z[z$V1 %in% rMATS_TestedExons$DEXSeqID,]) #keep dexseq exons that mapped to AS events
#listOfFiles <- lapply(listOfFiles, function(z) z[z$V1 %in% dexseq_with_AS_events$DEXSeqID,])

#save new counts files
sapply(1:length(listOfFiles), function(x) 
  write.table(listOfFiles[x], 
            countFiles[x], 
            row.names = FALSE, quote=FALSE, col.names = FALSE))

dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData = sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile = flattenedFile )

this_other <- as.data.frame(counts(dxd))
write.table(this_other, "this_others_counts.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

BPPARAM = MulticoreParam(8)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM )
#dxd = estimateDispersions( dxd, fitType='mean' )        ((keep commented as it was originally commented out))
dxd = testForDEU( dxd, BPPARAM=BPPARAM )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
dxr1 = DEXSeqResults( dxd )
dxr1 <- as.data.frame(dxr1)
rownames(dxr1) <- NULL
dxr1$transcripts <- vapply(dxr1$transcripts, paste, collapse = ", ", character(1L))
dxr1 <- unite(dxr1, col='countData.case_Bnaive', c('countData.Bnaive1','countData.Bnaive2','countData.Bnaive3', 'countData.Bnaive4', 'countData.Bnaive5'), sep=',')
dxr1 <- unite(dxr1, col='countData.case_CD8naive', c('countData.CD8naive1','countData.CD8naive2', 'countData.CD8naive3', 'countData.CD8naive4', 'countData.CD8naive5'), sep=',')
#write.table(dxr1,"/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/bonnal_dexseq_results_112402.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(dxr1,"/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/bonnal_dexseq_results_138997.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

dexseq_successful_pvalues <- dxr1 %>%
  filter(!is.na(pvalue) & is.na(padj))
dexseq_unsuccessful_pvalues <- dxr1 %>%
  filter(is.na(pvalue) & is.na(padj))
dexseq_successful_both <- dxr1 %>%
  filter(!is.na(pvalue) & !is.na(padj))

dexseq_successful <- rbind(dexseq_successful_both, dexseq_successful_pvalues)
dexseq_successful$corrected.padj <- p.adjust(dexseq_successful$pvalue, method="BH")
write.table(dexseq_successful,"bonnal_dexseq_corrected_padj.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
