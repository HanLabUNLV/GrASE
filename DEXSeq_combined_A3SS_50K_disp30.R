library(DEXSeq)
library(foreach)
library(doParallel)
library(BiocParallel)
library(tidyverse)

# cp dexseq annotation file to all simulation folders
#cat all_folders_50K.txt | xargs -n 1 -P 80 -I {} cp ENSG00000117222.14.gff simulations_disp3_SE/{}/ENSG00000117222.14.gff

my.cluster <- parallel::makeCluster(
  80, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster) 

folders_50K <- read.table("/mnt/storage/jaquino/polyester_sim_2/SE_A5SS/folders_ID_50K.txt", header=FALSE, sep="\t")

#foreach(i = 1:nrow(folders_50K)) %dopar% {
#  dexseq_counts <- list.files(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i]), pattern='_dexseq.cnts.txt$', full.names=TRUE)
#  for(j in dexseq_counts){
#    cnts <- read.table(j, header=FALSE, sep="\t")
#    cnts$V1 <- gsub("EN",paste0(folders_50K$V1[i], ".EN"), cnts$V1)
#    write.table(cnts, j, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#  }
#}


#foreach(i = 1:nrow(folders_50K)) %dopar% {
#  gff_file <- read.table(paste0('simulations_disp3_SE/simulated_reads_', folders_50K$V1[i],"/ENSG00000117222.14.gff"), header=FALSE, sep="\t")
#  gff_file$V9 <- gsub("ENSG",paste0(folders_50K$V1[i], ".ENSG"), gff_file$V9)
#  write.table(gff_file, paste0('simulations_disp3_SE/simulated_reads_', folders_50K$V1[i],"/ENSG00000117222.14.gff"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#}

#gff <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  gff_file <- read.table(paste0('simulations_disp3_SE/simulated_reads_', folders_50K$V1[i],"/ENSG00000117222.14.gff"), header=FALSE, sep="\t")
#  #gff_file$V9 <- gsub("ENSG",paste0(folders_50K$V1[i], ".ENSG"), gff_file$V9)
#  return(gff_file)
#}

#write.table(gff,"combined_ENSG00000117222.14.gff", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_01Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_01_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_02Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_02_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_03Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_03_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_04Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_04_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_05Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_05_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_06Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_06_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_07Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_07_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_08Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_08_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_09Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_09_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#
#cnts <- foreach(i = 1:nrow(folders_50K),.combine=rbind) %dopar% {
#  dexseq_counts <- read.table(paste0('simulations_disp30_SE/simulated_reads_', folders_50K$V1[i], "/sample_10Aligned.sortedByCoord.out.bam_dexseq.cnts.txt"), header=FALSE, sep="\t")
#  return(dexseq_counts)
#}
#
#write.table(cnts,"sample_10_dexseq.cnts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#parallel::stopCluster(cl = my.cluster)

flattenedFile = list.files("/mnt/storage/jaquino/polyester_sim_2/SE_A5SS", pattern="combined_ENSG*", full.names=TRUE) #read in dexseq annotation file
sampleTable = data.frame(
  row.names = c( "sample_01", "sample_02", "sample_03", "sample_04", "sample_05",
                 "sample_06", "sample_07", "sample_08", "sample_09", "sample_10"),
  condition = c(  "control", "control", "control", "control", "control", 
                  "case", "case", "case", "case", "case")) 
countFiles = list.files('/mnt/storage/jaquino/polyester_sim_2/SE_A5SS/dexseq_disp30', pattern="dexseq.cnts.txt$", full.names=TRUE)
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
BPPARAM = MulticoreParam(60)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM )
#dxd = estimateDispersions( dxd, fitType='mean' )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 = DEXSeqResults( dxd )
dxr1 <- as.data.frame(dxr1)
rownames(dxr1) <- NULL
dxr1$transcripts <- vapply(dxr1$transcripts, paste, collapse = ", ", character(1L))
dxr1 <- unite(dxr1, col='countData.control', c('countData.sample_01','countData.sample_02','countData.sample_03','countData.sample_04','countData.sample_05'), sep=',')
dxr1 <- unite(dxr1, col='countData.case', c('countData.sample_06','countData.sample_07','countData.sample_08','countData.sample_09','countData.sample_10'), sep=',')
write.table(dxr1,"dexseq_results_50K_sim.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

dexseq_successful_pvalues <- dxr1 %>% 
  filter(!is.na(pvalue) & is.na(padj))
dexseq_unsuccessful_pvalues <- dxr1 %>% 
  filter(is.na(pvalue) & is.na(padj))
dexseq_successful_both <- dxr1 %>% 
  filter(!is.na(pvalue) & !is.na(padj))

dexseq_successful <- rbind(dexseq_successful_both, dexseq_successful_pvalues)
dexseq_successful$corrected.padj <- p.adjust(dexseq_successful$pvalue, method="BH")
write.table(dexseq_successful,"dexseq_corrected_padj.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#length(dexseq_successful$corrected.padj[dexseq_successful$corrected.padj<0.05])
#length(dexseq_successful$padj[dexseq_successful$padj<0.05 & !is.na(dexseq_successful$padj)])
#length(dexseq_successful$pvalue[dexseq_successful$pvalue<0.05 & !is.na(dexseq_successful$padj)])

