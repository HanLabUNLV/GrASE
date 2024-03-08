library(ggsci)
library(RColorBrewer)
library(tidyverse)

folders_50K <- read.table("/mnt/storage/jaquino/polyester_sim_2/SE_A5SS/folders_ID_50K.txt", header=FALSE, sep="\t")

#A5SS from the RBBP5 gene
#Successful DEXSeq simulation runs
dxr1 <- read.table("/mnt/storage/jaquino/polyester_sim_2/SE_A5SS/dexseq_disp30/dexseq_corrected_padj.txt", header=TRUE, sep="\t")

dxr1 <- dxr1 %>% 
  separate_wider_delim(groupID, ".", names = c("ID", "groupID", "extra"))
dxr1 <- unite(dxr1, col='groupID', c('groupID','extra'), sep='.')

dxr1_exon3 <- dxr1 %>% 
  filter(featureID=="E002")

dxr1_exon3 <- dxr1_exon3 %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
dxr1_exon3$ReadDepth <- as.numeric(gsub("RD","", dxr1_exon3$ReadDepth))
dxr1_exon3$FoldChange <- gsub("FC","", dxr1_exon3$FoldChange)
dxr1_exon3$FoldChange <- gsub("minus","-", dxr1_exon3$FoldChange)
dxr1_exon3$FoldChange <- gsub("plus","", dxr1_exon3$FoldChange)
dxr1_exon3$FoldChange <- as.numeric(dxr1_exon3$FoldChange)

dexseq_FC_RD_freq <- data.frame(table(dxr1_exon3$FoldChange, dxr1_exon3$ReadDepth))
dexseq_FC_RD_freq$ID <- "DEXSeq"

jcec <- read.table("/mnt/storage/jaquino/polyester_sim_2/SE_A5SS/A5SS_disp30_rmats/A5SS.MATS.JCEC.txt", header=TRUE, sep="\t")
jcec <- jcec %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
jcec$ReadDepth <- as.numeric(gsub("RD","", jcec$ReadDepth))
jcec$FoldChange <- gsub("FCminus","-", jcec$FoldChange)
jcec$FoldChange <- gsub("FCplus","", jcec$FoldChange)
jcec$FoldChange <- gsub("FC","", jcec$FoldChange)
jcec$FoldChange <- as.numeric(jcec$FoldChange)

rMATS_FC_RD_freq <- data.frame(table(jcec$FoldChange, jcec$ReadDepth))
rMATS_FC_RD_freq$ID <- "rMATS"

jc <- read.table("/mnt/storage/jaquino/polyester_sim_2/SE_A5SS/A5SS_disp30_rmats/A5SS.MATS.JC.txt", header=TRUE, sep="\t")
jc <- jc %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
jc$ReadDepth <- as.numeric(gsub("RD","", jc$ReadDepth))
jc$FoldChange <- gsub("FCminus","-", jc$FoldChange)
jc$FoldChange <- gsub("FCplus","", jc$FoldChange)
jc$FoldChange <- gsub("FC","", jc$FoldChange)
jc$FoldChange <- as.numeric(jc$FoldChange)

########################################################################################################################
# dexseq exon 3 (skipped exon) table
dxr1_exon3$sig.pvalue <- NA
for(i in 1:nrow(dxr1_exon3)){
  if(is.na(dxr1_exon3$pvalue[i])){
    next
  }
  if(dxr1_exon3$pvalue[i] <= 0.05){
    dxr1_exon3$sig.pvalue[i] <- "+"
  }
  if(dxr1_exon3$pvalue[i] > 0.05){
    dxr1_exon3$sig.pvalue[i] <- "-"
  }
}

dxr1_exon3$sig.adjpvalue <- NA
for(i in 1:nrow(dxr1_exon3)){
  if(is.na(dxr1_exon3$corrected.padj[i])){
    next
  }
  if(dxr1_exon3$corrected.padj[i] <= 0.05){
    dxr1_exon3$sig.adjpvalue[i] <- "+"
  }
  if(dxr1_exon3$corrected.padj[i] > 0.05){
    dxr1_exon3$sig.adjpvalue[i] <- "-"
  }
}

# rmats jcec table

jcec$sig.pvalue <- NA
for(i in 1:nrow(jcec)){
  if(is.na(jcec$PValue[i])){
    next
  }
  if(jcec$PValue[i] <= 0.05){
    jcec$sig.pvalue[i] <- "+"
  }
  if(jcec$PValue[i] > 0.05){
    jcec$sig.pvalue[i] <- "-"
  }
}

jcec$sig.adjpvalue <- NA
for(i in 1:nrow(jcec)){
  if(is.na(jcec$FDR[i])){
    next
  }
  if(jcec$FDR[i] <= 0.05){
    jcec$sig.adjpvalue[i] <- "+"
  }
  if(jcec$FDR[i] > 0.05){
    jcec$sig.adjpvalue[i] <- "-"
  }
}

# rmats jc table

jc$sig.pvalue <- NA
for(i in 1:nrow(jc)){
  if(is.na(jc$PValue[i])){
    next
  }
  if(jc$PValue[i] <= 0.05){
    jc$sig.pvalue[i] <- "+"
  }
  if(jc$PValue[i] > 0.05){
    jc$sig.pvalue[i] <- "-"
  }
}

jc$sig.adjpvalue <- NA
for(i in 1:nrow(jc)){
  if(is.na(jc$FDR[i])){
    next
  }
  if(jc$FDR[i] <= 0.05){
    jc$sig.adjpvalue[i] <- "+"
  }
  if(jc$FDR[i] > 0.05){
    jc$sig.adjpvalue[i] <- "-"
  }
}

dexseq_jcec_merge <- merge(dxr1_exon3, jcec, by=c("FoldChange", "ReadDepth", "SimRun"))
dexseq_jc_merge <- merge(dxr1_exon3, jc, by=c("FoldChange", "ReadDepth", "SimRun"))


###########################################################################################################################
#SE false positives and True Negatives
#specificity
DEXSeqTN=c()
DEXSeqFP=c()
rMATS_jcecTN=c()
rMATS_jcecFP=c()
rMATS_jcTN=c()
rMATS_jcFP=c()

for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
  print(i)
  merged_filt <- filter(dexseq_jc_merge, FoldChange == 0 & ReadDepth == i)
  
  dexseq_2way_tbl <- table(merged_filt$sig.pvalue.x)
  if(length(dexseq_2way_tbl)==2){
    DEXSeqFP <- c(DEXSeqFP, dexseq_2way_tbl[1])
    DEXSeqTN <- c(DEXSeqTN, dexseq_2way_tbl[2])
  }
  
  if(length(dexseq_2way_tbl)==1){
    if(names(dexseq_2way_tbl)[1]=="+"){
      DEXSeqFP<- c(DEXSeqFP, dexseq_2way_tbl[1])
      DEXSeqTN<- c(DEXSeqTN, 0)
    }
    if(names(dexseq_2way_tbl)[1]=="-"){
      DEXSeqTN<- c(DEXSeqTN, dexseq_2way_tbl[1])
      DEXSeqFP<- c(DEXSeqFP, 0)
    }
  }
  
  rmatsjc_2way_tbl <- table(merged_filt$sig.pvalue.y)
  if(length(rmatsjc_2way_tbl)==2){
    rMATS_jcFP <- c(rMATS_jcFP, rmatsjc_2way_tbl[1])
    rMATS_jcTN <- c(rMATS_jcTN, rmatsjc_2way_tbl[2])
  }
  
  if(length(rmatsjc_2way_tbl)==1){
    if(names(rmatsjc_2way_tbl)[1]=="+"){
      rMATS_jcFP<- c(rMATS_jcFP, rmatsjc_2way_tbl[1])
      rMATS_jcTN<- c(rMATS_jcTN, 0)
    }
    if(names(rmatsjc_2way_tbl)[1]=="-"){
      rMATS_jcTN<- c(rMATS_jcTN, rmatsjc_2way_tbl[1])
      rMATS_jcFP<- c(rMATS_jcFP, 0)
    }
  }
  
  merged_filt <- filter(dexseq_jcec_merge, FoldChange == 0 & ReadDepth == i)
  rmatsjcec_2way_tbl <- table(merged_filt$sig.pvalue.y)
  if(length(rmatsjcec_2way_tbl)==2){
    rMATS_jcecFP <- c(rMATS_jcecFP, rmatsjcec_2way_tbl[1])
    rMATS_jcecTN <- c(rMATS_jcecTN, rmatsjcec_2way_tbl[2])
  }
  
  if(length(rmatsjcec_2way_tbl)==1){
    if(names(rmatsjcec_2way_tbl)[1]=="+"){
      rMATS_jcecFP<- c(rMATS_jcecFP, rmatsjcec_2way_tbl[1])
      rMATS_jcecTN<- c(rMATS_jcecTN, 0)
    }
    if(names(rmatsjcec_2way_tbl)[1]=="-"){
      rMATS_jcecTN<- c(rMATS_jcecTN, rmatsjcec_2way_tbl[1])
      rMATS_jcecFP<- c(rMATS_jcecFP, 0)
    }
  }
}

FP_TN_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                       DEXSeqFP,
                       DEXSeqTN,
                       rMATS_jcecFP,
                       rMATS_jcecTN,
                       rMATS_jcFP,
                       rMATS_jcTN)


#get specificity
specificity_tbl <- FP_TN_df %>% 
  mutate(DEXSeqSpec = DEXSeqTN/(DEXSeqTN+DEXSeqFP),
         rMATS_jcecSpec = rMATS_jcecTN/(rMATS_jcecTN+rMATS_jcecFP),
         rMATS_jcSpec = rMATS_jcTN/(rMATS_jcTN+rMATS_jcFP))
specificity_tbl <- specificity_tbl %>% 
  gather(var, val, 8:ncol(specificity_tbl))
specificity_tbl$ReadDepth <- factor(specificity_tbl$ReadDepth, levels=c(1:10,50,100,500,1000))

specificity_tbl_SE_A3SS <- specificity_tbl
specificity_tbl_SE_A3SS$gene <- "RBBP5"

#SE false negatives and true positives
Foldchanges <- c(-8, -4, -2, -1, 8, 4, 2, 1)
sens_tbl <- data.frame(matrix(ncol = 8, nrow = 0))
for(FC in Foldchanges){
  DEXSeqFN=c()
  DEXSeqTP=c()
  rMATS_jcecFN=c()
  rMATS_jcecTP=c()
  rMATS_jcFN=c()
  rMATS_jcTP=c()
  for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
    print(i)
    merged_filt <- filter(dexseq_jc_merge, FoldChange == FC & ReadDepth == i)
    
    dexseq_2way_tbl <- table(merged_filt$sig.pvalue.x)
    if(length(dexseq_2way_tbl)==2){
      DEXSeqTP <- c(DEXSeqTP, dexseq_2way_tbl[1])
      DEXSeqFN <- c(DEXSeqFN, dexseq_2way_tbl[2])
    }
    
    if(length(dexseq_2way_tbl)==1){
      if(names(dexseq_2way_tbl)[1]=="+"){
        DEXSeqTP<- c(DEXSeqTP, dexseq_2way_tbl[1])
        DEXSeqFN<- c(DEXSeqFN, 0)
      }
      if(names(dexseq_2way_tbl)[1]=="-"){
        DEXSeqFN<- c(DEXSeqFN, dexseq_2way_tbl[1])
        DEXSeqTP<- c(DEXSeqTP, 0)
      }
    }
    
    rmatsjc_2way_tbl <- table(merged_filt$sig.pvalue.y)
    if(length(rmatsjc_2way_tbl)==2){
      rMATS_jcTP <- c(rMATS_jcTP, rmatsjc_2way_tbl[1])
      rMATS_jcFN <- c(rMATS_jcFN, rmatsjc_2way_tbl[2])
    }
    
    if(length(rmatsjc_2way_tbl)==1){
      if(names(rmatsjc_2way_tbl)[1]=="+"){
        rMATS_jcTP<- c(rMATS_jcTP, rmatsjc_2way_tbl[1])
        rMATS_jcFN<- c(rMATS_jcFN, 0)
      }
      if(names(rmatsjc_2way_tbl)[1]=="-"){
        rMATS_jcFN<- c(rMATS_jcFN, rmatsjc_2way_tbl[1])
        rMATS_jcTP<- c(rMATS_jcTP, 0)
      }
    }
    
    merged_filt <- filter(dexseq_jcec_merge, FoldChange == FC & ReadDepth == i)
    rmatsjcec_2way_tbl <- table(merged_filt$sig.pvalue.y)
    if(length(rmatsjcec_2way_tbl)==2){
      rMATS_jcecTP <- c(rMATS_jcecTP, rmatsjcec_2way_tbl[1])
      rMATS_jcecFN <- c(rMATS_jcecFN, rmatsjcec_2way_tbl[2])
    }
    
    if(length(rmatsjcec_2way_tbl)==1){
      if(names(rmatsjcec_2way_tbl)[1]=="+"){
        rMATS_jcecTP<- c(rMATS_jcecTP, rmatsjcec_2way_tbl[1])
        rMATS_jcecFN<- c(rMATS_jcecFN, 0)
      }
      if(names(rmatsjcec_2way_tbl)[1]=="-"){
        rMATS_jcecFN<- c(rMATS_jcecFN, rmatsjcec_2way_tbl[1])
        rMATS_jcecTP<- c(rMATS_jcecTP, 0)
      }
    }
  }
  
  
  FN_TP_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                         FoldChange = c(rep(FC, 14)),
                         DEXSeqFN,
                         DEXSeqTP,
                         rMATS_jcecFN,
                         rMATS_jcecTP,
                         rMATS_jcFN,
                         rMATS_jcTP)
  sens_tbl <- rbind(sens_tbl, FN_TP_df)
}

#get sensitivity
sens_tbl <- sens_tbl %>% 
  mutate(DEXSeqSens = DEXSeqTP/(DEXSeqTP+DEXSeqFN),
         rMATS_jcecSens = rMATS_jcecTP/(rMATS_jcecTP+rMATS_jcecFN),
         rMATS_jcSens = rMATS_jcTP/(rMATS_jcTP+rMATS_jcFN))
sens_tbl <- sens_tbl %>% 
  gather(var, val, 9:ncol(sens_tbl))
sens_tbl$ReadDepth <- factor(sens_tbl$ReadDepth, levels=c(1:10,50,100,500,1000))
sens_tbl$FoldChange <- factor(sens_tbl$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))

sens_tbl_SE_A3SS <- sens_tbl
sens_tbl_SE_A3SS$gene <- "RBBP5"

##################################################################################################################
# specificity and sensitivity for adjusted pvalues

sens_tbl <- data.frame(matrix(ncol = 8, nrow = 0))
for(FC in Foldchanges){
  DEXSeqFN=c()
  DEXSeqTP=c()
  rMATS_jcecFN=c()
  rMATS_jcecTP=c()
  rMATS_jcFN=c()
  rMATS_jcTP=c()
  for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
    print(i)
    merged_filt <- filter(dexseq_jc_merge, FoldChange == FC & ReadDepth == i)
    
    dexseq_2way_tbl <- table(merged_filt$sig.adjpvalue.x)
    if(length(dexseq_2way_tbl)==2){
      DEXSeqTP <- c(DEXSeqTP, dexseq_2way_tbl[1])
      DEXSeqFN <- c(DEXSeqFN, dexseq_2way_tbl[2])
    }
    
    if(length(dexseq_2way_tbl)==1){
      if(names(dexseq_2way_tbl)[1]=="+"){
        DEXSeqTP<- c(DEXSeqTP, dexseq_2way_tbl[1])
        DEXSeqFN<- c(DEXSeqFN, 0)
      }
      if(names(dexseq_2way_tbl)[1]=="-"){
        DEXSeqFN<- c(DEXSeqFN, dexseq_2way_tbl[1])
        DEXSeqTP<- c(DEXSeqTP, 0)
      }
    }
    
    rmatsjc_2way_tbl <- table(merged_filt$sig.adjpvalue.y)
    if(length(rmatsjc_2way_tbl)==2){
      rMATS_jcTP <- c(rMATS_jcTP, rmatsjc_2way_tbl[1])
      rMATS_jcFN <- c(rMATS_jcFN, rmatsjc_2way_tbl[2])
    }
    
    if(length(rmatsjc_2way_tbl)==1){
      if(names(rmatsjc_2way_tbl)[1]=="+"){
        rMATS_jcTP<- c(rMATS_jcTP, rmatsjc_2way_tbl[1])
        rMATS_jcFN<- c(rMATS_jcFN, 0)
      }
      if(names(rmatsjc_2way_tbl)[1]=="-"){
        rMATS_jcFN<- c(rMATS_jcFN, rmatsjc_2way_tbl[1])
        rMATS_jcTP<- c(rMATS_jcTP, 0)
      }
    }
    
    merged_filt <- filter(dexseq_jcec_merge, FoldChange == FC & ReadDepth == i)
    rmatsjcec_2way_tbl <- table(merged_filt$sig.adjpvalue.y)
    if(length(rmatsjcec_2way_tbl)==2){
      rMATS_jcecTP <- c(rMATS_jcecTP, rmatsjcec_2way_tbl[1])
      rMATS_jcecFN <- c(rMATS_jcecFN, rmatsjcec_2way_tbl[2])
    }
    
    if(length(rmatsjcec_2way_tbl)==1){
      if(names(rmatsjcec_2way_tbl)[1]=="+"){
        rMATS_jcecTP<- c(rMATS_jcecTP, rmatsjcec_2way_tbl[1])
        rMATS_jcecFN<- c(rMATS_jcecFN, 0)
      }
      if(names(rmatsjcec_2way_tbl)[1]=="-"){
        rMATS_jcecFN<- c(rMATS_jcecFN, rmatsjcec_2way_tbl[1])
        rMATS_jcecTP<- c(rMATS_jcecTP, 0)
      }
    }
  }
  
  
  FN_TP_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                         FoldChange = c(rep(FC, 14)),
                         DEXSeqFN,
                         DEXSeqTP,
                         rMATS_jcecFN,
                         rMATS_jcecTP,
                         rMATS_jcFN,
                         rMATS_jcTP)
  sens_tbl <- rbind(sens_tbl, FN_TP_df)
}

#get sensitivity
sens_tbl <- sens_tbl %>% 
  mutate(DEXSeqSens = DEXSeqTP/(DEXSeqTP+DEXSeqFN),
         rMATS_jcecSens = rMATS_jcecTP/(rMATS_jcecTP+rMATS_jcecFN),
         rMATS_jcSens = rMATS_jcTP/(rMATS_jcTP+rMATS_jcFN))
sens_tbl <- sens_tbl %>% 
  gather(var, val, 9:ncol(sens_tbl))
sens_tbl$ReadDepth <- factor(sens_tbl$ReadDepth, levels=c(1:10,50,100,500,1000))
sens_tbl$FoldChange <- factor(sens_tbl$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))

sens_tbl_SE_A3SS_adj <- sens_tbl
sens_tbl_SE_A3SS_adj$gene <- "RBBP5"

#specificity
DEXSeqTN=c()
DEXSeqFP=c()
rMATS_jcecTN=c()
rMATS_jcecFP=c()
rMATS_jcTN=c()
rMATS_jcFP=c()

for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
  print(i)
  merged_filt <- filter(dexseq_jc_merge, FoldChange == 0 & ReadDepth == i)
  
  dexseq_2way_tbl <- table(merged_filt$sig.adjpvalue.x)
  if(length(dexseq_2way_tbl)==2){
    DEXSeqFP <- c(DEXSeqFP, dexseq_2way_tbl[1])
    DEXSeqTN <- c(DEXSeqTN, dexseq_2way_tbl[2])
  }
  
  if(length(dexseq_2way_tbl)==1){
    if(names(dexseq_2way_tbl)[1]=="+"){
      DEXSeqFP<- c(DEXSeqFP, dexseq_2way_tbl[1])
      DEXSeqTN<- c(DEXSeqTN, 0)
    }
    if(names(dexseq_2way_tbl)[1]=="-"){
      DEXSeqTN<- c(DEXSeqTN, dexseq_2way_tbl[1])
      DEXSeqFP<- c(DEXSeqFP, 0)
    }
  }
  
  rmatsjc_2way_tbl <- table(merged_filt$sig.adjpvalue.y)
  if(length(rmatsjc_2way_tbl)==2){
    rMATS_jcFP <- c(rMATS_jcFP, rmatsjc_2way_tbl[1])
    rMATS_jcTN <- c(rMATS_jcTN, rmatsjc_2way_tbl[2])
  }
  
  if(length(rmatsjc_2way_tbl)==1){
    if(names(rmatsjc_2way_tbl)[1]=="+"){
      rMATS_jcFP<- c(rMATS_jcFP, rmatsjc_2way_tbl[1])
      rMATS_jcTN<- c(rMATS_jcTN, 0)
    }
    if(names(rmatsjc_2way_tbl)[1]=="-"){
      rMATS_jcTN<- c(rMATS_jcTN, rmatsjc_2way_tbl[1])
      rMATS_jcFP<- c(rMATS_jcFP, 0)
    }
  }
  
  merged_filt <- filter(dexseq_jcec_merge, FoldChange == 0 & ReadDepth == i)
  rmatsjcec_2way_tbl <- table(merged_filt$sig.adjpvalue.y)
  if(length(rmatsjcec_2way_tbl)==2){
    rMATS_jcecFP <- c(rMATS_jcecFP, rmatsjcec_2way_tbl[1])
    rMATS_jcecTN <- c(rMATS_jcecTN, rmatsjcec_2way_tbl[2])
  }
  
  if(length(rmatsjcec_2way_tbl)==1){
    if(names(rmatsjcec_2way_tbl)[1]=="+"){
      rMATS_jcecFP<- c(rMATS_jcecFP, rmatsjcec_2way_tbl[1])
      rMATS_jcecTN<- c(rMATS_jcecTN, 0)
    }
    if(names(rmatsjcec_2way_tbl)[1]=="-"){
      rMATS_jcecTN<- c(rMATS_jcecTN, rmatsjcec_2way_tbl[1])
      rMATS_jcecFP<- c(rMATS_jcecFP, 0)
    }
  }
}

FP_TN_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                       DEXSeqFP,
                       DEXSeqTN,
                       rMATS_jcecFP,
                       rMATS_jcecTN,
                       rMATS_jcFP,
                       rMATS_jcTN)


#get specificity
specificity_tbl <- FP_TN_df %>% 
  mutate(DEXSeqSpec = DEXSeqTN/(DEXSeqTN+DEXSeqFP),
         rMATS_jcecSpec = rMATS_jcecTN/(rMATS_jcecTN+rMATS_jcecFP),
         rMATS_jcSpec = rMATS_jcTN/(rMATS_jcTN+rMATS_jcFP))
specificity_tbl <- specificity_tbl %>% 
  gather(var, val, 8:ncol(specificity_tbl))
specificity_tbl$ReadDepth <- factor(specificity_tbl$ReadDepth, levels=c(1:10,50,100,500,1000))

specificity_tbl_SE_A3SS_adj <- specificity_tbl
specificity_tbl_SE_A3SS_adj$gene <- "RBBP5"
##################################################################################################################
# A5SS for TRERF1 gene (ENSG00000117400.18)
dxr1 <- read.table("/mnt/storage/jaquino/polyester_sim_2/A5SS_A3SS/disp30/A5SS_dexseq_corrected_padj.txt", header=TRUE, sep="\t")

dxr1 <- dxr1 %>% 
  separate_wider_delim(groupID, ".", names = c("ID", "groupID", "extra"))
dxr1 <- unite(dxr1, col='groupID', c('groupID','extra'), sep='.')

dxr1_exon6 <- dxr1 %>% 
  filter(featureID=="E015")

dxr1_exon6 <- dxr1_exon6 %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
dxr1_exon6$ReadDepth <- as.numeric(gsub("RD","", dxr1_exon6$ReadDepth))
dxr1_exon6$FoldChange <- gsub("FC","", dxr1_exon6$FoldChange)
dxr1_exon6$FoldChange <- gsub("minus","-", dxr1_exon6$FoldChange)
dxr1_exon6$FoldChange <- gsub("plus","", dxr1_exon6$FoldChange)
dxr1_exon6$FoldChange <- as.numeric(dxr1_exon6$FoldChange)

dexseq_FC_RD_freq <- data.frame(table(dxr1_exon6$FoldChange, dxr1_exon6$ReadDepth))
dexseq_FC_RD_freq$ID <- "DEXSeq"

#rMATS simulation runs - A5SS JCEC & JC
jcec <- read.table("/mnt/storage/jaquino/polyester_sim_2/A5SS_A3SS/disp30/A5SS.MATS.JCEC.txt", header=TRUE, sep="\t")
jcec <- jcec %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
jcec$ReadDepth <- as.numeric(gsub("RD","", jcec$ReadDepth))
jcec$FoldChange <- gsub("FCminus","-", jcec$FoldChange)
jcec$FoldChange <- gsub("FCplus","", jcec$FoldChange)
jcec$FoldChange <- gsub("FC","", jcec$FoldChange)
jcec$FoldChange <- as.numeric(jcec$FoldChange)

rMATS_FC_RD_freq <- data.frame(table(jcec$FoldChange, jcec$ReadDepth))
rMATS_FC_RD_freq$ID <- "rMATS"

jc <- read.table("/mnt/storage/jaquino/polyester_sim_2/A5SS_A3SS/disp30/A5SS.MATS.JCEC.txt", header=TRUE, sep="\t")
jc <- jc %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
jc$ReadDepth <- as.numeric(gsub("RD","", jc$ReadDepth))
jc$FoldChange <- gsub("FCminus","-", jc$FoldChange)
jc$FoldChange <- gsub("FCplus","", jc$FoldChange)
jc$FoldChange <- gsub("FC","", jc$FoldChange)
jc$FoldChange <- as.numeric(jc$FoldChange)

########################################################################################################################
# dexseq exon 6 (A3SS exon) table
dxr1_exon6$sig.pvalue <- NA
for(i in 1:nrow(dxr1_exon6)){
  if(is.na(dxr1_exon6$pvalue[i])){
    next
  }
  if(dxr1_exon6$pvalue[i] <= 0.05){
    dxr1_exon6$sig.pvalue[i] <- "+"
  }
  if(dxr1_exon6$pvalue[i] > 0.05){
    dxr1_exon6$sig.pvalue[i] <- "-"
  }
}

dxr1_exon6$sig.adjpvalue <- NA
for(i in 1:nrow(dxr1_exon6)){
  if(is.na(dxr1_exon6$corrected.padj[i])){
    next
  }
  if(dxr1_exon6$corrected.padj[i] <= 0.05){
    dxr1_exon6$sig.adjpvalue[i] <- "+"
  }
  if(dxr1_exon6$corrected.padj[i] > 0.05){
    dxr1_exon6$sig.adjpvalue[i] <- "-"
  }
}

# rmats jcec table

jcec$sig.pvalue <- NA
for(i in 1:nrow(jcec)){
  if(is.na(jcec$PValue[i])){
    next
  }
  if(jcec$PValue[i] <= 0.05){
    jcec$sig.pvalue[i] <- "+"
  }
  if(jcec$PValue[i] > 0.05){
    jcec$sig.pvalue[i] <- "-"
  }
}

jcec$sig.adjpvalue <- NA
for(i in 1:nrow(jcec)){
  if(is.na(jcec$FDR[i])){
    next
  }
  if(jcec$FDR[i] <= 0.05){
    jcec$sig.adjpvalue[i] <- "+"
  }
  if(jcec$FDR[i] > 0.05){
    jcec$sig.adjpvalue[i] <- "-"
  }
}

# rmats jc table

jc$sig.pvalue <- NA
for(i in 1:nrow(jc)){
  if(is.na(jc$PValue[i])){
    next
  }
  if(jc$PValue[i] <= 0.05){
    jc$sig.pvalue[i] <- "+"
  }
  if(jc$PValue[i] > 0.05){
    jc$sig.pvalue[i] <- "-"
  }
}

jc$sig.adjpvalue <- NA
for(i in 1:nrow(jc)){
  if(is.na(jc$FDR[i])){
    next
  }
  if(jc$FDR[i] <= 0.05){
    jc$sig.adjpvalue[i] <- "+"
  }
  if(jc$FDR[i] > 0.05){
    jc$sig.adjpvalue[i] <- "-"
  }
}

dexseq_jcec_merge <- merge(dxr1_exon6, jcec, by=c("FoldChange", "ReadDepth", "SimRun"))
dexseq_jc_merge <- merge(dxr1_exon6, jc, by=c("FoldChange", "ReadDepth", "SimRun"))

###########################################################################################################################
#SE false positives and True Negatives
#specificity
DEXSeqTN=c()
DEXSeqFP=c()
rMATS_jcecTN=c()
rMATS_jcecFP=c()
rMATS_jcTN=c()
rMATS_jcFP=c()

for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
  print(i)
  merged_filt <- filter(dexseq_jc_merge, FoldChange == 0 & ReadDepth == i)
  
  dexseq_2way_tbl <- table(merged_filt$sig.pvalue.x)
  if(length(dexseq_2way_tbl)==2){
    DEXSeqFP <- c(DEXSeqFP, dexseq_2way_tbl[1])
    DEXSeqTN <- c(DEXSeqTN, dexseq_2way_tbl[2])
  }
  
  if(length(dexseq_2way_tbl)==1){
    if(names(dexseq_2way_tbl)[1]=="+"){
      DEXSeqFP<- c(DEXSeqFP, dexseq_2way_tbl[1])
      DEXSeqTN<- c(DEXSeqTN, 0)
    }
    if(names(dexseq_2way_tbl)[1]=="-"){
      DEXSeqTN<- c(DEXSeqTN, dexseq_2way_tbl[1])
      DEXSeqFP<- c(DEXSeqFP, 0)
    }
  }
  
  rmatsjc_2way_tbl <- table(merged_filt$sig.pvalue.y)
  if(length(rmatsjc_2way_tbl)==2){
    rMATS_jcFP <- c(rMATS_jcFP, rmatsjc_2way_tbl[1])
    rMATS_jcTN <- c(rMATS_jcTN, rmatsjc_2way_tbl[2])
  }
  
  if(length(rmatsjc_2way_tbl)==1){
    if(names(rmatsjc_2way_tbl)[1]=="+"){
      rMATS_jcFP<- c(rMATS_jcFP, rmatsjc_2way_tbl[1])
      rMATS_jcTN<- c(rMATS_jcTN, 0)
    }
    if(names(rmatsjc_2way_tbl)[1]=="-"){
      rMATS_jcTN<- c(rMATS_jcTN, rmatsjc_2way_tbl[1])
      rMATS_jcFP<- c(rMATS_jcFP, 0)
    }
  }
  
  merged_filt <- filter(dexseq_jcec_merge, FoldChange == 0 & ReadDepth == i)
  rmatsjcec_2way_tbl <- table(merged_filt$sig.pvalue.y)
  if(length(rmatsjcec_2way_tbl)==2){
    rMATS_jcecFP <- c(rMATS_jcecFP, rmatsjcec_2way_tbl[1])
    rMATS_jcecTN <- c(rMATS_jcecTN, rmatsjcec_2way_tbl[2])
  }
  
  if(length(rmatsjcec_2way_tbl)==1){
    if(names(rmatsjcec_2way_tbl)[1]=="+"){
      rMATS_jcecFP<- c(rMATS_jcecFP, rmatsjcec_2way_tbl[1])
      rMATS_jcecTN<- c(rMATS_jcecTN, 0)
    }
    if(names(rmatsjcec_2way_tbl)[1]=="-"){
      rMATS_jcecTN<- c(rMATS_jcecTN, rmatsjcec_2way_tbl[1])
      rMATS_jcecFP<- c(rMATS_jcecFP, 0)
    }
  }
}

FP_TN_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                       DEXSeqFP,
                       DEXSeqTN,
                       rMATS_jcecFP,
                       rMATS_jcecTN,
                       rMATS_jcFP,
                       rMATS_jcTN)


#get specificity
specificity_tbl <- FP_TN_df %>% 
  mutate(DEXSeqSpec = DEXSeqTN/(DEXSeqTN+DEXSeqFP),
         rMATS_jcecSpec = rMATS_jcecTN/(rMATS_jcecTN+rMATS_jcecFP),
         rMATS_jcSpec = rMATS_jcTN/(rMATS_jcTN+rMATS_jcFP))
specificity_tbl <- specificity_tbl %>% 
  gather(var, val, 8:ncol(specificity_tbl))
specificity_tbl$ReadDepth <- factor(specificity_tbl$ReadDepth, levels=c(1:10,50,100,500,1000))

specificity_tbl_A3SS_A5SS <- specificity_tbl
specificity_tbl_A3SS_A5SS$gene <- "TRERF1"

#SE false negatives and true positives
Foldchanges <- c(-8, -4, -2, -1, 8, 4, 2, 1)
sens_tbl <- data.frame(matrix(ncol = 8, nrow = 0))
for(FC in Foldchanges){
  DEXSeqFN=c()
  DEXSeqTP=c()
  rMATS_jcecFN=c()
  rMATS_jcecTP=c()
  rMATS_jcFN=c()
  rMATS_jcTP=c()
  for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
    print(i)
    merged_filt <- filter(dexseq_jc_merge, FoldChange == FC & ReadDepth == i)
    
    dexseq_2way_tbl <- table(merged_filt$sig.pvalue.x)
    if(length(dexseq_2way_tbl)==2){
      DEXSeqTP <- c(DEXSeqTP, dexseq_2way_tbl[1])
      DEXSeqFN <- c(DEXSeqFN, dexseq_2way_tbl[2])
    }
    
    if(length(dexseq_2way_tbl)==1){
      if(names(dexseq_2way_tbl)[1]=="+"){
        DEXSeqTP<- c(DEXSeqTP, dexseq_2way_tbl[1])
        DEXSeqFN<- c(DEXSeqFN, 0)
      }
      if(names(dexseq_2way_tbl)[1]=="-"){
        DEXSeqFN<- c(DEXSeqFN, dexseq_2way_tbl[1])
        DEXSeqTP<- c(DEXSeqTP, 0)
      }
    }
    
    rmatsjc_2way_tbl <- table(merged_filt$sig.pvalue.y)
    if(length(rmatsjc_2way_tbl)==2){
      rMATS_jcTP <- c(rMATS_jcTP, rmatsjc_2way_tbl[1])
      rMATS_jcFN <- c(rMATS_jcFN, rmatsjc_2way_tbl[2])
    }
    
    if(length(rmatsjc_2way_tbl)==1){
      if(names(rmatsjc_2way_tbl)[1]=="+"){
        rMATS_jcTP<- c(rMATS_jcTP, rmatsjc_2way_tbl[1])
        rMATS_jcFN<- c(rMATS_jcFN, 0)
      }
      if(names(rmatsjc_2way_tbl)[1]=="-"){
        rMATS_jcFN<- c(rMATS_jcFN, rmatsjc_2way_tbl[1])
        rMATS_jcTP<- c(rMATS_jcTP, 0)
      }
    }
    
    merged_filt <- filter(dexseq_jcec_merge, FoldChange == FC & ReadDepth == i)
    rmatsjcec_2way_tbl <- table(merged_filt$sig.pvalue.y)
    if(length(rmatsjcec_2way_tbl)==2){
      rMATS_jcecTP <- c(rMATS_jcecTP, rmatsjcec_2way_tbl[1])
      rMATS_jcecFN <- c(rMATS_jcecFN, rmatsjcec_2way_tbl[2])
    }
    
    if(length(rmatsjcec_2way_tbl)==1){
      if(names(rmatsjcec_2way_tbl)[1]=="+"){
        rMATS_jcecTP<- c(rMATS_jcecTP, rmatsjcec_2way_tbl[1])
        rMATS_jcecFN<- c(rMATS_jcecFN, 0)
      }
      if(names(rmatsjcec_2way_tbl)[1]=="-"){
        rMATS_jcecFN<- c(rMATS_jcecFN, rmatsjcec_2way_tbl[1])
        rMATS_jcecTP<- c(rMATS_jcecTP, 0)
      }
    }
  }
  
  
  FN_TP_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                         FoldChange = c(rep(FC, 14)),
                         DEXSeqFN,
                         DEXSeqTP,
                         rMATS_jcecFN,
                         rMATS_jcecTP,
                         rMATS_jcFN,
                         rMATS_jcTP)
  sens_tbl <- rbind(sens_tbl, FN_TP_df)
}

#get sensitivity
sens_tbl <- sens_tbl %>% 
  mutate(DEXSeqSens = DEXSeqTP/(DEXSeqTP+DEXSeqFN),
         rMATS_jcecSens = rMATS_jcecTP/(rMATS_jcecTP+rMATS_jcecFN),
         rMATS_jcSens = rMATS_jcTP/(rMATS_jcTP+rMATS_jcFN))
sens_tbl <- sens_tbl %>% 
  gather(var, val, 9:ncol(sens_tbl))
sens_tbl$ReadDepth <- factor(sens_tbl$ReadDepth, levels=c(1:10,50,100,500,1000))
sens_tbl$FoldChange <- factor(sens_tbl$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))

sens_tbl_A3SS_A5SS <- sens_tbl
sens_tbl_A3SS_A5SS$gene <- "TRERF1"

#######################################################################################################################
# specificity and sensitivity for adjusted pvalues

sens_tbl <- data.frame(matrix(ncol = 8, nrow = 0))
for(FC in Foldchanges){
  DEXSeqFN=c()
  DEXSeqTP=c()
  rMATS_jcecFN=c()
  rMATS_jcecTP=c()
  rMATS_jcFN=c()
  rMATS_jcTP=c()
  for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
    print(i)
    merged_filt <- filter(dexseq_jc_merge, FoldChange == FC & ReadDepth == i)
    
    dexseq_2way_tbl <- table(merged_filt$sig.adjpvalue.x)
    if(length(dexseq_2way_tbl)==2){
      DEXSeqTP <- c(DEXSeqTP, dexseq_2way_tbl[1])
      DEXSeqFN <- c(DEXSeqFN, dexseq_2way_tbl[2])
    }
    
    if(length(dexseq_2way_tbl)==1){
      if(names(dexseq_2way_tbl)[1]=="+"){
        DEXSeqTP<- c(DEXSeqTP, dexseq_2way_tbl[1])
        DEXSeqFN<- c(DEXSeqFN, 0)
      }
      if(names(dexseq_2way_tbl)[1]=="-"){
        DEXSeqFN<- c(DEXSeqFN, dexseq_2way_tbl[1])
        DEXSeqTP<- c(DEXSeqTP, 0)
      }
    }
    
    rmatsjc_2way_tbl <- table(merged_filt$sig.adjpvalue.y)
    if(length(rmatsjc_2way_tbl)==2){
      rMATS_jcTP <- c(rMATS_jcTP, rmatsjc_2way_tbl[1])
      rMATS_jcFN <- c(rMATS_jcFN, rmatsjc_2way_tbl[2])
    }
    
    if(length(rmatsjc_2way_tbl)==1){
      if(names(rmatsjc_2way_tbl)[1]=="+"){
        rMATS_jcTP<- c(rMATS_jcTP, rmatsjc_2way_tbl[1])
        rMATS_jcFN<- c(rMATS_jcFN, 0)
      }
      if(names(rmatsjc_2way_tbl)[1]=="-"){
        rMATS_jcFN<- c(rMATS_jcFN, rmatsjc_2way_tbl[1])
        rMATS_jcTP<- c(rMATS_jcTP, 0)
      }
    }
    
    merged_filt <- filter(dexseq_jcec_merge, FoldChange == FC & ReadDepth == i)
    rmatsjcec_2way_tbl <- table(merged_filt$sig.adjpvalue.y)
    if(length(rmatsjcec_2way_tbl)==2){
      rMATS_jcecTP <- c(rMATS_jcecTP, rmatsjcec_2way_tbl[1])
      rMATS_jcecFN <- c(rMATS_jcecFN, rmatsjcec_2way_tbl[2])
    }
    
    if(length(rmatsjcec_2way_tbl)==1){
      if(names(rmatsjcec_2way_tbl)[1]=="+"){
        rMATS_jcecTP<- c(rMATS_jcecTP, rmatsjcec_2way_tbl[1])
        rMATS_jcecFN<- c(rMATS_jcecFN, 0)
      }
      if(names(rmatsjcec_2way_tbl)[1]=="-"){
        rMATS_jcecFN<- c(rMATS_jcecFN, rmatsjcec_2way_tbl[1])
        rMATS_jcecTP<- c(rMATS_jcecTP, 0)
      }
    }
  }
  
  
  FN_TP_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                         FoldChange = c(rep(FC, 14)),
                         DEXSeqFN,
                         DEXSeqTP,
                         rMATS_jcecFN,
                         rMATS_jcecTP,
                         rMATS_jcFN,
                         rMATS_jcTP)
  sens_tbl <- rbind(sens_tbl, FN_TP_df)
}

#get sensitivity
sens_tbl <- sens_tbl %>% 
  mutate(DEXSeqSens = DEXSeqTP/(DEXSeqTP+DEXSeqFN),
         rMATS_jcecSens = rMATS_jcecTP/(rMATS_jcecTP+rMATS_jcecFN),
         rMATS_jcSens = rMATS_jcTP/(rMATS_jcTP+rMATS_jcFN))
sens_tbl <- sens_tbl %>% 
  gather(var, val, 9:ncol(sens_tbl))
sens_tbl$ReadDepth <- factor(sens_tbl$ReadDepth, levels=c(1:10,50,100,500,1000))
sens_tbl$FoldChange <- factor(sens_tbl$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))

sens_tbl_A3SS_A5SS_adj <- sens_tbl
sens_tbl_A3SS_A5SS_adj$gene <- "TRERF1"

#specificity
DEXSeqTN=c()
DEXSeqFP=c()
rMATS_jcecTN=c()
rMATS_jcecFP=c()
rMATS_jcTN=c()
rMATS_jcFP=c()

for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
  print(i)
  merged_filt <- filter(dexseq_jc_merge, FoldChange == 0 & ReadDepth == i)
  
  dexseq_2way_tbl <- table(merged_filt$sig.adjpvalue.x)
  if(length(dexseq_2way_tbl)==2){
    DEXSeqFP <- c(DEXSeqFP, dexseq_2way_tbl[1])
    DEXSeqTN <- c(DEXSeqTN, dexseq_2way_tbl[2])
  }
  
  if(length(dexseq_2way_tbl)==1){
    if(names(dexseq_2way_tbl)[1]=="+"){
      DEXSeqFP<- c(DEXSeqFP, dexseq_2way_tbl[1])
      DEXSeqTN<- c(DEXSeqTN, 0)
    }
    if(names(dexseq_2way_tbl)[1]=="-"){
      DEXSeqTN<- c(DEXSeqTN, dexseq_2way_tbl[1])
      DEXSeqFP<- c(DEXSeqFP, 0)
    }
  }
  
  rmatsjc_2way_tbl <- table(merged_filt$sig.adjpvalue.y)
  if(length(rmatsjc_2way_tbl)==2){
    rMATS_jcFP <- c(rMATS_jcFP, rmatsjc_2way_tbl[1])
    rMATS_jcTN <- c(rMATS_jcTN, rmatsjc_2way_tbl[2])
  }
  
  if(length(rmatsjc_2way_tbl)==1){
    if(names(rmatsjc_2way_tbl)[1]=="+"){
      rMATS_jcFP<- c(rMATS_jcFP, rmatsjc_2way_tbl[1])
      rMATS_jcTN<- c(rMATS_jcTN, 0)
    }
    if(names(rmatsjc_2way_tbl)[1]=="-"){
      rMATS_jcTN<- c(rMATS_jcTN, rmatsjc_2way_tbl[1])
      rMATS_jcFP<- c(rMATS_jcFP, 0)
    }
  }
  
  merged_filt <- filter(dexseq_jcec_merge, FoldChange == 0 & ReadDepth == i)
  rmatsjcec_2way_tbl <- table(merged_filt$sig.adjpvalue.y)
  if(length(rmatsjcec_2way_tbl)==2){
    rMATS_jcecFP <- c(rMATS_jcecFP, rmatsjcec_2way_tbl[1])
    rMATS_jcecTN <- c(rMATS_jcecTN, rmatsjcec_2way_tbl[2])
  }
  
  if(length(rmatsjcec_2way_tbl)==1){
    if(names(rmatsjcec_2way_tbl)[1]=="+"){
      rMATS_jcecFP<- c(rMATS_jcecFP, rmatsjcec_2way_tbl[1])
      rMATS_jcecTN<- c(rMATS_jcecTN, 0)
    }
    if(names(rmatsjcec_2way_tbl)[1]=="-"){
      rMATS_jcecTN<- c(rMATS_jcecTN, rmatsjcec_2way_tbl[1])
      rMATS_jcecFP<- c(rMATS_jcecFP, 0)
    }
  }
}

FP_TN_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                       DEXSeqFP,
                       DEXSeqTN,
                       rMATS_jcecFP,
                       rMATS_jcecTN,
                       rMATS_jcFP,
                       rMATS_jcTN)


#get specificity
specificity_tbl <- FP_TN_df %>% 
  mutate(DEXSeqSpec = DEXSeqTN/(DEXSeqTN+DEXSeqFP),
         rMATS_jcecSpec = rMATS_jcecTN/(rMATS_jcecTN+rMATS_jcecFP),
         rMATS_jcSpec = rMATS_jcTN/(rMATS_jcTN+rMATS_jcFP))
specificity_tbl <- specificity_tbl %>% 
  gather(var, val, 8:ncol(specificity_tbl))
specificity_tbl$ReadDepth <- factor(specificity_tbl$ReadDepth, levels=c(1:10,50,100,500,1000))

specificity_tbl_A3SS_A5SS_adj <- specificity_tbl
specificity_tbl_A3SS_A5SS_adj$gene <- "TRERF1"

#combine the two specificity tables from RBBP5 and TRERF1 for A5SS for raw p-values
specificity_tbl <- rbind(specificity_tbl_SE_A3SS, specificity_tbl_A3SS_A5SS)
specificity_tbl$var2 <- paste0(specificity_tbl$var, specificity_tbl$gene)
write.table(specificity_tbl, "/mnt/storage/jaquino/polyester_sim_2/A5SS_RBBP5_TRERF1_res/Spec_pvalue_A5SS_RBBP5_TRERF1_disp30.txt", 
            sep="\t", quote=FALSE, row.names = FALSE)

#plot specificity for raw p-values
ggplot(specificity_tbl, aes(group=var2, y=val, x=ReadDepth)) +
  geom_line(aes(linetype=gene, color=var),size = 1, position = position_dodge(width = 0.1)) +
  labs(title="SE Specificity raw p-value", x="Read Depth", y="Specificity") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1)+
  theme_bw(base_size = 20)
ggsave("/mnt/storage/jaquino/polyester_sim_2/A5SS_RBBP5_TRERF1_res/Spec_pvalue_lineplot_A5SS_RBBP5_TRERF1_disp30.pdf", width = 10, height = 6, dpi = 100, units = "in", device='pdf')

#combine the two sensitivity tables from RBBP5 and TRERF1 for A5SS for raw p-values
sens_tbl <- rbind(sens_tbl_SE_A3SS, sens_tbl_A3SS_A5SS)
sens_tbl$var2 <- paste0(sens_tbl$var, sens_tbl$gene)
write.table(sens_tbl, "/mnt/storage/jaquino/polyester_sim_2/A5SS_RBBP5_TRERF1_res/Sens_pvalue_A5SS_RBBP5_TRERF1_disp30.txt", 
            sep="\t", quote=FALSE, row.names = FALSE)

#plot sensitivity for raw p-values
ggplot(sens_tbl, aes(group=var2, y=val, x=ReadDepth)) +
  geom_line(aes(linetype=gene, color=var), size =1, position = position_dodge(width = 0.2)) +
  labs(title="SE Sensitivity raw p-value", x="Read Depth", y="Sensitivity") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1) +
  theme_bw(base_size = 20) +
  facet_wrap(~FoldChange, nrow=2)
ggsave("/mnt/storage/jaquino/polyester_sim_2/A5SS_RBBP5_TRERF1_res/Sens_pvalue_lineplot_A5SS_RBBP5_TRERF1_disp30.pdf", width = 28, height = 12, dpi = 100, units = "in", device='pdf')

#combine the two specificity tables from RBBP5 and TRERF1 for A5SS for adjusted p-values
specificity_tbl <- rbind(specificity_tbl_SE_A3SS_adj, specificity_tbl_A3SS_A5SS_adj)
specificity_tbl$var2 <- paste0(specificity_tbl$var, specificity_tbl$gene)
write.table(specificity_tbl, "/mnt/storage/jaquino/polyester_sim_2/A5SS_RBBP5_TRERF1_res/Spec_pvalue_A5SS_RBBP5_TRERF1_disp30_adj.txt", 
            sep="\t", quote=FALSE, row.names = FALSE)

#plot specificity for adjusted p-values
ggplot(specificity_tbl, aes(group=var2, y=val, x=ReadDepth)) +
  geom_line(aes(linetype=gene, color=var),size = 1, position = position_dodge(width = 0.1)) +
  labs(title="SE Specificity adjusted p-value", x="Read Depth", y="Specificity") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1)+
  theme_bw(base_size = 20)
ggsave("/mnt/storage/jaquino/polyester_sim_2/A5SS_RBBP5_TRERF1_res/Spec_pvalue_lineplot_A5SS_RBBP5_TRERF1_adj_disp30.pdf", width = 10, height = 6, dpi = 100, units = "in", device='pdf')

#combine the two sensitivity tables from RBBP5 and TRERF1 for A5SS for adjusted p-values
sens_tbl <- rbind(sens_tbl_SE_A3SS_adj, sens_tbl_A3SS_A5SS_adj)
sens_tbl$var2 <- paste0(sens_tbl$var, sens_tbl$gene)
write.table(sens_tbl, "/mnt/storage/jaquino/polyester_sim_2/A5SS_RBBP5_TRERF1_res/Sens_pvalue_A5SS_RBBP5_TRERF1_disp30_adj.txt", 
            sep="\t", quote=FALSE, row.names = FALSE)

#plot sensitivity for adjusted p-values
ggplot(sens_tbl, aes(group=var2, y=val, x=ReadDepth)) +
  geom_line(aes(linetype=gene, color=var), size =1, position = position_dodge(width = 0.2)) +
  labs(title="SE Sensitivity adjusted p-value", x="Read Depth", y="Sensitivity") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1) +
  theme_bw(base_size = 20) +
  facet_wrap(~FoldChange, nrow=2)
ggsave("/mnt/storage/jaquino/polyester_sim_2/A5SS_RBBP5_TRERF1_res/Sens_pvalue_lineplot_A5SS_RBBP5_TRERF1_adj_disp30.pdf", width = 28, height = 12, dpi = 100, units = "in", device='pdf')

