#library(ggsci)
library(RColorBrewer)
library(tidyverse)
library(tidyr)
library(gtools)

#DEXSeq simulation runs - combined A3SS
dxr1 <- read.table("/scratch/han_lab/dwito/A3SS_A5SS_D30/dexseq/A3SS/A3SS_dexseq_corrected_padj.txt", header=TRUE, sep="\t")

dxr1 <- dxr1 %>% 
  separate_wider_delim(groupID, ".", names = c("ID", "groupID", "extra"))
dxr1 <- unite(dxr1, col='groupID', c('groupID','extra'), sep='.')

dxr1_exon6 <- dxr1 %>% 
  filter(featureID=="E006")

dxr1_exon6 <- dxr1_exon6 %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
dxr1_exon6$ReadDepth <- as.numeric(gsub("RD","", dxr1_exon6$ReadDepth))
dxr1_exon6$FoldChange <- gsub("FC","", dxr1_exon6$FoldChange)
dxr1_exon6$FoldChange <- gsub("minus","-", dxr1_exon6$FoldChange)
dxr1_exon6$FoldChange <- gsub("plus","", dxr1_exon6$FoldChange)
dxr1_exon6$FoldChange <- as.numeric(dxr1_exon6$FoldChange)

dexseq_FC_RD_freq <- data.frame(table(dxr1_exon6$FoldChange, dxr1_exon6$ReadDepth))
dexseq_FC_RD_freq$ID <- "DEXSeq"

#rMATS simulation runs - A3SS JCEC & JC
jcec <- read.table("/scratch/han_lab/dwito/A3SS_A5SS_D30/rmats/A3SS/A3SS.MATS.JCEC.txt", header=TRUE, sep="\t")
jcec <- jcec %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
jcec$ReadDepth <- as.numeric(gsub("RD","", jcec$ReadDepth))
jcec$FoldChange <- gsub("FCminus","-", jcec$FoldChange)
jcec$FoldChange <- gsub("FCplus","", jcec$FoldChange)
jcec$FoldChange <- gsub("FC","", jcec$FoldChange)
jcec$FoldChange <- as.numeric(jcec$FoldChange)

rMATS_FC_RD_freq <- data.frame(table(jcec$FoldChange, jcec$ReadDepth))
rMATS_FC_RD_freq$ID <- "rMATS"

jc <- read.table("/scratch/han_lab/dwito/A3SS_A5SS_D30/rmats/A3SS/A3SS.MATS.JC.txt", header=TRUE, sep="\t")
jc <- jc %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
jc$ReadDepth <- as.numeric(gsub("RD","", jc$ReadDepth))
jc$FoldChange <- gsub("FCminus","-", jc$FoldChange)
jc$FoldChange <- gsub("FCplus","", jc$FoldChange)
jc$FoldChange <- gsub("FC","", jc$FoldChange)
jc$FoldChange <- as.numeric(jc$FoldChange)

#total simulations
folders_50K_A3SS <- read.table("/scratch/han_lab/dwito/A3SS_A5SS_D30/folders_50K_A3SS.txt")
total_sim <- folders_50K_A3SS %>% 
  separate_wider_delim(V1, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
total_sim$ReadDepth <- as.numeric(gsub("RD","", total_sim$ReadDepth))
total_sim$FoldChange <- gsub("FC","", total_sim$FoldChange)
total_sim$FoldChange <- gsub("minus","-", total_sim$FoldChange)
total_sim$FoldChange <- gsub("plus","", total_sim$FoldChange)
total_sim$FoldChange <- as.numeric(total_sim$FoldChange)
total_FC_RD_freq <- data.frame(table(total_sim$FoldChange, total_sim$ReadDepth))
total_FC_RD_freq$ID <- "Total"

#plot that shows simulation distributions for total, DEXSeq, and rMATS
combined_sim <- rbind(total_FC_RD_freq, rMATS_FC_RD_freq, dexseq_FC_RD_freq)
combined_sim$ID <- factor(combined_sim$ID, levels = c("Total", "DEXSeq", "rMATS"))
ggplot(combined_sim , aes(x=Var2, y=Freq, fill=ID)) + 
  geom_bar(stat="identity", position="identity", alpha=0.5) +
  facet_wrap(~Var1, scales="free_y") + 
  labs(title="Total Simulation Runs vs Successful DEXSeq and rMATS Runs", x="Read Depth", fill = "Simulations") +
  theme_bw() +
  scale_fill_manual(values=c("cyan3", "cornflowerblue", "darkorchid2")) 
ggsave("simulation_dist.pdf", width = 14, height = 8, dpi = 100, units = "in", device='pdf') 

########################################################################################################################
# dexseq exon 6 (spliced exon) table
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
#specificity and sensitivity for pvalue

for(i in unique(mixedsort(dexseq_jcec_merge$ReadDepth))){
  merged_filt <- filter(dexseq_jcec_merge, FoldChange == 0 & ReadDepth == i)
  print(i)
  print(table(merged_filt$sig.pvalue.y, merged_filt$sig.pvalue.x, dnn=c("rMATS", "DEXSeq")))
}

for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
  merged_filt <- filter(dexseq_jc_merge, FoldChange == 0 & ReadDepth == i)
  print(i)
  print(table(merged_filt$sig.pvalue.y, merged_filt$sig.pvalue.x, dnn=c("rMATS", "DEXSeq")))
}

#A3SS false positives and True Negatives
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

#plot specificity
ggplot(specificity_tbl, aes(group=var, y=val, x=ReadDepth)) +
  geom_line(aes(color=var), position = position_dodge(width = 0.1)) +
  labs(title="A3SS Specificity raw pvalue", x="Read Depth", y="Specificity") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1)+
  theme_bw() 
ggsave("Spec_pvalue_lineplot.pdf", width = 6, height = 5, dpi = 100, units = "in", device='pdf')

#A3SS false negatives and true positives
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

#plot sensitivity
ggplot(sens_tbl, aes(group=var, y=val, x=ReadDepth)) +
  geom_line(aes(color=var), position = position_dodge(width = 0.2)) +
  labs(title="A3SS Sensitivity raw pvalue", x="Read Depth", y="Sensitivity") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1)+
  theme_bw() +
  facet_wrap(~FoldChange, nrow=2)
ggsave("Sens_lineplot.pdf", width = 18, height = 7, dpi = 100, units = "in", device='pdf')

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

#plot sensitivity
ggplot(sens_tbl, aes(group=var, y=val, x=ReadDepth)) +
  geom_line(aes(color=var), position = position_dodge(width = 0.2)) +
  labs(title="A3SS Sensitivity adjusted pvalue", x="Read Depth", y="Sensitivity") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1)+
  theme_bw() +
  facet_wrap(~FoldChange, nrow=2)
ggsave("Sens_adjpvalue_lineplot.pdf", width = 18, height = 7, dpi = 100, units = "in", device='pdf')

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

#plot specificity
ggplot(specificity_tbl, aes(group=var, y=val, x=ReadDepth)) +
  geom_line(aes(color=var), position = position_dodge(width = 0.1)) +
  labs(title="A3SS Specificity adjusted pvalue", x="Read Depth", y="Specificity") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1)+
  theme_bw()
ggsave("Spec_adjpvalue_lineplot.pdf", width = 6, height = 5, dpi = 100, units = "in", device='pdf')
