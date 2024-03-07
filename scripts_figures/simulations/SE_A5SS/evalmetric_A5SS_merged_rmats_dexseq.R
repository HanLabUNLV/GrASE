library(tidyverse)

setwd("/mnt/storage/jaquino/a5SS_a3SS_polyester_sim/a3ss")

#Successful DEXSeq simulation runs
dxr1 <- read.table("/mnt/storage/jaquino/a5SS_a3SS_polyester_sim/a3ss/dexseq_corrected_padj.txt", header=TRUE, sep="\t")

dxr1 <- dxr1 %>% 
  separate_wider_delim(groupID, ".", names = c("ID", "groupID", "extra"))
dxr1 <- unite(dxr1, col='groupID', c('groupID','extra'), sep='.')

dxr1_exon2 <- dxr1 %>% 
  filter(featureID=="E002")

dxr1_exon2 <- dxr1_exon2 %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
dxr1_exon2$ReadDepth <- as.numeric(gsub("RD","", dxr1_exon2$ReadDepth))
dxr1_exon2$FoldChange <- gsub("FC","", dxr1_exon2$FoldChange)
dxr1_exon2$FoldChange <- gsub("minus","-", dxr1_exon2$FoldChange)
dxr1_exon2$FoldChange <- gsub("plus","", dxr1_exon2$FoldChange)
dxr1_exon2$FoldChange <- as.numeric(dxr1_exon2$FoldChange)

jcec <- read.table("/mnt/storage/jaquino/a5SS_a3SS_polyester_sim/a3ss/A3SS.MATS.JCEC.txt", header=TRUE, sep="\t")
jcec <- jcec %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
jcec$ReadDepth <- as.numeric(gsub("RD","", jcec$ReadDepth))
jcec$FoldChange <- gsub("FCminus","-", jcec$FoldChange)
jcec$FoldChange <- gsub("FCplus","", jcec$FoldChange)
jcec$FoldChange <- gsub("FC","", jcec$FoldChange)
jcec$FoldChange <- as.numeric(jcec$FoldChange)

jc <- read.table("/mnt/storage/jaquino/a5SS_a3SS_polyester_sim/a3ss/A3SS.MATS.JC.txt", header=TRUE, sep="\t")
jc <- jc %>% 
  separate_wider_delim(ID, "_", names = c("FoldChange", "ReadDepth", "SimRun"))
jc$ReadDepth <- as.numeric(gsub("RD","", jc$ReadDepth))
jc$FoldChange <- gsub("FCminus","-", jc$FoldChange)
jc$FoldChange <- gsub("FCplus","", jc$FoldChange)
jc$FoldChange <- gsub("FC","", jc$FoldChange)
jc$FoldChange <- as.numeric(jc$FoldChange)

# dexseq exon 3 (skipped exon) table
dxr1_exon2$predicted.pvalue <- NA
for(i in 1:nrow(dxr1_exon2)){
  if(is.na(dxr1_exon2$pvalue[i])){
    next
  }
  if(dxr1_exon2$pvalue[i] <= 0.05){
    dxr1_exon2$predicted.pvalue[i] <- "+"
  }
  if(dxr1_exon2$pvalue[i] > 0.05){
    dxr1_exon2$predicted.pvalue[i] <- "-"
  }
}

dxr1_exon2$predicted.adjpvalue <- NA
for(i in 1:nrow(dxr1_exon2)){
  if(is.na(dxr1_exon2$corrected.padj[i])){
    next
  }
  if(dxr1_exon2$corrected.padj[i] <= 0.05){
    dxr1_exon2$predicted.adjpvalue[i] <- "+"
  }
  if(dxr1_exon2$corrected.padj[i] > 0.05){
    dxr1_exon2$predicted.adjpvalue[i] <- "-"
  }
}

# rmats jcec table

jcec$predicted.pvalue <- NA
for(i in 1:nrow(jcec)){
  if(is.na(jcec$PValue[i])){
    next
  }
  if(jcec$PValue[i] <= 0.05){
    jcec$predicted.pvalue[i] <- "+"
  }
  if(jcec$PValue[i] > 0.05){
    jcec$predicted.pvalue[i] <- "-"
  }
}

jcec$predicted.adjpvalue <- NA
for(i in 1:nrow(jcec)){
  if(is.na(jcec$FDR[i])){
    next
  }
  if(jcec$FDR[i] <= 0.05){
    jcec$predicted.adjpvalue[i] <- "+"
  }
  if(jcec$FDR[i] > 0.05){
    jcec$predicted.adjpvalue[i] <- "-"
  }
}

# rmats jc table

jc$predicted.pvalue <- NA
for(i in 1:nrow(jc)){
  if(is.na(jc$PValue[i])){
    next
  }
  if(jc$PValue[i] <= 0.05){
    jc$predicted.pvalue[i] <- "+"
  }
  if(jc$PValue[i] > 0.05){
    jc$predicted.pvalue[i] <- "-"
  }
}

jc$predicted.adjpvalue <- NA
for(i in 1:nrow(jc)){
  if(is.na(jc$FDR[i])){
    next
  }
  if(jc$FDR[i] <= 0.05){
    jc$predicted.adjpvalue[i] <- "+"
  }
  if(jc$FDR[i] > 0.05){
    jc$predicted.adjpvalue[i] <- "-"
  }
}

#actual pvalue and adjusted pvalue for each fold change
dxr1_exon2$actual.pvalue <- NA
for(i in 1:nrow(dxr1_exon2)){
  if(dxr1_exon2$FoldChange[i] == 0){
    dxr1_exon2$actual.pvalue[i] <- "-"
  }
  if(dxr1_exon2$FoldChange[i] > 0){
    dxr1_exon2$actual.pvalue[i] <- "+"
  }
  if(dxr1_exon2$FoldChange[i] < 0){
    dxr1_exon2$actual.pvalue[i] <- "+"
  }
}

jcec$actual.pvalue <- NA
for(i in 1:nrow(jcec)){
  if(jcec$FoldChange[i] == 0){
    jcec$actual.pvalue[i] <- "-"
  }
  if(jcec$FoldChange[i] > 0){
    jcec$actual.pvalue[i] <- "+"
  }
  if(jcec$FoldChange[i] < 0){
    jcec$actual.pvalue[i] <- "+"
  }
}

jc$actual.pvalue <- NA
for(i in 1:nrow(jc)){
  if(jc$FoldChange[i] == 0){
    jc$actual.pvalue[i] <- "-"
  }
  if(jc$FoldChange[i] > 0){
    jc$actual.pvalue[i] <- "+"
  }
  if(jc$FoldChange[i] < 0){
    jc$actual.pvalue[i] <- "+"
  }
}

dexseq_jcec_merge <- merge(dxr1_exon2, jcec, by=c("FoldChange", "ReadDepth", "SimRun"))
dexseq_jc_merge <- merge(dxr1_exon2, jc, by=c("FoldChange", "ReadDepth", "SimRun"))

table(dexseq_jcec_merge$predicted.pvalue.x, dexseq_jcec_merge$actual.pvalue.x, dnn=c("predicted", "actual"))
dexseq_spec_raw_pval <- 4479/(4479+27)
dexseq_sens_raw_pval <- 235/(235+8923)

table(dexseq_jcec_merge$predicted.adjpvalue.x, dexseq_jcec_merge$actual.pvalue.x, dnn=c("predicted", "actual"))
dexseq_spec_adjpval <- 4506/(4506+0)
dexseq_sens_adjpval <- 63/(63+9095)

table(dexseq_jcec_merge$predicted.pvalue.y, dexseq_jcec_merge$actual.pvalue.x, dnn=c("predicted", "actual"))
rmats_spec_raw_pval <- 4487/(4487+19)
rmats_sens_raw_pval <- 170/(170+8988)

table(dexseq_jcec_merge$predicted.adjpvalue.y, dexseq_jcec_merge$actual.pvalue.x, dnn=c("predicted", "actual"))
rmats_spec_adjpval <- 4500/(4500+6)
rmats_sens_adjpval <- 105/(105+9053)

value <- c(dexseq_spec_raw_pval, dexseq_sens_raw_pval, dexseq_spec_adjpval, dexseq_sens_adjpval,
           rmats_spec_raw_pval, rmats_sens_raw_pval, rmats_spec_adjpval, rmats_sens_adjpval)
software <- c(rep("DEXSeq", 4), rep("rMATS", 4))
evalmetric <- c(rep(c("Specificity", "Sensitivity"), 4))
pvalue <- c("raw pvalue", "raw pvalue", "adjusted pvalue", "adjusted pvalue", 
            "raw pvalue", "raw pvalue", "adjusted pvalue", "adjusted pvalue")

df1 <- data.frame(value, software, evalmetric, pvalue)

ggplot(df1, aes(fill=evalmetric, y=value, x=software)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~pvalue) + 
  theme_bw() +
  labs(fill = "Evaluation Metric", y="Sensitivity/Specificity")
ggsave("SensSpec_merged_rmats_dexseq_all.pdf", width = 6, height = 3, dpi = 100, units = "in", device='pdf')


#plot confusion matrix heatmap

cm <- confusionMatrix(factor(dexseq_jcec_merge$predicted.pvalue.x), factor(dexseq_jcec_merge$actual.pvalue.x), dnn = c("Prediction", "Reference"))

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=c("+", "-"))
plt$Reference <- factor(plt$Reference, levels=c("+", "-"))
plotTable <- plt %>%
  mutate(goodbad = ifelse(plt$Prediction == plt$Reference, "good", "bad"),
         software = "DEXSeq",
         pvalue = "raw pvalue") %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

cm <- confusionMatrix(factor(dexseq_jcec_merge$predicted.adjpvalue.x), factor(dexseq_jcec_merge$actual.pvalue.x), dnn = c("Prediction", "Reference"))

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=c("+", "-"))
plt$Reference <- factor(plt$Reference, levels=c("+", "-"))
plotTable2 <- plt %>%
  mutate(goodbad = ifelse(plt$Prediction == plt$Reference, "good", "bad"),
         software = "DEXSeq",
         pvalue = "adjusted pvalue") %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

cm <- confusionMatrix(factor(dexseq_jcec_merge$predicted.pvalue.y), factor(dexseq_jcec_merge$actual.pvalue.x), dnn = c("Prediction", "Reference"))

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=c("+", "-"))
plt$Reference <- factor(plt$Reference, levels=c("+", "-"))
plotTable3 <- plt %>%
  mutate(goodbad = ifelse(plt$Prediction == plt$Reference, "good", "bad"),
         software = "rMATS",
         pvalue = "raw pvalue") %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

cm <- confusionMatrix(factor(dexseq_jcec_merge$predicted.adjpvalue.y), factor(dexseq_jcec_merge$actual.pvalue.x), dnn = c("Prediction", "Reference"))

plt <- as.data.frame(cm$table)
plt$Prediction <- factor(plt$Prediction, levels=c("+", "-"))
plt$Reference <- factor(plt$Reference, levels=c("+", "-"))
plotTable4 <- plt %>%
  mutate(goodbad = ifelse(plt$Prediction == plt$Reference, "good", "bad"),
         software = "rMATS",
         pvalue = "adjusted pvalue") %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

comb_plotTable <- rbind(plotTable, plotTable2, plotTable3, plotTable4)

# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
ggplot(data = comb_plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "green", bad = "red")) +
  theme_bw() +
  labs(x = "Actual",y = "Predicted") +
  ylim(rev(levels(plt$Prediction))) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 15)) +
  facet_grid(rows=vars(software), cols=vars(pvalue))
ggsave("CM_merged_rmats_dexseq_all.pdf", width = 6, height = 4, dpi = 100, units = "in", device='pdf')

#############################################################################################################################

#stacked barplots for each folchange and read depth for raw pvalues

#FC 0 FP and TN counts
DEXSeqTN=c()
DEXSeqFP=c()
rMATS_jcecTN=c()
rMATS_jcecFP=c()
rMATS_jcTN=c()
rMATS_jcFP=c()

for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
  print(i)
  merged_filt <- filter(dexseq_jc_merge, FoldChange == 0 & ReadDepth == i)
  
  dexseq_2way_tbl <- table(merged_filt$predicted.pvalue.x)
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
  
  rmatsjc_2way_tbl <- table(merged_filt$predicted.pvalue.y)
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
  rmatsjcec_2way_tbl <- table(merged_filt$predicted.pvalue.y)
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

FP_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                    DEXSeqFP,
                    rMATS_jcecFP,
                    rMATS_jcFP)
FP_df <- FP_df %>% 
  gather(var, val, 2:ncol(FP_df ))
FP_df$ReadDepth <- factor(FP_df$ReadDepth, levels=c(1:10,50,100,500,1000))

ggplot(FP_df, aes(fill=var, y=val, x=ReadDepth)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="A3SS raw pvalue false positives", x="Read Depth", y="Count") +
  scale_fill_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  theme_bw()
ggsave("A3SS_FP_groupedbarchart_raw_pvalue.pdf", width = 6, height = 5, dpi = 100, units = "in", device='pdf')

TN_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                    DEXSeqTN,
                    rMATS_jcecTN,
                    rMATS_jcTN)
TN_df <- TN_df %>% 
  gather(var, val, 2:ncol(TN_df ))
TN_df$ReadDepth <- factor(TN_df$ReadDepth, levels=c(1:10,50,100,500,1000))

ggplot(TN_df, aes(fill=var, y=val, x=ReadDepth)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="A3SS raw pvalue true negatives", x="Read Depth", y="Count") +
  scale_fill_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  theme_bw()
ggsave("A3SS_TN_groupedbarchart_raw_pvalue.pdf", width = 6, height = 5, dpi = 100, units = "in", device='pdf')

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
    
    dexseq_2way_tbl <- table(merged_filt$predicted.pvalue.x)
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
    
    rmatsjc_2way_tbl <- table(merged_filt$predicted.pvalue.y)
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
    rmatsjcec_2way_tbl <- table(merged_filt$predicted.pvalue.y)
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

TP_df <- sens_tbl %>% 
  select(-c("DEXSeqFN", "rMATS_jcecFN", "rMATS_jcFN"))

TP_df <- TP_df %>% 
  gather(var, val, 3:ncol(TP_df ))
TP_df$ReadDepth <- factor(TP_df$ReadDepth, levels=c(1:10,50,100,500,1000))
TP_df$FoldChange <- factor(TP_df$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))

ggplot(TP_df, aes(fill=var, y=val, x=ReadDepth)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="A3SS raw pvalue true positives", x="Read Depth", y="Count") +
  scale_fill_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  theme_bw() + 
  facet_wrap(~FoldChange, nrow=2)
ggsave("A3SS_TP_groupedbarchart_raw_pvalue.pdf", width = 18, height = 7, dpi = 100, units = "in", device='pdf')

FN_df <- sens_tbl %>% 
  select(-c("DEXSeqTP", "rMATS_jcecTP", "rMATS_jcTP"))

FN_df <- FN_df %>% 
  gather(var, val, 3:ncol(FN_df ))
FN_df$ReadDepth <- factor(FN_df$ReadDepth, levels=c(1:10,50,100,500,1000))
FN_df$FoldChange <- factor(FN_df$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))

ggplot(FN_df, aes(fill=var, y=val, x=ReadDepth)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="A3SS raw pvalue false negatives", x="Read Depth", y="Count") +
  scale_fill_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  theme_bw() + 
  facet_wrap(~FoldChange, nrow=2)
ggsave("A3SS_FN_groupedbarchart_raw_pvalue.pdf", width = 18, height = 7, dpi = 100, units = "in", device='pdf')

#############################################################################################################################

#stacked barplots for each folchange and read depth for adjusted pvalues

#FC 0 FP and TN counts
DEXSeqTN=c()
DEXSeqFP=c()
rMATS_jcecTN=c()
rMATS_jcecFP=c()
rMATS_jcTN=c()
rMATS_jcFP=c()

for(i in unique(mixedsort(dexseq_jc_merge$ReadDepth))){
  print(i)
  merged_filt <- filter(dexseq_jc_merge, FoldChange == 0 & ReadDepth == i)
  
  dexseq_2way_tbl <- table(merged_filt$predicted.adjpvalue.x)
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
  
  rmatsjc_2way_tbl <- table(merged_filt$predicted.adjpvalue.y)
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
  rmatsjcec_2way_tbl <- table(merged_filt$predicted.adjpvalue.y)
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

FP_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                    DEXSeqFP,
                    rMATS_jcecFP,
                    rMATS_jcFP)
FP_df <- FP_df %>% 
  gather(var, val, 2:ncol(FP_df ))
FP_df$ReadDepth <- factor(FP_df$ReadDepth, levels=c(1:10,50,100,500,1000))

ggplot(FP_df, aes(fill=var, y=val, x=ReadDepth)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="A3SS adjusted pvalue false positives", x="Read Depth", y="Count") +
  scale_fill_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  theme_bw()
ggsave("A3SS_FP_groupedbarchart_adjpvalue.pdf", width = 6, height = 5, dpi = 100, units = "in", device='pdf')

TN_df <- data.frame(ReadDepth = c(1:10,50,100,500,1000),
                    DEXSeqTN,
                    rMATS_jcecTN,
                    rMATS_jcTN)
TN_df <- TN_df %>% 
  gather(var, val, 2:ncol(TN_df ))
TN_df$ReadDepth <- factor(TN_df$ReadDepth, levels=c(1:10,50,100,500,1000))

ggplot(TN_df, aes(fill=var, y=val, x=ReadDepth)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="A3SS adjusted pvalue true negatives", x="Read Depth", y="Count") +
  scale_fill_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  theme_bw()
ggsave("A3SS_TN_groupedbarchart_adjpvalue.pdf", width = 6, height = 5, dpi = 100, units = "in", device='pdf')

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
    
    dexseq_2way_tbl <- table(merged_filt$predicted.adjpvalue.x)
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
    
    rmatsjc_2way_tbl <- table(merged_filt$predicted.adjpvalue.y)
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
    rmatsjcec_2way_tbl <- table(merged_filt$predicted.adjpvalue.y)
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

TP_df <- sens_tbl %>% 
  select(-c("DEXSeqFN", "rMATS_jcecFN", "rMATS_jcFN"))

TP_df <- TP_df %>% 
  gather(var, val, 3:ncol(TP_df ))
TP_df$ReadDepth <- factor(TP_df$ReadDepth, levels=c(1:10,50,100,500,1000))
TP_df$FoldChange <- factor(TP_df$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))

ggplot(TP_df, aes(fill=var, y=val, x=ReadDepth)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="A3SS adjusted pvalue true positives", x="Read Depth", y="Count") +
  scale_fill_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  theme_bw() + 
  facet_wrap(~FoldChange, nrow=2)
ggsave("A3SS_TP_groupedbarchart_adjpvalue.pdf", width = 18, height = 7, dpi = 100, units = "in", device='pdf')

FN_df <- sens_tbl %>% 
  select(-c("DEXSeqTP", "rMATS_jcecTP", "rMATS_jcTP"))

FN_df <- FN_df %>% 
  gather(var, val, 3:ncol(FN_df ))
FN_df$ReadDepth <- factor(FN_df$ReadDepth, levels=c(1:10,50,100,500,1000))
FN_df$FoldChange <- factor(FN_df$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))

ggplot(FN_df, aes(fill=var, y=val, x=ReadDepth)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title="A3SS adjusted pvalue false negatives", x="Read Depth", y="Count") +
  scale_fill_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  theme_bw() + 
  facet_wrap(~FoldChange, nrow=2)
ggsave("A3SS_FN_groupedbarchart_adjpvalue.pdf", width = 18, height = 7, dpi = 100, units = "in", device='pdf')

