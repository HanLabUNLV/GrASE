#3/12/25 removing genes that show opposite results of rmats and dexseq
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
setwd("/home/jessicaholder1/master_thesis/Monaco_paper_data/B_vs_T/grase_results_NvB_vs_CD8T/results") #directory with grase results
##############################################################################
##############################################################################
##############################################################################
######---DEXSeq significant exons heatmap---######
dexseq_res_all <- read.table("Mapped.ExonsToEvents.txt",
                             header=TRUE, sep="\t")
#726256
dexseq_res_all_subset <- dexseq_res_all %>% dplyr::select(groupID, rMATS_ID, featureID, geneSymbol, pvalue, padj, log2fold_Nv_CD8_T_Nv_B,
                                                        countData.case_Nv_B, countData.case_Nv_CD8_T, padj,
                                                        FDR, IncLevel1, IncLevel2, IncLevelDifference)

dexseq_res_all_filt <- dexseq_res_all_subset %>% filter(padj <= 0.05 & FDR <= 0.05)
#121
dexseq_res_all_filt <- dexseq_res_all_filt %>% filter(abs(log2fold_Nv_CD8_T_Nv_B) > 1.0) #we want about 100 exon labels after filtering (changed from 2.6 which gave about 26 exons to 1.0 which gives about 86 exons in the final df)
#86
dexseq_res_all_filt <- dexseq_res_all_filt %>% separate(rMATS_ID, c("Event", "ID"))
dexseq_res_all_filt <- dexseq_res_all_filt %>%
  mutate(cell_type_dexseq = case_when(
    log2fold_Nv_CD8_T_Nv_B > 0 ~ "T",
    log2fold_Nv_CD8_T_Nv_B < 0 ~ "B"
  ))
dexseq_res_all_filt <- dexseq_res_all_filt %>%
  mutate(cell_type_rmats = case_when(
    IncLevelDifference > 0 ~ "B",
    IncLevelDifference < 0 ~ "T"
  ))
dexseq_res_all_filt_no_opp <- dexseq_res_all_filt[dexseq_res_all_filt$cell_type_rmats == dexseq_res_all_filt$cell_type_dexseq, ]
#74
dexseq_res_all_filt_final <- dexseq_res_all_filt_no_opp %>% group_by(across(featureID:geneSymbol)) %>% dplyr::slice(which.min(FDR))
###when rmats detects two events for one exonic part, this line says to keep the one with the lowest FDR
#67
dexseq_res_all_filt_final$DEXSeqID <- paste0(dexseq_res_all_filt_final$geneSymbol, "_", 
                                             dexseq_res_all_filt_final$featureID, "_",
                                             dexseq_res_all_filt_final$Event)
#67
dexseq_res_all_filt_final <- dexseq_res_all_filt_final %>% 
  separate(countData.case_Nv_B,
           c("countData.B1", "countData.B2","countData.B3", "countData.B4"))
dexseq_res_all_filt_final <- dexseq_res_all_filt_final %>% 
  separate(countData.case_Nv_CD8_T, 
           c("countData.T1", "countData.T2","countData.T3","countData.T4" ))
#A3SS
dexseq_A3SS <- dexseq_res_all_filt_final %>% 
  filter(Event == "A3SS") %>%
  ungroup() %>% 
  dplyr::select(DEXSeqID, 9:16) #columns where the counts are found
dexseq_A3SS <- as.data.frame(dexseq_A3SS)
rownames(dexseq_A3SS) <- dexseq_A3SS$DEXSeqID
dexseq_A3SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_A3SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_A3SS <- dexseq_A3SS[order(cutree(hclust_object, k = 2)), ]

#A5SS
dexseq_A5SS <- dexseq_res_all_filt_final %>% 
  filter(Event == "A5SS") %>%
  ungroup() %>% 
  dplyr::select(DEXSeqID, 9:16)
dexseq_A5SS <- as.data.frame(dexseq_A5SS)
rownames(dexseq_A5SS) <- dexseq_A5SS$DEXSeqID
dexseq_A5SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_A5SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_A5SS <- dexseq_A5SS[order(cutree(hclust_object, k = 2)), ]

#RI
dexseq_RI <- dexseq_res_all_filt_final %>% 
  filter(Event == "RI") %>%
  ungroup() %>% 
  dplyr::select(DEXSeqID, 9:16)
dexseq_RI <- as.data.frame(dexseq_RI)
rownames(dexseq_RI) <- dexseq_RI$DEXSeqID
dexseq_RI$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_RI, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_RI <- dexseq_RI[order(cutree(hclust_object, k = 2)), ]

#SE
dexseq_SE <- dexseq_res_all_filt_final %>% 
  filter(Event == "SE") %>%
  ungroup() %>% 
  dplyr::select(DEXSeqID, 9:16)
dexseq_SE <- as.data.frame(dexseq_SE)
rownames(dexseq_SE) <- dexseq_SE$DEXSeqID
dexseq_SE$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_SE, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_SE <- dexseq_SE[order(cutree(hclust_object, k = 2)), ]

dexseq_combined <- rbind(dexseq_A3SS, dexseq_A5SS, dexseq_RI, dexseq_SE)
#67
dexseq_combined[dexseq_combined==0] <- .00001 #set 0s to very small number so you can use the log2 function
dexseq_combined <- mutate_all(dexseq_combined, function(x) as.numeric(as.character(x)))
dexseq_combined <- log2(dexseq_combined)

dexseq_combined_2 <- dexseq_combined
dexseq_combined_2$ID <- rownames(dexseq_combined) #counts with ID as a column 
dexseq_combined_2_tmp <- merge(dexseq_combined_2, dexseq_res_all_filt_final, 
                               by.x = "ID", by.y = "DEXSeqID") #counts with information from the full DEXSeq results table
rownames(dexseq_combined_2_tmp) <- dexseq_combined_2_tmp$ID
AS_row_dexseq <- dexseq_combined_2_tmp %>% dplyr::select(Event) #need this to label splicing events on heatmap
col_annot <- data.frame(samples = c(rep("NvB", 4), rep("CD8T", 4)))
rownames(col_annot) <- colnames(dexseq_combined) #need this to label samples on heatmap


pdf("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/dexseq_sig_exonic_heatmap_no_opp.pdf", width=27, height=3.4)
pheatmap(t(dexseq_combined), 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot,
         annotation_col = AS_row_dexseq,
         angle_col = "270")
dev.off()
################################################################################
################################################################################
################################################################################

######---rMATS significant events heatmap---######
rMATS_res_all <- read.table("Mapped.EventsToExons.txt",
                             header=TRUE, sep="\t")
#173270
rMATS_res_all_filt <- rMATS_res_all %>% dplyr::select(GeneID, ID, DexseqFragment, geneSymbol, pvalue, padj, log2fold_Nv_CD8_T_Nv_B,
                                                      countData.case_Nv_B, countData.case_Nv_CD8_T, padj,
                                                      FDR, IncLevel1, IncLevel2, IncLevelDifference)
#173270
rMATS_res_all_filt <- rMATS_res_all_filt %>% filter(padj <= 0.05 & FDR <= 0.05)
#121
rMATS_res_all_filt <- rMATS_res_all_filt %>% filter(abs(IncLevelDifference) >= 0.1) #changed from .43 which gave 9 exons to .1 to give me 79 significant exons in the final dataframe
#79
rMATS_res_all_filt <- rMATS_res_all_filt %>% separate(ID, c("Event", "ID"))
rMATS_res_all_filt <- rMATS_res_all_filt %>%
  mutate(cell_type_dexseq = case_when(
    log2fold_Nv_CD8_T_Nv_B > 0 ~ "T",
    log2fold_Nv_CD8_T_Nv_B < 0 ~ "B"
  ))
rMATS_res_all_filt <- rMATS_res_all_filt %>%
  mutate(cell_type_rmats = case_when(
    IncLevelDifference > 0 ~ "B",
    IncLevelDifference < 0 ~ "T"
  ))
rMATS_res_all_filt_no_opp <- rMATS_res_all_filt[rMATS_res_all_filt$cell_type_rmats == rMATS_res_all_filt$cell_type_dexseq, ]
rMATS_res_all_filt_final <- rMATS_res_all_filt_no_opp %>% group_by(across(Event:ID)) %>% dplyr::slice(which.min(padj))
#57

rMATS_res_all_filt_final$DEXSeqID <- paste(rMATS_res_all_filt_final$geneSymbol, 
                                           rMATS_res_all_filt_final$DexseqFragment,
                                           rMATS_res_all_filt_final$ID,
                                           sep="_")
rMATS_res_all_filt_final <- rMATS_res_all_filt_final %>% 
  separate(IncLevel1, 
           c("IncLevel.B1","IncLevel.B2","IncLevel.B3","IncLevel.B4"), 
           sep=",")
rMATS_res_all_filt_final <- rMATS_res_all_filt_final %>% 
  separate(IncLevel2,
           c("IncLevel.T1", "IncLevel.T2","IncLevel.T3","IncLevel.T4"),
           sep=",")

#A3SS
rMATS_A3SS <- rMATS_res_all_filt_final %>% 
  filter(Event == "A3SS") %>% 
  ungroup() %>% 
  dplyr::select(DEXSeqID, 12:19) #columns where the inclusion levels are 
rMATS_A3SS <- as.data.frame(rMATS_A3SS)
rownames(rMATS_A3SS) <- rMATS_A3SS$DEXSeqID
rMATS_A3SS$DEXSeqID <- NULL
#Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(rMATS_A3SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
rMATS_A3SS <- rMATS_A3SS[order(cutree(hclust_object, k = 2)), ]

#A5SS
rMATS_A5SS <- rMATS_res_all_filt_final %>% 
  filter(Event == "A5SS") %>% 
  ungroup() %>% 
  dplyr::select(DEXSeqID, 12:19)
rMATS_A5SS <- as.data.frame(rMATS_A5SS)
rownames(rMATS_A5SS) <- rMATS_A5SS$DEXSeqID
rMATS_A5SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(rMATS_A5SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
rMATS_A5SS <- rMATS_A5SS[order(cutree(hclust_object, k = 2)), ]

#RI
rMATS_RI <- rMATS_res_all_filt_final %>% 
  filter(Event == "RI") %>% 
  ungroup() %>% 
  dplyr::select(DEXSeqID, 12:19)
rMATS_RI <- as.data.frame(rMATS_RI)
rownames(rMATS_RI) <- rMATS_RI$DEXSeqID
rMATS_RI$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(rMATS_RI, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
rMATS_RI <- rMATS_RI[order(cutree(hclust_object, k = 2)), ]

#SE
rMATS_SE <- rMATS_res_all_filt_final %>% 
  filter(Event == "SE") %>% 
  ungroup() %>% 
  dplyr::select(DEXSeqID, 12:19)
rMATS_SE <- as.data.frame(rMATS_SE)
rownames(rMATS_SE) <- rMATS_SE$DEXSeqID
rMATS_SE$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(rMATS_SE, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
rMATS_SE <- rMATS_SE[order(cutree(hclust_object, k = 2)), ]

rMATS_combined <- rbind(rMATS_A3SS, rMATS_A5SS, rMATS_RI, rMATS_SE)
#57
rMATS_combined <- rMATS_combined %>% dplyr::select(IncLevel.T4, IncLevel.T3, IncLevel.T2, IncLevel.T1,
                                                   IncLevel.B4, IncLevel.B3, IncLevel.B2, IncLevel.B1,
                                                   )
rMATS_combined <- mutate_all(rMATS_combined, function(x) as.numeric(as.character(x)))
rMATS_combined_2 <- rMATS_combined
rMATS_combined_2$ID <- rownames(rMATS_combined) #counts with ID as a column 
rMATS_combined_2_tmp <- merge(rMATS_combined_2, rMATS_res_all_filt_final, 
                               by.x = "ID", by.y = "DEXSeqID") #counts with information from the full DEXSeq results table
rownames(rMATS_combined_2_tmp) <- rMATS_combined_2_tmp$ID
AS_row_rMATS <- rMATS_combined_2_tmp %>% dplyr::select(Event) #need this to label splicing events on heatmap
col_annot_rmats <- data.frame(samples = c(rep("CD8T", 4), rep("NvB", 4))) 
rownames(col_annot_rmats) <- colnames(rMATS_combined) #need this to label samples on heatmap


pdf("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/rmats_sig_exonic_heatmap_no_opp.pdf", width=30, height=3.1)
pheatmap(t(rMATS_combined), 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rmats,
         annotation_col = AS_row_rMATS,
         annotation_colors = list(
           samples = c(`NvB`="violet", `CD8T`="salmon"),
           Event = c( A3SS= "darkkhaki", A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         angle_col = "270")
dev.off()

################################################################################
################################################################################
################################################################################
#obtaining heatmap with same IDs in both rMATS and DEXSeq

DEXSeq_sig_IDs <- data.frame(ID = rownames(dexseq_combined))
#67
DEXSeq_sig_IDs <- DEXSeq_sig_IDs %>% separate(ID, c("geneSymbol", "DEXSeqExonFragment", "Event"), sep = "_")
rMATS_sig_IDs <- data.frame(ID = rownames(rMATS_combined))
#57
rMATS_sig_IDs <- rMATS_sig_IDs %>% separate(ID, c("geneSymbol", "DEXSeqExonFragment", "rMATSID"), sep = "_")

sameIDs <- merge(DEXSeq_sig_IDs, rMATS_sig_IDs, by=c("geneSymbol", "DEXSeqExonFragment"))
#43
sameIDs$DEXSeqID <- paste0(sameIDs$geneSymbol, "_", sameIDs$DEXSeqExonFragment, "_", sameIDs$Event)
sameIDs$rMATSID <- paste0(sameIDs$geneSymbol, "_", sameIDs$DEXSeqExonFragment, "_", sameIDs$rMATSID)

# rMATS_sameIDs <- rMATS_combined_2_tmp[rownames(rMATS_combined_2_tmp) %in% sameIDs$rMATSID, ]
# # rMATS_sameIDs_filt <- subset(rMATS_sameIDs_full_df, (log2fold_Nv_CD8_T_Nv_B > 0 & IncLevelDifference > 0) | (log2fold_Nv_CD8_T_Nv_B < 0 & IncLevelDifference < 0)) #remove events that are opposite between rMATS and DEXSeq
# sameIDs <- merge(sameIDs, data.frame(ID = rMATS_sameIDs_filt$ID), by.x="rMATSID", by.y="ID") 
# Don't need this code? it was only for Raji and Jurkat data?
rMATS_sameIDs <- rMATS_combined[rownames(rMATS_combined) %in% sameIDs$rMATSID, ] 
AS_row_rMATS_sameIDs <- subset(AS_row_rMATS, (rownames(AS_row_rMATS) %in% sameIDs$rMATSID))
col_annot_rMATS <- col_annot_rmats
rownames(col_annot_rMATS) <- colnames(rMATS_sameIDs)
rMATS_sameIDs <- na.omit(rMATS_sameIDs)

pdf("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/heatmaps_used/rmats_sameIDs_heatmap_no_opp.pdf", width=12, height=3.4)
pheatmap(t(rMATS_sameIDs), ####################################################################### this heatmap is used in paper
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rMATS,
         annotation_col = AS_row_rMATS_sameIDs,
         annotation_colors = list(
           samples = c(`NvB`="violet", `CD8T`="salmon"),
           Event = c( A3SS= "darkkhaki", A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         color = brewer.pal(n = 9, name = "Blues"),
         angle_col = "270")

dev.off()
###need to check labels for this heatmap


rMATS_sameIDs$IDs <- rownames(rMATS_sameIDs)
rMATS_sameIDs <- rMATS_sameIDs %>% separate(IDs, c("geneSymbol", "DEXSeqExonFragment", "rMATSID"), sep = "_")
rMATS_sameIDs$DEXSeqID <- paste0(rMATS_sameIDs$geneSymbol, "_", rMATS_sameIDs$DEXSeqExonFragment)

dexseq_combined_2_tmp <- dexseq_combined_2_tmp %>% unite(col = "DEXSeqID", geneSymbol, featureID, sep= "_")
dexseq_combined_2_tmp <- dexseq_combined_2_tmp[dexseq_combined_2_tmp$DEXSeqID %in% rMATS_sameIDs$DEXSeqID, ]
dexseq_sameIDs <- dexseq_combined_2_tmp %>% dplyr::select(2:9) #selecting the count data
colnames(dexseq_sameIDs) <- colnames(dexseq_combined)

#match rows of rMATS sig to DEXSeq sig
dexseq_sameIDs$DEXSeqID <- rownames(dexseq_sameIDs)
dexseq_sameIDs <- dexseq_sameIDs %>% separate(DEXSeqID, c("geneSymbol", "DEXSeqExonFragment", "Event"), sep= "_")
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", geneSymbol, DEXSeqExonFragment, sep= "_")
dexseq_sameIDs <- left_join(rMATS_sameIDs, dexseq_sameIDs, by = "DEXSeqID")
dexseq_sameIDs <- dexseq_sameIDs[,12:21] # DexseqID, countData, and events columns
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", DEXSeqID, Event, sep= "_")
dexseq_sameIDs <- unique(dexseq_sameIDs)
rownames(dexseq_sameIDs) <- dexseq_sameIDs$DEXSeqID
dexseq_sameIDs$DEXSeqID <- NULL

AS_row_dexseq_sameIDs <- subset(AS_row_dexseq, (rownames(AS_row_dexseq) %in% rownames(dexseq_sameIDs)))
col_annot_dexseq <- col_annot
rownames(col_annot_dexseq) <- colnames(dexseq_sameIDs)
pdf("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/dexseq_sameIDs_heatmap_no_opp.pdf", width=8, height=3)
pheatmap(t(dexseq_sameIDs),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_sameIDs,
         annotation_colors = list(
           samples = c(`NvB`="violet", `CD8T`="salmon"),
           Event = c(A3SS= "darkkhaki", A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         color = brewer.pal(n = 9, name = "YlOrRd"),
         angle_col = "270")
dev.off()

################################################################################
################################################################################
################################################################################
#DEXSeq only significant exons

dexseq_only <- dexseq_combined[!(rownames(dexseq_combined) %in% sameIDs$DEXSeqID), ]
AS_row_dexseq_only <- subset(AS_row_dexseq, !(rownames(AS_row_dexseq) %in% rownames(AS_row_dexseq_sameIDs)))
pdf("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/dexseq_only_heatmap_no_opp.pdf", width=10, height=3)
pheatmap(t(dexseq_only),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_only,
         annotation_colors = list(
           samples = c(`NvB`="violet", `CD8T`="salmon"),
           Event = c(A3SS="darkkhaki", RI="cornflowerblue", SE="green3", A5SS="cyan2")),
         color = brewer.pal(n = 9, name = "YlOrRd"),
         angle_col = "270")
dev.off()

#rMATS only significant events

rMATS_only <- rMATS_combined[!(rownames(rMATS_combined) %in% sameIDs$rMATSID), ]
AS_row_rMATS_only <- subset(AS_row_rMATS, !(rownames(AS_row_rMATS) %in% rownames(AS_row_rMATS_sameIDs)))
pdf("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/heatmaps_used/rmats_only_heatmap_no_opp.pdf", width=15, height=3.4)
pheatmap(t(rMATS_only), ##############################################################################
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rMATS,
         annotation_col = AS_row_rMATS_only,
         annotation_colors = list(
           samples = c(`NvB`="violet",`CD8T`="salmon"),
           Event = c(A3SS="darkkhaki",  A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         color = brewer.pal(n = 9, name = "Blues"),
         angle_col = "270")
dev.off()

################################################################################
################################################################################



#using this/others counts instead of this counts-rerun dexseq to count(dxd) object found in dexseq_count.dxd.R script (copy over script from /home/jessicaholder1/master_thesis/Monaco_paper_data/Comapare_CM_LD.Neu/Heatmaps/dexseq_counts.dxd.R)
this_others<- read.table("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/dexseq_counts.txt")
this_others$DEXSeqID <- rownames(this_others)
this_others <- this_others %>% separate(DEXSeqID, c("groupID", "featureID"), sep=":")
this_others_merged <- merge(dexseq_res_all_filt_final, this_others, by = c("groupID", "featureID"))
write.csv(this_others_merged, file = "/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/this_others_merged.csv")
this_others <- this_others %>% mutate(B1 = V1/(V1+V9),
                                      B2 = V2/(V2+V10),
                                      B3 = V3/(V3+V11),
                                      B4 = V4/(V4+V12),
                                      T1 = V5/(V5+V13),
                                      T2 = V6/(V6+V14),
                                      T3 = V7/(V6+V15),
                                      T4 = V8/(V6+V16))
this_others <- this_others[,17:26] #the columns you mutated in the code above and the groupID
this_others$DEXSeqID <- rownames(this_others)
this_others <- this_others %>% separate(DEXSeqID, c("groupID", "featureID"), sep=":") #really just moves the groupID and featureID columns to the last columns 
this_others_merged <- merge(dexseq_res_all_filt_final, this_others, by = c("groupID", "featureID"))

#A3SS
dexseq_A3SS <- this_others_merged %>% filter(Event == "A3SS") %>% dplyr::select(DEXSeqID, 24:31) #columns were counts
  rownames(dexseq_A3SS) <- dexseq_A3SS$DEXSeqID
dexseq_A3SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_A3SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_A3SS <- dexseq_A3SS[order(cutree(hclust_object, k = 2)), ]

#A5SS
dexseq_A5SS <- this_others_merged %>% filter(Event == "A5SS") %>% dplyr::select(DEXSeqID, 24:31)
rownames(dexseq_A5SS) <- dexseq_A5SS$DEXSeqID
dexseq_A5SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_A5SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_A5SS <- dexseq_A5SS[order(cutree(hclust_object, k = 2)), ]

#RI
dexseq_RI <- this_others_merged %>% filter(Event == "RI") %>% dplyr::select(DEXSeqID, 24:31)
rownames(dexseq_RI) <- dexseq_RI$DEXSeqID
dexseq_RI$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_RI, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_RI <- dexseq_RI[order(cutree(hclust_object, k = 2)), ]

#SE
dexseq_SE <- this_others_merged %>% filter(Event == "SE") %>% dplyr::select(DEXSeqID, 24:31)
rownames(dexseq_SE) <- dexseq_SE$DEXSeqID
dexseq_SE$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_SE, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_SE <- dexseq_SE[order(cutree(hclust_object, k = 2)), ]

dexseq_combined <- rbind(dexseq_A3SS, dexseq_A5SS, dexseq_RI, dexseq_SE)
#dexseq_combined <- dexseq_combined %>% dplyr::select(countData.CD8naive1, countData.CD8naive2, countData.CD8naive3, countData.CD8naive4, countData.CD8naive5,
#                                                     countData.Bnaive1, countData.Bnaive2, countData.Bnaive3, countData.Bnaive4, countData.Bnaive5) 
#dexseq_combined[dexseq_combined==0] <- .00001 #set 0s to very small number so you can use the log2 function
dexseq_combined <- mutate_all(dexseq_combined, function(x) as.numeric(as.character(x)))
dexseq_combined <- log2(dexseq_combined)

dexseq_combined_2 <- dexseq_combined
dexseq_combined_2$ID <- rownames(dexseq_combined) #counts with ID as a column 
dexseq_combined_2_tmp <- merge(dexseq_combined_2, dexseq_res_all_filt_final, 
                               by.x = "ID", by.y = "DEXSeqID") #counts with information from the full DEXSeq results table
rownames(dexseq_combined_2_tmp) <- dexseq_combined_2_tmp$ID
AS_row_dexseq <- dexseq_combined_2_tmp %>% dplyr::select(Event) #need this to label splicing events on heatmap
col_annot <- data.frame(samples = c(rep("NvB", 4), rep("CD8T", 4))) 
rownames(col_annot) <- colnames(dexseq_combined) #need this to label samples on heatmap


pdf("dexseq_sig_exonic_heatmap_this_others_no_opp.pdf", width=18, height=3.4)
pheatmap(t(dexseq_combined),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot,
         annotation_col = AS_row_dexseq,
         color = brewer.pal(n = 9, name = "YlOrRd"),
         angle_col = "270")
dev.off()

dexseq_sameIDs <- dexseq_combined[(rownames(dexseq_combined) %in% sameIDs$DEXSeqID), ]
dexseq_sameIDs$DEXSeqID <- rownames(dexseq_sameIDs)
dexseq_sameIDs <- dexseq_sameIDs %>% separate(DEXSeqID, c("geneSymbol", "DEXSeqExonFragment", "Event"), sep="_")
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", geneSymbol, DEXSeqExonFragment, sep="_")
dexseq_sameIDs <- left_join(rMATS_sameIDs, dexseq_sameIDs, by = "DEXSeqID")
#now has 43
dexseq_sameIDs <- dexseq_sameIDs[,12:21] #countData, DexseqID, and events columns
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", DEXSeqID, Event, sep="_")
dexseq_sameIDs <- unique(dexseq_sameIDs)
#then changes to 40
rownames(dexseq_sameIDs) <- dexseq_sameIDs$DEXSeqID
dexseq_sameIDs$DEXSeqID <- NULL

AS_row_dexseq_sameIDs <- subset(AS_row_dexseq, (rownames(AS_row_dexseq) %in% rownames(dexseq_sameIDs)))
col_annot_dexseq <- col_annot
rownames(col_annot_dexseq) <- colnames(dexseq_sameIDs)
library(dplyr)
dexseq_sameIDs <- na.omit(dexseq_sameIDs)
pdf("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/heatmaps_used/dexseq_sameIDs_heatmap_this_others_no_opp.pdf", width=12, height=3.4)
pheatmap(t(dexseq_sameIDs), ############################################################################### this heatmap is used in the paper
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_sameIDs,
         annotation_colors = list(
           samples = c(`NvB`="violet", `CD8T`="salmon"),
           Event = c(A3SS="darkkhaki",  A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(25),
         angle_col = "270")
#### THE LABELS ARE CORRECT, I think it was changing the col_annot
write.csv(dexseq_sameIDs, file = "/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/heatmaps_used/dexseq_heatmap_sameIDs_no_opp.csv")
dev.off()

dexseq_only <- dexseq_combined[!(rownames(dexseq_combined) %in% sameIDs$DEXSeqID), ]
AS_row_dexseq_only <- subset(AS_row_dexseq, !(rownames(AS_row_dexseq) %in% rownames(AS_row_dexseq_sameIDs)))
#dexseq_only <- dexseq_only %>% relocate(countData.T1, countData.T2, countData.T3, countData.T4) #move those columns to the front so orientation of samples is correct for heatmaps
pdf("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/heatmaps_used/dexseq_only_heatmap_this_others_no_opp.pdf", width=15, height=3)
pheatmap(t(dexseq_only), ####################################################################################
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_only,
         annotation_colors = list(
           samples = c(`NvB`="violet", `CD8T`="salmon"),
           Event = c(A3SS="darkkhaki",  A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(25),
         angle_col = "270")
#### 
write.csv(dexseq_only, file = "/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps/heatmaps_used/dexseq_heatmap_only_no_opp.csv")
dev.off()
##
##when it says this_others (this/others) I mean this_total(this/others+this) 
setwd("/home/jessicaholder1/master_thesis/Monaco_paper_data/Compare_B_T/Heatmaps")
write.csv(rMATS_sameIDs, "Monaco_rmats_sameIDs_no_opp.csv")
write.csv(rMATS_only, "Monaco_rmats_only_no_opp.csv")
write.csv(dexseq_sameIDs, "Monaco_dexseq_sameIDs_no_opp.csv")
write.csv(dexseq_only, "Monaco_dexseq_only_no_opp.csv")

