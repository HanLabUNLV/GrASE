library(tidyverse)
library(pheatmap)
library(RColorBrewer)
setwd("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/grase_results_v6/results") #directory with grase results
##############################################################################
##############################################################################
##############################################################################
######---DEXSeq significant exons heatmap---######
dexseq_res_all <- read.table("Mapped.ExonsToEvents.txt",
                             header=TRUE, sep="\t")
dexseq_res_all_filt <- dexseq_res_all %>% dplyr::select(groupID, rMATS_ID, featureID, geneSymbol, pvalue, padj, log2fold_Raji_Jurkat,
                                                        countData.case_raji, countData.case_jurkat, padj,
                                                        FDR, IncLevel1, IncLevel2, IncLevelDifference)

dexseq_res_all_filt <- dexseq_res_all_filt %>% filter(padj <= 0.05 & FDR <= 0.05) ##filtered dexseq by
dexseq_res_all_filt <- dexseq_res_all_filt %>% filter(abs(log2fold_Raji_Jurkat) > 2.4) #we want about 100 exon labels after filtering (changed from 1.2 which gave about 350 exons to 2.4 which gives about 101 exons in the final df)
dexseq_res_all_filt <- dexseq_res_all_filt %>% separate(rMATS_ID, c("Event", "ID"))

dexseq_res_all_filt_final <- dexseq_res_all_filt %>% group_by(across(featureID:geneSymbol)) %>% slice(which.min(FDR))

dexseq_res_all_filt_final$DEXSeqID <- paste0(dexseq_res_all_filt_final$geneSymbol, "_", 
                                             dexseq_res_all_filt_final$featureID, "_",
                                             dexseq_res_all_filt_final$Event)  
dexseq_res_all_filt_final <- dexseq_res_all_filt_final %>% 
  separate(countData.case_raji,
           c("countData.raji1", "countData.raji2","countData.raji3"))
dexseq_res_all_filt_final <- dexseq_res_all_filt_final %>% 
  separate(countData.case_jurkat, 
           c("countData.jurkat1", "countData.jurkat2","countData.jurkat3"))

#A3SS
dexseq_A3SS <- dexseq_res_all_filt_final %>% 
  filter(Event == "A3SS") %>%
  ungroup() %>% 
  dplyr::select(DEXSeqID, 9:14) #columns where the counts are found
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
  dplyr::select(DEXSeqID, 9:14)
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
  dplyr::select(DEXSeqID, 9:14)
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
  dplyr::select(DEXSeqID, 9:14)
dexseq_SE <- as.data.frame(dexseq_SE)
rownames(dexseq_SE) <- dexseq_SE$DEXSeqID
dexseq_SE$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_SE, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_SE <- dexseq_SE[order(cutree(hclust_object, k = 2)), ]

dexseq_combined <- rbind(dexseq_A3SS, dexseq_A5SS, dexseq_RI, dexseq_SE)
dexseq_combined[dexseq_combined==0] <- .00001 #set 0s to very small number so you can use the log2 function
dexseq_combined <- mutate_all(dexseq_combined, function(x) as.numeric(as.character(x)))
dexseq_combined <- log2(dexseq_combined)

dexseq_combined_2 <- dexseq_combined
dexseq_combined_2$ID <- rownames(dexseq_combined) #counts with ID as a column 
dexseq_combined_2_tmp <- merge(dexseq_combined_2, dexseq_res_all_filt_final, 
                               by.x = "ID", by.y = "DEXSeqID") #counts with information from the full DEXSeq results table
rownames(dexseq_combined_2_tmp) <- dexseq_combined_2_tmp$ID
AS_row_dexseq <- dexseq_combined_2_tmp %>% dplyr::select(Event) #need this to label splicing events on heatmap
#col_annot <- data.frame(samples = c(rep("CD8 naive", 5), rep("B naive", 5))) 
col_annot <- data.frame(samples = c(rep("Raji", 3), rep("Jurkat", 3)))
rownames(col_annot) <- colnames(dexseq_combined) #need this to label samples on heatmap


pdf("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/heatmaps/dexseq_sig_exonic_heatmap.pdf", width=27, height=3.4)
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
rMATS_res_all_filt <- rMATS_res_all %>% dplyr::select(GeneID, ID, DexseqFragment, geneSymbol, pvalue, padj, log2fold_Raji_Jurkat,
                                                      countData.case_raji, countData.case_jurkat, padj,
                                                      FDR, IncLevel1, IncLevel2, IncLevelDifference)
rMATS_res_all_filt <- rMATS_res_all_filt %>% filter(padj <= 0.05 & FDR <= 0.05) #filtered rmats by
rMATS_res_all_filt <- rMATS_res_all_filt %>% filter(abs(IncLevelDifference) >= 0.38) #changed to .38 to give me around 100 significant exons in the final dataframe
rMATS_res_all_filt <- rMATS_res_all_filt %>% separate(ID, c("Event", "ID"))
rMATS_res_all_filt_final <- rMATS_res_all_filt %>% group_by(across(Event:ID)) %>% slice(which.min(padj))


rMATS_res_all_filt_final$DEXSeqID <- paste(rMATS_res_all_filt_final$geneSymbol, 
                                           rMATS_res_all_filt_final$DexseqFragment,
                                           rMATS_res_all_filt_final$ID,
                                           sep="_")
rMATS_res_all_filt_final <- rMATS_res_all_filt_final %>% 
  separate(IncLevel1, 
           c("IncLevel.raji1", "IncLevel.raji2","IncLevel.raji3"), 
           sep=",")
rMATS_res_all_filt_final <- rMATS_res_all_filt_final %>% 
  separate(IncLevel2,
           c("IncLevel.jurkat1", "IncLevel.jurkat2","IncLevel.jurkat3"),
           sep=",")


#A3SS
rMATS_A3SS <- rMATS_res_all_filt_final %>% 
  filter(Event == "A3SS") %>% 
  ungroup() %>% 
  dplyr::select(DEXSeqID, 12:17) #columns where the inclusion levels are 
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
  dplyr::select(DEXSeqID, 12:17)
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
  dplyr::select(DEXSeqID, 12:17)
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
  dplyr::select(DEXSeqID, 12:17)
rMATS_SE <- as.data.frame(rMATS_SE)
rownames(rMATS_SE) <- rMATS_SE$DEXSeqID
rMATS_SE$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(rMATS_SE, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
rMATS_SE <- rMATS_SE[order(cutree(hclust_object, k = 2)), ]

rMATS_combined <- rbind(rMATS_A3SS, rMATS_A5SS, rMATS_RI, rMATS_SE)
rMATS_combined <- rMATS_combined %>% dplyr::select(IncLevel.jurkat3, IncLevel.jurkat2, IncLevel.jurkat1,
                                                  IncLevel.raji3, IncLevel.raji2, IncLevel.raji1,
                                                   )
rMATS_combined <- mutate_all(rMATS_combined, function(x) as.numeric(as.character(x)))
rMATS_combined_2 <- rMATS_combined
rMATS_combined_2$ID <- rownames(rMATS_combined) #counts with ID as a column 
rMATS_combined_2_tmp <- merge(rMATS_combined_2, rMATS_res_all_filt_final, 
                               by.x = "ID", by.y = "DEXSeqID") #counts with information from the full DEXSeq results table
rownames(rMATS_combined_2_tmp) <- rMATS_combined_2_tmp$ID
AS_row_rMATS <- rMATS_combined_2_tmp %>% dplyr::select(Event) #need this to label splicing events on heatmap
col_annot_rmats <- data.frame(samples = c(rep("Jurkat", 3), rep("Raji", 3))) 
rownames(col_annot_rmats) <- colnames(rMATS_combined) #need this to label samples on heatmap


pdf("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/heatmaps/rmats_sig_exonic_heatmap.pdf", width=30, height=3.1)
pheatmap(t(rMATS_combined), 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rmats,
         annotation_col = AS_row_rMATS,
         annotation_colors = list(
           samples = c(`Raji`="violet", `Jurkat`="salmon"),
           Event = c( A3SS= "darkkhaki", A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         angle_col = "270")
dev.off()
#####this one is the correct orientation: Raji is on the bottom and it is blue because its inc level is lower than Jurkat which is red because its inc level is 1
################################################################################
################################################################################
################################################################################
#obtaining heatmap with same IDs in both rMATS and DEXSeq

DEXSeq_sig_IDs <- data.frame(ID = rownames(dexseq_combined))
DEXSeq_sig_IDs <- DEXSeq_sig_IDs %>% separate(ID, c("geneSymbol", "DEXSeqExonFragment", "Event"), sep = "_")
rMATS_sig_IDs <- data.frame(ID = rownames(rMATS_combined))
rMATS_sig_IDs <- rMATS_sig_IDs %>% separate(ID, c("geneSymbol", "DEXSeqExonFragment", "rMATSID"), sep = "_")

sameIDs <- merge(DEXSeq_sig_IDs, rMATS_sig_IDs, by=c("geneSymbol", "DEXSeqExonFragment"))
sameIDs$DEXSeqID <- paste0(sameIDs$geneSymbol, "_", sameIDs$DEXSeqExonFragment, "_", sameIDs$Event)
sameIDs$rMATSID <- paste0(sameIDs$geneSymbol, "_", sameIDs$DEXSeqExonFragment, "_", sameIDs$rMATSID)

rMATS_sameIDs_full_df <- rMATS_combined_2_tmp[rownames(rMATS_combined_2_tmp) %in% sameIDs$rMATSID, ]
rMATS_sameIDs_filt <- subset(rMATS_sameIDs_full_df, (log2fold_Raji_Jurkat > 0 & IncLevelDifference > 0) | (log2fold_Raji_Jurkat < 0 & IncLevelDifference < 0)) #remove events that are opposite between rMATS and DEXSeq
sameIDs <- merge(sameIDs, data.frame(ID = rMATS_sameIDs_filt$ID), by.x="rMATSID", by.y="ID") 

rMATS_sameIDs <- rMATS_combined[rownames(rMATS_combined) %in% sameIDs$rMATSID, ]
AS_row_rMATS_sameIDs <- subset(AS_row_rMATS, (rownames(AS_row_rMATS) %in% sameIDs$rMATSID))
col_annot_rMATS <- col_annot_rmats
rownames(col_annot_rMATS) <- colnames(rMATS_sameIDs)

pdf("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/heatmaps/heatmaps_used/rmats_sameIDs_heatmap.pdf", width=12, height=3.4)
pheatmap(t(rMATS_sameIDs), #######################################################################
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rMATS,
         annotation_col = AS_row_rMATS_sameIDs,
         annotation_colors = list(
           samples = c(`Raji`="violet", `Jurkat`="salmon"),
           Event = c( A3SS= "darkkhaki", A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         color = brewer.pal(n = 9, name = "Blues"),
         angle_col = "270")
dev.off()
##this is the correct orientation: raji is blue because it has a low inc level difference (found in the rmats_res_all_filt) and should be on the bottom so when inverted will be on the top


rMATS_sameIDs$IDs <- rownames(rMATS_sameIDs)
rMATS_sameIDs <- rMATS_sameIDs %>% separate(IDs, c("geneSymbol", "DEXSeqExonFragment", "rMATSID"), sep = "_")
rMATS_sameIDs$DEXSeqID <- paste0(rMATS_sameIDs$geneSymbol, "_", rMATS_sameIDs$DEXSeqExonFragment)

dexseq_combined_2_tmp <- dexseq_combined_2_tmp %>% unite(col = "DEXSeqID", geneSymbol, featureID, sep="_")
dexseq_combined_2_tmp <- dexseq_combined_2_tmp[dexseq_combined_2_tmp$DEXSeqID %in% rMATS_sameIDs$DEXSeqID, ]
dexseq_sameIDs <- dexseq_combined_2_tmp %>% dplyr::select(2:7) 
colnames(dexseq_sameIDs) <- colnames(dexseq_combined)

#match rows of rMATS sig to DEXSeq sig
dexseq_sameIDs$DEXSeqID <- rownames(dexseq_sameIDs)
dexseq_sameIDs <- dexseq_sameIDs %>% separate(DEXSeqID, c("geneSymbol", "DEXSeqExonFragment", "Event"))
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", geneSymbol, DEXSeqExonFragment, sep="_")
dexseq_sameIDs <- left_join(rMATS_sameIDs, dexseq_sameIDs, by = "DEXSeqID")
dexseq_sameIDs <- dexseq_sameIDs[,10:17] #count, DexseqID, and events columns
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", DEXSeqID, Event, sep="_")
dexseq_sameIDs <- unique(dexseq_sameIDs)
rownames(dexseq_sameIDs) <- dexseq_sameIDs$DEXSeqID
dexseq_sameIDs$DEXSeqID <- NULL

AS_row_dexseq_sameIDs <- subset(AS_row_dexseq, (rownames(AS_row_dexseq) %in% rownames(dexseq_sameIDs)))
col_annot_dexseq <- col_annot
rownames(col_annot_dexseq) <- colnames(dexseq_sameIDs)
pdf("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/heatmaps/dexseq_sameIDs_heatmap.pdf", width=8, height=3)
pheatmap(t(dexseq_sameIDs),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_sameIDs,
         annotation_colors = list(
           samples = c(`Raji`="violet", `Jurkat`="salmon"),
           Event = c(A5SS="cyan2", SE="green3")),
         color = brewer.pal(n = 9, name = "YlOrRd"),
         angle_col = "270")
dev.off()

################################################################################
################################################################################
################################################################################
#DEXSeq only significant exons

dexseq_only <- dexseq_combined[!(rownames(dexseq_combined) %in% sameIDs$DEXSeqID), ]
AS_row_dexseq_only <- subset(AS_row_dexseq, !(rownames(AS_row_dexseq) %in% rownames(AS_row_dexseq_sameIDs)))
pdf("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/heatmaps/dexseq_only_heatmap.pdf", width=10, height=3)
pheatmap(t(dexseq_only),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_only,
         annotation_colors = list(
           samples = c(`Raji`="violet", `Jurkat`="salmon"),
           Event = c(A3SS="darkkhaki", RI="cornflowerblue", SE="green3", A5SS="cyan2")),
         color = brewer.pal(n = 9, name = "YlOrRd"),
         angle_col = "270")
dev.off()

#rMATS only significant events

rMATS_only <- rMATS_combined[!(rownames(rMATS_combined) %in% sameIDs$rMATSID), ]
AS_row_rMATS_only <- subset(AS_row_rMATS, !(rownames(AS_row_rMATS) %in% rownames(AS_row_rMATS_sameIDs)))
pdf("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/heatmaps/heatmaps_used/rmats_only_heatmap.pdf", width=15, height=3.4)
pheatmap(t(rMATS_only), ##############################################################################
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rMATS,
         annotation_col = AS_row_rMATS_only,
         annotation_colors = list(
           samples = c(`Raji`="violet",`Jurkat`="salmon"),
           Event = c(A3SS="darkkhaki",  A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         color = brewer.pal(n = 9, name = "Blues"),
         angle_col = "270")
dev.off()
##this is correct Raji is on the bottom and is blue because the inc level is very low
################################################################################
################################################################################
#using this/others counts instead of this counts-rerun dexseq to count(dxd) object found in dexseq_count.dxd.R script
this_others_RJ <- read.table("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/dexseq_counts.txt")
this_others_RJ$DEXSeqID <- rownames(this_others_RJ)
this_others_RJ <- this_others_RJ %>% separate(DEXSeqID, c("groupID", "featureID"), sep=":")
this_others_merged <- merge(dexseq_res_all_filt_final, this_others, by = c("groupID", "featureID"))
write.csv(this_others_merged, file = "this_others_merged.csv")
this_others <- this_others %>% mutate(Jurkat1 = V1/(V1+V7), #V1 is that sample counts and V1+V7 are the other counts
                                      Jurkat2 = V2/(V2+V8),
                                      Jurkat3 = V3/(V3+V9),
                                      Raji1 = V4/(V4+V10),
                                      Raji2 = V5/(V5+V11),
                                      Raji3 = V6/(V6+V12))
this_others <- this_others[,13:20] #the columns you mutated in the code above and the groupID
this_others$DEXSeqID <- rownames(this_others)
this_others <- this_others %>% separate(DEXSeqID, c("groupID", "featureID"), sep=":")
this_others_merged <- merge(dexseq_res_all_filt_final, this_others, by = c("groupID", "featureID"))

#A3SS
dexseq_A3SS <- this_others_merged %>% filter(Event == "A3SS") %>% dplyr::select(DEXSeqID, 19:25)
  rownames(dexseq_A3SS) <- dexseq_A3SS$DEXSeqID
dexseq_A3SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_A3SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_A3SS <- dexseq_A3SS[order(cutree(hclust_object, k = 2)), ]

#A5SS
dexseq_A5SS <- this_others_merged %>% filter(Event == "A5SS") %>% dplyr::select(DEXSeqID, 19:25)
rownames(dexseq_A5SS) <- dexseq_A5SS$DEXSeqID
dexseq_A5SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_A5SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_A5SS <- dexseq_A5SS[order(cutree(hclust_object, k = 2)), ]

#RI
dexseq_RI <- this_others_merged %>% filter(Event == "RI") %>% dplyr::select(DEXSeqID, 19:25)
rownames(dexseq_RI) <- dexseq_RI$DEXSeqID
dexseq_RI$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_RI, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_RI <- dexseq_RI[order(cutree(hclust_object, k = 2)), ]

#SE
dexseq_SE <- this_others_merged %>% filter(Event == "SE") %>% dplyr::select(DEXSeqID, 19:25)
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
col_annot <- data.frame(samples = c(rep("Jurkat", 3), rep("Raji", 3))) 
rownames(col_annot) <- colnames(dexseq_combined) #need this to label samples on heatmap


pdf("dexseq_sig_exonic_heatmap_this_others.pdf", width=18, height=3.4)
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
dexseq_sameIDs <- dexseq_sameIDs %>% separate(DEXSeqID, c("geneSymbol", "DEXSeqExonFragment", "Event"))
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", geneSymbol, DEXSeqExonFragment, sep="_")
dexseq_sameIDs <- left_join(rMATS_sameIDs, dexseq_sameIDs, by = "DEXSeqID")
dexseq_sameIDs <- dexseq_sameIDs[,10:17]
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", DEXSeqID, Event, sep="_")
dexseq_sameIDs <- unique(dexseq_sameIDs)
rownames(dexseq_sameIDs) <- dexseq_sameIDs$DEXSeqID
dexseq_sameIDs$DEXSeqID <- NULL

AS_row_dexseq_sameIDs <- subset(AS_row_dexseq, (rownames(AS_row_dexseq) %in% rownames(dexseq_sameIDs)))
col_annot_dexseq <- col_annot
rownames(col_annot_dexseq) <- colnames(dexseq_sameIDs)
library(dplyr)
dexseq_sameIDs <- dexseq_sameIDs %>% relocate(Raji1,Raji2, Raji3) #move those columns to the front so orientation of samples is correct for heatmaps
pdf("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/heatmaps/heatmaps_used/dexseq_sameIDs_heatmap_this_others.pdf", width=12, height=3.4)
pheatmap(t(dexseq_sameIDs), ###############################################################################
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_sameIDs,
         annotation_colors = list(
           samples = c(`Raji`="violet", `Jurkat`="salmon"),
           Event = c(A5SS="cyan2", SE="green3")),
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(25),
         angle_col = "270")
#write.csv(dexseq_sameIDs, file = "dexseq_heatmap_sameIDs.csv")
dev.off()
##this is correct
dexseq_only <- dexseq_combined[!(rownames(dexseq_combined) %in% sameIDs$DEXSeqID), ]
AS_row_dexseq_only <- subset(AS_row_dexseq, !(rownames(AS_row_dexseq) %in% rownames(AS_row_dexseq_sameIDs)))
dexseq_only <- dexseq_only %>% relocate(Raji1,Raji2, Raji3) #move those columns to the front so orientation of samples is correct for heatmaps
pdf("/home/jessicaholder1/master_thesis/Raji_Jurkat_Cells/heatmaps/heatmaps_used/dexseq_only_heatmap_this_others.pdf", width=15, height=3)
pheatmap(t(dexseq_only), ####################################################################################
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_only,
         annotation_colors = list(
           samples = c(`Raji`="violet", `Jurkat`="salmon"),
           Event = c(A3SS="darkkhaki",  A5SS="cyan2", RI="cornflowerblue", SE="green3")),
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(25),
         angle_col = "270")
#write.csv(dexseq_only, file = "dexseq_heatmap_only.csv")
dev.off()
##this is correct
##when it says this_others (this/others) I mean this_total(this/others+this) 
