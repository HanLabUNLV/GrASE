library(tidyverse)
library(pheatmap)
setwd("/mnt/storage/jaquino/GrASE_results_bonnal/grase_results_v2")

#DEXSeq significant exons heatmap
dexseq_res_all <- read.table("Mapped.ExonsToEvents.txt",
                             header=TRUE, sep="\t")
dexseq_res_all_filt <- dexseq_res_all %>% dplyr::select(groupID, rMATS_ID, featureID, geneSymbol, pvalue, padj, log2fold_CD8naive_Bnaive,
                                                 countData.case_Bnaive, countData.case_CD8naive, padj,
                                                 FDR, IncLevel1, IncLevel2, IncLevelDifference)

dexseq_res_all_filt <- dexseq_res_all_filt %>% filter(padj <= 0.05 & FDR <= 0.05)
dexseq_res_all_filt <- dexseq_res_all_filt %>% filter(abs(log2fold_CD8naive_Bnaive) > 1.2)
dexseq_res_all_filt <- dexseq_res_all_filt %>% separate(rMATS_ID, c("Event", "ID"))

dexseq_res_all_filt_final <- dexseq_res_all_filt %>% group_by(across(featureID:geneSymbol)) %>% slice(which.min(FDR))

dexseq_res_all_filt_final$DEXSeqID <- paste0(dexseq_res_all_filt_final$geneSymbol, "_", 
                                             dexseq_res_all_filt_final$featureID, "_",
                                             dexseq_res_all_filt_final$Event)  
dexseq_res_all_filt_final <- dexseq_res_all_filt_final %>% 
  separate(countData.case_Bnaive,
           c("countData.Bnaive1", "countData.Bnaive2","countData.Bnaive3", "countData.Bnaive4", "countData.Bnaive5"))
dexseq_res_all_filt_final <- dexseq_res_all_filt_final %>% 
  separate(countData.case_CD8naive, 
           c("countData.CD8naive1", "countData.CD8naive2","countData.CD8naive3", "countData.CD8naive4", "countData.CD8naive5"))

#A3SS
dexseq_A3SS <- dexseq_res_all_filt_final %>% 
  filter(Event == "A3SS") %>%
  ungroup() %>% 
  dplyr::select(DEXSeqID, 9:18)
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
  dplyr::select(DEXSeqID, 9:18)
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
  dplyr::select(DEXSeqID, 9:18)
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
  dplyr::select(DEXSeqID, 9:18)
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
col_annot <- data.frame(samples = c(rep("B naive", 5), rep("CD8 naive", 5))) 
rownames(col_annot) <- colnames(dexseq_combined) #need this to label samples on heatmap


pdf("dexseq_sig_exonic_heatmap.pdf", width=18, height=3.4)
pheatmap(t(dexseq_combined),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot,
         annotation_col = AS_row_dexseq,
         angle_col = "270")
dev.off()
#######################################################################################################################################################
#rMATS significant events heatmap
rMATS_res_all <- read.table("Mapped.EventsToExons.txt",
                             header=TRUE, sep="\t")
rMATS_res_all_filt <- rMATS_res_all %>% dplyr::select(GeneID, ID, DexseqFragment, geneSymbol, pvalue, padj, log2fold_CD8naive_Bnaive,
                                                      countData.case_Bnaive, countData.case_CD8naive, padj,
                                                      FDR, IncLevel1, IncLevel2, IncLevelDifference)
rMATS_res_all_filt <- rMATS_res_all_filt %>% filter(padj <= 0.05 & FDR <= 0.05)
rMATS_res_all_filt <- rMATS_res_all_filt %>% filter(abs(IncLevelDifference) >= 0.22)
rMATS_res_all_filt <- rMATS_res_all_filt %>% separate(ID, c("Event", "ID"))

rMATS_res_all_filt_final <- rMATS_res_all_filt %>% group_by(across(Event:ID)) %>% slice(which.min(padj))


rMATS_res_all_filt_final$DEXSeqID <- paste(rMATS_res_all_filt_final$geneSymbol, 
                                           rMATS_res_all_filt_final$DexseqFragment,
                                           rMATS_res_all_filt_final$ID,
                                           sep="_")
rMATS_res_all_filt_final <- rMATS_res_all_filt_final %>% 
  separate(IncLevel1, 
           c("IncLevel.Bnaive1", "IncLevel.Bnaive2","IncLevel.Bnaive3", "IncLevel.Bnaive4", "IncLevel.Bnaive5"), 
           sep=",")
rMATS_res_all_filt_final <- rMATS_res_all_filt_final %>% 
  separate(IncLevel2,
           c("IncLevel.CD8naive1", "IncLevel.CD8naive2","IncLevel.CD8naive3", "IncLevel.CD8naive4", "IncLevel.CD8naive5"),
           sep=",")


#A3SS
rMATS_A3SS <- rMATS_res_all_filt_final %>% 
  filter(Event == "A3SS") %>% 
  ungroup() %>% 
  dplyr::select(DEXSeqID, 12:21)
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
  dplyr::select(DEXSeqID, 12:21)
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
  dplyr::select(DEXSeqID, 12:21)
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
  dplyr::select(DEXSeqID, 12:21)
rMATS_SE <- as.data.frame(rMATS_SE)
rownames(rMATS_SE) <- rMATS_SE$DEXSeqID
rMATS_SE$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(rMATS_SE, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
rMATS_SE <- rMATS_SE[order(cutree(hclust_object, k = 2)), ]

rMATS_combined <- rbind(rMATS_A3SS, rMATS_A5SS, rMATS_RI, rMATS_SE)
rMATS_combined <- rMATS_combined %>% dplyr::select(IncLevel.CD8naive5, IncLevel.CD8naive4, IncLevel.CD8naive3, IncLevel.CD8naive2, IncLevel.CD8naive1,
                                                   IncLevel.Bnaive5, IncLevel.Bnaive4, IncLevel.Bnaive3, IncLevel.Bnaive2, IncLevel.Bnaive1)
rMATS_combined <- mutate_all(rMATS_combined, function(x) as.numeric(as.character(x)))
rMATS_combined_2 <- rMATS_combined
rMATS_combined_2$ID <- rownames(rMATS_combined) #counts with ID as a column 
rMATS_combined_2_tmp <- merge(rMATS_combined_2, rMATS_res_all_filt_final, 
                               by.x = "ID", by.y = "DEXSeqID") #counts with information from the full DEXSeq results table
rownames(rMATS_combined_2_tmp) <- rMATS_combined_2_tmp$ID
AS_row_rMATS <- rMATS_combined_2_tmp %>% dplyr::select(Event) #need this to label splicing events on heatmap
col_annot_rmats <- data.frame(samples = c(rep("CD8 naive", 5), rep("B naive", 5))) 
rownames(col_annot_rmats) <- colnames(rMATS_combined) #need this to label samples on heatmap


pdf("rmats_sig_exonic_heatmap.pdf", width=18, height=3.4)
pheatmap(t(rMATS_combined),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rmats,
         annotation_col = AS_row_rMATS,
         angle_col = "270")
dev.off()
##################################################################################################################
#obtaining heatmap with same IDs in both rMATS and DEXSeq

DEXSeq_sig_IDs <- data.frame(ID = rownames(dexseq_combined))
DEXSeq_sig_IDs <- DEXSeq_sig_IDs %>% separate(ID, c("geneSymbol", "DEXSeqExonFragment", "Event"), sep = "_")
rMATS_sig_IDs <- data.frame(ID = rownames(rMATS_combined))
rMATS_sig_IDs <- rMATS_sig_IDs %>% separate(ID, c("geneSymbol", "DEXSeqExonFragment", "rMATSID"), sep = "_")

sameIDs <- merge(DEXSeq_sig_IDs, rMATS_sig_IDs, by=c("geneSymbol", "DEXSeqExonFragment"))
sameIDs$DEXSeqID <- paste0(sameIDs$geneSymbol, "_", sameIDs$DEXSeqExonFragment, "_", sameIDs$Event)
sameIDs$rMATSID <- paste0(sameIDs$geneSymbol, "_", sameIDs$DEXSeqExonFragment, "_", sameIDs$rMATSID)

rMATS_sameIDs <- rMATS_combined[rownames(rMATS_combined) %in% sameIDs$rMATSID, ]
AS_row_rMATS_sameIDs <- subset(AS_row_rMATS, (rownames(AS_row_rMATS) %in% sameIDs$rMATSID))
col_annot_rMATS <- col_annot_rmats
rownames(col_annot_rMATS) <- colnames(rMATS_sameIDs)
pdf("rmats_sameIDs_heatmap.pdf", width=9, height=3.4)
pheatmap(t(rMATS_sameIDs),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rMATS,
         annotation_col = AS_row_rMATS_sameIDs,
         angle_col = "270")
dev.off()



rMATS_sameIDs$IDs <- rownames(rMATS_sameIDs)
rMATS_sameIDs <- rMATS_sameIDs %>% separate(IDs, c("geneSymbol", "DEXSeqExonFragment", "rMATSID"), sep = "_")
rMATS_sameIDs$DEXSeqID <- paste0(rMATS_sameIDs$geneSymbol, "_", rMATS_sameIDs$DEXSeqExonFragment)

dexseq_combined_2_tmp <- dexseq_combined_2_tmp %>% unite(col = "DEXSeqID", geneSymbol, featureID, sep="_")
dexseq_combined_2_tmp <- dexseq_combined_2_tmp[dexseq_combined_2_tmp$DEXSeqID %in% rMATS_sameIDs$DEXSeqID, ]
dexseq_sameIDs <- dexseq_combined_2_tmp %>% dplyr::select(2:11)
colnames(dexseq_sameIDs) <- colnames(dexseq_combined)

#match rows of rMATS sig to DEXSeq sig
dexseq_sameIDs$DEXSeqID <- rownames(dexseq_sameIDs)
dexseq_sameIDs <- dexseq_sameIDs %>% separate(DEXSeqID, c("geneSymbol", "DEXSeqExonFragment", "Event"))
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", geneSymbol, DEXSeqExonFragment, sep="_")
dexseq_sameIDs <- left_join(rMATS_sameIDs, dexseq_sameIDs, by = "DEXSeqID")
dexseq_sameIDs <- dexseq_sameIDs[,14:25]
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", DEXSeqID, Event, sep="_")
dexseq_sameIDs <- unique(dexseq_sameIDs)
rownames(dexseq_sameIDs) <- dexseq_sameIDs$DEXSeqID
dexseq_sameIDs$DEXSeqID <- NULL

AS_row_dexseq_sameIDs <- subset(AS_row_dexseq, (rownames(AS_row_dexseq) %in% rownames(dexseq_sameIDs)))
col_annot_dexseq <- col_annot
rownames(col_annot_dexseq) <- colnames(dexseq_sameIDs)
pdf("dexseq_sameIDs_heatmap.pdf", width=8, height=3)
pheatmap(t(dexseq_sameIDs),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_sameIDs,
         angle_col = "270")
dev.off()

####################################################################################################################
#DEXSeq only significant exons

dexseq_only <- dexseq_combined[!(rownames(dexseq_combined) %in% sameIDs$DEXSeqID), ]
AS_row_dexseq_only <- subset(AS_row_dexseq, !(rownames(AS_row_dexseq) %in% rownames(AS_row_dexseq_sameIDs)))
pdf("dexseq_only_heatmap.pdf", width=10, height=3)
pheatmap(t(dexseq_only),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_only,
         #annotation_colors = list(
           #samples = c(`B naive`="violet", `CD8 naive`="salmon"),
           #Event = c(A3SS="darkkhaki", RI="cornflowerblue", SE="green3")),
         angle_col = "270")
dev.off()

#rMATS only significant events

rMATS_only <- rMATS_combined[!(rownames(rMATS_combined) %in% sameIDs$rMATSID), ]
AS_row_rMATS_only <- subset(AS_row_rMATS, !(rownames(AS_row_rMATS) %in% rownames(AS_row_rMATS_sameIDs)))
pdf("rmats_only_heatmap.pdf", width=10, height=3.4)
pheatmap(t(rMATS_only),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_rMATS,
         annotation_col = AS_row_rMATS_only,
         angle_col = "270")
dev.off()

########################################################################################################################
#using this/others counts instead of this counts
this_others <- read.table("/mnt/storage/jaquino/GrASE_results_bonnal/this_others_counts.txt",
                          header=TRUE, sep="\t")
this_others <- this_others %>% mutate(Bnaive1 = V1/V11,
                                      Bnaive2 = V2/V12,
                                      Bnaive3 = V3/V13,
                                      Bnaive4 = V4/V14,
                                      Bnaive5 = V5/V15,
                                      CD8naive1 = V6/V16,
                                      CD8naive2 = V7/V17,
                                      CD8naive3 = V8/V18,
                                      CD8naive4 = V9/V19,
                                      CD8naive5 = V10/V20)
this_others <- this_others[,21:30]
this_others$DEXSeqID <- rownames(this_others)
this_others <- this_others %>% separate(DEXSeqID, c("groupID", "featureID"), sep=":")
this_others_merged <- merge(dexseq_res_all_filt_final, this_others, by = c("groupID", "featureID"))

#A3SS
dexseq_A3SS <- this_others_merged %>% filter(Event == "A3SS") %>% dplyr::select(DEXSeqID, 24:33)
rownames(dexseq_A3SS) <- dexseq_A3SS$DEXSeqID
dexseq_A3SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_A3SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_A3SS <- dexseq_A3SS[order(cutree(hclust_object, k = 2)), ]

#A5SS
dexseq_A5SS <- this_others_merged %>% filter(Event == "A5SS") %>% dplyr::select(DEXSeqID, 24:33)
rownames(dexseq_A5SS) <- dexseq_A5SS$DEXSeqID
dexseq_A5SS$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_A5SS, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_A5SS <- dexseq_A5SS[order(cutree(hclust_object, k = 2)), ]

#RI
dexseq_RI <- this_others_merged %>% filter(Event == "RI") %>% dplyr::select(DEXSeqID, 24:33)
rownames(dexseq_RI) <- dexseq_RI$DEXSeqID
dexseq_RI$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_RI, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_RI <- dexseq_RI[order(cutree(hclust_object, k = 2)), ]

#SE
dexseq_SE <- this_others_merged %>% filter(Event == "SE") %>% dplyr::select(DEXSeqID, 24:33)
rownames(dexseq_SE) <- dexseq_SE$DEXSeqID
dexseq_SE$DEXSeqID <- NULL
# Perform hierarchical clustering on rows (genes)
dist_matrix <- dist(dexseq_SE, method = "euclidean")
hclust_object <- hclust(dist_matrix, method = "complete")
dexseq_SE <- dexseq_SE[order(cutree(hclust_object, k = 2)), ]

dexseq_combined <- rbind(dexseq_A3SS, dexseq_A5SS, dexseq_RI, dexseq_SE)
#dexseq_combined <- dexseq_combined %>% dplyr::select(countData.CD8naive1, countData.CD8naive2, countData.CD8naive3, countData.CD8naive4, countData.CD8naive5,
#                                                     countData.Bnaive1, countData.Bnaive2, countData.Bnaive3, countData.Bnaive4, countData.Bnaive5) 
dexseq_combined[dexseq_combined==0] <- .00001 #set 0s to very small number so you can use the log2 function
dexseq_combined <- mutate_all(dexseq_combined, function(x) as.numeric(as.character(x)))
dexseq_combined <- log2(dexseq_combined)

dexseq_combined_2 <- dexseq_combined
dexseq_combined_2$ID <- rownames(dexseq_combined) #counts with ID as a column 
dexseq_combined_2_tmp <- merge(dexseq_combined_2, dexseq_res_all_filt_final, 
                               by.x = "ID", by.y = "DEXSeqID") #counts with information from the full DEXSeq results table
rownames(dexseq_combined_2_tmp) <- dexseq_combined_2_tmp$ID
AS_row_dexseq <- dexseq_combined_2_tmp %>% dplyr::select(Event) #need this to label splicing events on heatmap
col_annot <- data.frame(samples = c(rep("B naive", 5), rep("CD8 naive", 5))) 
rownames(col_annot) <- colnames(dexseq_combined) #need this to label samples on heatmap


pdf("dexseq_sig_exonic_heatmap_this_others.pdf", width=18, height=3.4)
pheatmap(t(dexseq_combined),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot,
         annotation_col = AS_row_dexseq,
         angle_col = "270")
dev.off()

dexseq_sameIDs <- dexseq_combined[(rownames(dexseq_combined) %in% sameIDs$DEXSeqID), ]
dexseq_sameIDs$DEXSeqID <- rownames(dexseq_sameIDs)
dexseq_sameIDs <- dexseq_sameIDs %>% separate(DEXSeqID, c("geneSymbol", "DEXSeqExonFragment", "Event"))
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", geneSymbol, DEXSeqExonFragment, sep="_")
dexseq_sameIDs <- left_join(rMATS_sameIDs, dexseq_sameIDs, by = "DEXSeqID")
dexseq_sameIDs <- dexseq_sameIDs[,14:25]
dexseq_sameIDs <- dexseq_sameIDs %>% unite(col = "DEXSeqID", DEXSeqID, Event, sep="_")
dexseq_sameIDs <- unique(dexseq_sameIDs)
rownames(dexseq_sameIDs) <- dexseq_sameIDs$DEXSeqID
dexseq_sameIDs$DEXSeqID <- NULL

AS_row_dexseq_sameIDs <- subset(AS_row_dexseq, (rownames(AS_row_dexseq) %in% rownames(dexseq_sameIDs)))
col_annot_dexseq <- col_annot
rownames(col_annot_dexseq) <- colnames(dexseq_sameIDs)
pdf("dexseq_sameIDs_heatmap_this_others.pdf", width=8, height=3)
pheatmap(t(dexseq_sameIDs),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_sameIDs,
         angle_col = "270")
dev.off()

dexseq_only <- dexseq_combined[!(rownames(dexseq_combined) %in% sameIDs$DEXSeqID), ]
AS_row_dexseq_only <- subset(AS_row_dexseq, !(rownames(AS_row_dexseq) %in% rownames(AS_row_dexseq_sameIDs)))
pdf("dexseq_only_heatmap_this_others.pdf", width=10, height=3)
pheatmap(t(dexseq_only),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = col_annot_dexseq,
         annotation_col = AS_row_dexseq_only,
         #annotation_colors = list(
           #samples = c(`B naive`="violet", `CD8 naive`="salmon"),
           #Event = c(A3SS="darkkhaki", RI="cornflowerblue", SE="green3")),
         angle_col = "270")
dev.off()
