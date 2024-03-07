dexseq_res_all <- read.table("/mnt/storage/jaquino/GrASE_results_bonnal/grase_results_v2/Mapped.ExonsToEvents.txt",
                             header=TRUE, sep="\t")

#All of DEXSeq significant exons = 3647
dexseq_all_sig <- dexseq_res_all %>% 
  filter(padj <= 0.05) %>% 
  group_by(groupID, featureID) %>% 
  slice(which.min(padj))

#DEXSeq sig, exon has AS event, has rMATS FDR, but not sig in rMATS = 616
dexseq_sig_rMATS_tested_not_sig <- dexseq_all_sig %>% 
  filter(!(is.na(rMATS_ID)) | rMATS_ID != "") %>% 
  filter(FDR > 0.05) 

#DEXSeq sig, exon has AS event, sig in rMATS = 177
dexseq_sig_rMATS_tested_sig <- dexseq_all_sig %>% 
  filter(!(is.na(rMATS_ID)) | rMATS_ID != "") %>% 
  filter(FDR <= 0.05) 

#DEXSeq sig, has AS event, no rMATS FDR = 124
dexseq_sig_rMATS_detected_not_sig <- dexseq_all_sig %>% 
  filter(!(is.na(rMATS_ID)), rMATS_ID != "") %>% 
  filter(is.na(FDR))

#rMATS sig = 7432
rMATS_sig <- dexseq_res_all %>% 
  filter(FDR <= 0.05) %>% 
  group_by(groupID, featureID) %>% 
  slice(which.min(FDR))

rMATS_res_all <- read.table("/mnt/storage/jaquino/GrASE_results_bonnal/grase_results_v2/Mapped.EventsToExons.txt",
                             header=TRUE, sep="\t")
#rMATS detected
length(unique(rMATS_res_all$ID)) #135880

#rMATS tested
rMATS_tested <- rMATS_res_all %>% filter(FDR <= 1) 
length(unique(rMATS_tested$ID)) #92063

#rMATS sig
rMATS_sig <- rMATS_tested %>% filter(FDR <= 0.05)
length(unique(rMATS_sig$ID)) #3504

#rMATS sig, DEXSeq sig
rMATS_sig_DEXSeq_sig <- rMATS_sig %>% filter(padj <= 0.05)
length(unique(rMATS_sig_DEXSeq_sig$ID)) #268

#rMATS tested, not sig, but DEXSeq sig
rMATS_tested_DEXSeq_sig <- rMATS_tested %>% filter(FDR > 0.05 & padj <= 0.05)
length(unique(rMATS_tested_DEXSeq_sig$ID)) #1033

#rMATS detected, no FDR, but DEXSeq sig
rMATS_detected_DEXSeq_sig <- rMATS_res_all %>% 
  filter(is.na(FDR)) %>%
  filter(padj <=0.05)
length(unique(rMATS_detected_DEXSeq_sig$ID)) #172

#DEXSeq tested after filtering #133970
DEXSeq_filtered_res <- read.table("/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/bonnal_dexseq_results_138997.txt",
                                  header=TRUE, sep="\t")

DEXSeq_filtered_res_rMATS_ID <- merge(DEXSeq_filtered_res, dexseq_res_all, by = c("groupID", "featureID"))
DEXSeq_filtered_res_rMATS_ID$DEXSeqID <- paste0(DEXSeq_filtered_res_rMATS_ID$groupID,"_",DEXSeq_filtered_res_rMATS_ID$featureID)
length(unique(DEXSeq_filtered_res_rMATS_ID$DEXSeqID)) #133970 exons with AS events

DEXSeq_all_sig <- DEXSeq_filtered_res_rMATS_ID %>% 
  filter(padj.x <= 0.05) %>% 
  group_by(groupID, featureID) %>% 
  slice(which.min(padj.x))
length(unique(DEXSeq_all_sig$DEXSeqID)) #645

#DEXSeq sig, exon has AS event, has rMATS FDR, but not sig in rMATS 
dexseq_sig_rMATS_tested_not_sig <- DEXSeq_all_sig %>% 
  filter(!(is.na(rMATS_ID)) | rMATS_ID != "") %>% 
  filter(FDR > 0.05) 
length(unique(dexseq_sig_rMATS_tested_not_sig$DEXSeqID)) #422

#DEXSeq sig, exon has AS event, sig in rMATS 
dexseq_sig_rMATS_tested_sig <- DEXSeq_all_sig %>% 
  filter(!(is.na(rMATS_ID)) | rMATS_ID != "") %>% 
  filter(FDR <= 0.05) 
length(unique(dexseq_sig_rMATS_tested_sig$DEXSeqID)) #139

#DEXSeq sig, has AS event, no rMATS FDR 
dexseq_sig_rMATS_detected_not_sig <- DEXSeq_all_sig %>% 
  filter(!(is.na(rMATS_ID)), rMATS_ID != "") %>% 
  filter(is.na(FDR))
length(unique(dexseq_sig_rMATS_detected_not_sig$DEXSeqID)) #84

#rMATS sig
rMATS_sig <- DEXSeq_filtered_res_rMATS_ID %>% 
  filter(FDR <= 0.05) %>% 
  group_by(groupID, featureID) %>% 
  slice(which.min(FDR))
length(unique(rMATS_sig$DEXSeqID)) #7123

#rMATS tested = 91738
rMATS_tested <- DEXSeq_filtered_res_rMATS_ID %>% 
  filter(FDR <= 1) %>% 
  group_by(groupID, featureID) %>% 
  slice(which.min(FDR))


