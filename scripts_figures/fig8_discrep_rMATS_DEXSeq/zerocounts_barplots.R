library(tidyverse)

rmats_res_all <- read.table("Mapped.EventsToExons.txt",
                            header=TRUE, sep="\t")

###############################################################################################################
#rMATS zero
#sig in both dexseq and rmats
dexseq_sig_rmats_sig <- rmats_res_all %>% filter(padj <= 0.05 & FDR <=0.05)
dexseq_sig_rmats_sig_tmp <- dexseq_sig_rmats_sig %>%
  group_by(across(ID)) %>% 
  summarise(padj = min(padj), .groups = 'drop')

dexseq_sig_rmats_sig <- merge(dexseq_sig_rmats_sig_tmp, dexseq_sig_rmats_sig, by=c("ID", "padj"))
dexseq_sig_rmats_sig$IJC_SAMPLE_1_sum <- sapply(dexseq_sig_rmats_sig$IJC_SAMPLE_1, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_sig$SJC_SAMPLE_1_sum <- sapply(dexseq_sig_rmats_sig$SJC_SAMPLE_1, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_sig$IJC_SAMPLE_2_sum <- sapply(dexseq_sig_rmats_sig$IJC_SAMPLE_2, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_sig$SJC_SAMPLE_2_sum <- sapply(dexseq_sig_rmats_sig$SJC_SAMPLE_2, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))


dexseq_sig_rmats_sig_filt <- dexseq_sig_rmats_sig %>% filter(IJC_SAMPLE_1_sum == 0 | SJC_SAMPLE_1_sum ==0 |
                                                               IJC_SAMPLE_2_sum == 0 | SJC_SAMPLE_2_sum ==0)
#14/268 have at least 1 zero

#sig in dexseq, but not in rmats
dexseq_sig_rmats_not <- rmats_res_all %>% filter(padj <= 0.05 & FDR > 0.05)
dexseq_sig_rmats_not_tmp <- dexseq_sig_rmats_not %>%
  group_by(across(ID)) %>% 
  summarise(padj = min(padj), .groups = 'drop')

dexseq_sig_rmats_not <- merge(dexseq_sig_rmats_not_tmp, dexseq_sig_rmats_not, by=c("ID", "padj"))
dexseq_sig_rmats_not$IJC_SAMPLE_1_sum <- sapply(dexseq_sig_rmats_not$IJC_SAMPLE_1, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_not$SJC_SAMPLE_1_sum <- sapply(dexseq_sig_rmats_not$SJC_SAMPLE_1, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_not$IJC_SAMPLE_2_sum <- sapply(dexseq_sig_rmats_not$IJC_SAMPLE_2, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_not$SJC_SAMPLE_2_sum <- sapply(dexseq_sig_rmats_not$SJC_SAMPLE_2, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))

dexseq_sig_rmats_not_filt <- dexseq_sig_rmats_not %>% filter(IJC_SAMPLE_1_sum == 0 | SJC_SAMPLE_1_sum ==0 |
                                                               IJC_SAMPLE_2_sum == 0 | SJC_SAMPLE_2_sum ==0)
#442/1033 have at least 1 zero

#sig in rmats, but not in dexseq
dexseq_not_rmats_sig <- rmats_res_all %>% filter(padj > 0.05 & FDR <= 0.05)
dexseq_not_rmats_sig_tmp <- dexseq_not_rmats_sig %>%
  group_by(across(ID)) %>% 
  summarise(padj = min(padj), .groups = 'drop')

dexseq_not_rmats_sig <- merge(dexseq_not_rmats_sig_tmp, dexseq_not_rmats_sig, by=c("ID", "padj"))
dexseq_not_rmats_sig$IJC_SAMPLE_1_sum <- sapply(dexseq_not_rmats_sig$IJC_SAMPLE_1, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_not_rmats_sig$SJC_SAMPLE_1_sum <- sapply(dexseq_not_rmats_sig$SJC_SAMPLE_1, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_not_rmats_sig$IJC_SAMPLE_2_sum <- sapply(dexseq_not_rmats_sig$IJC_SAMPLE_2, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_not_rmats_sig$SJC_SAMPLE_2_sum <- sapply(dexseq_not_rmats_sig$SJC_SAMPLE_2, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))

dexseq_not_rmats_sig_filt <- dexseq_not_rmats_sig %>% filter(IJC_SAMPLE_1_sum == 0 | SJC_SAMPLE_1_sum ==0 |
                                                               IJC_SAMPLE_2_sum == 0 | SJC_SAMPLE_2_sum ==0)
#249/3209 have at least 1 zero

###############################################################################################################
#DEXSeq zeros
#sig in both dexseq and rmats
dexseq_sig_rmats_sig$DEXSeq_Bcell_sum <- sapply(dexseq_sig_rmats_sig$countData.case_Bnaive, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_sig$DEXSeq_Tcell_sum <- sapply(dexseq_sig_rmats_sig$countData.case_CD8naive, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_sig_filt_2 <- dexseq_sig_rmats_sig %>% filter(DEXSeq_Bcell_sum == 0 | DEXSeq_Tcell_sum == 0)

#1/268

#sig in dexseq, but not in rmats
dexseq_sig_rmats_not$DEXSeq_Bcell_sum <- sapply(dexseq_sig_rmats_not$countData.case_Bnaive, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_not$DEXSeq_Tcell_sum <- sapply(dexseq_sig_rmats_not$countData.case_CD8naive, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_sig_rmats_not_filt_2 <- dexseq_sig_rmats_not %>% filter(DEXSeq_Bcell_sum == 0 | DEXSeq_Tcell_sum == 0)

#12/1033


#sig in rmats, but not in dexseq
dexseq_not_rmats_sig$DEXSeq_Bcell_sum <- sapply(dexseq_not_rmats_sig$countData.case_Bnaive, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_not_rmats_sig$DEXSeq_Tcell_sum <- sapply(dexseq_not_rmats_sig$countData.case_CD8naive, 
                                                function(x) sum(as.numeric(unlist(str_split(x, ",")))))
dexseq_not_rmats_sig_filt_2 <- dexseq_not_rmats_sig %>% filter(DEXSeq_Bcell_sum == 0 | DEXSeq_Tcell_sum == 0)

#12/3209


compare_groups <- data.frame(groups = c(rep(c("sig DEXSeq & sig rMATS", 
                                              "not sig DEXSeq & sig rMATS", 
                                              "sig DEXSeq & not sig rMATS"), 2)),
                             one.zero = c(14/268, 249/3209, 442/1033, 1/268, 12/3209, 12/1033),
                             software = c(rep("rMATS", 3), rep("DEXSeq", 3)))

pdf("onezero_both.pdf", width=11, height=5) 
ggplot(compare_groups, aes(x=factor(groups, levels =c("sig DEXSeq & sig rMATS", 
                                                      "not sig DEXSeq & sig rMATS", 
                                                      "sig DEXSeq & not sig rMATS")), 
                           y=one.zero, fill=software)) +
  geom_bar(stat="identity", position = 'dodge') +
  theme_classic(base_size = 22) +
  labs(title= "At least one zero",
       x = "", 
       y = "normalized no. events") +
  ylim(0,1) +
  scale_fill_manual(values = c("darkgreen", "lightpink1")) +
  coord_flip()
dev.off()

compare_groups <- data.frame(groups = c(rep(c("sig DEXSeq & sig rMATS", 
                                              "not sig DEXSeq & sig rMATS", 
                                              "sig DEXSeq & not sig rMATS"), 2)),
                             one.zero = c(14/(268*4), 249/(3209*4), 442/(1033*4), 1/(2682*2), 12/(3209*2), 12/(1033*2)),
                             software = c(rep("rMATS", 3), rep("DEXSeq", 3)))

pdf("totalzero_both.pdf", width=11, height=5) 
ggplot(compare_groups, aes(x=factor(groups, levels =c("sig DEXSeq & sig rMATS", 
                                                      "not sig DEXSeq & sig rMATS", 
                                                      "sig DEXSeq & not sig rMATS")), 
                           y=one.zero, fill=software)) +
  geom_bar(stat="identity", position = 'dodge') +
  theme_classic(base_size = 22) +
  labs(title= "Total number of zeros",
       x = "", 
       y = "normalized no. events") +
  ylim(0,1) +
  scale_fill_manual(values = c("darkgreen", "lightpink1")) +
  coord_flip()
dev.off()

