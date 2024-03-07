library(tidyverse)
#GrASE results: exploded exons mapped to AS events
dexseq_res_all <- read.table("/mnt/storage/jaquino/GrASE_results_bonnal/grase_results_v2/Mapped.ExonsToEvents.txt",
                             header=TRUE, sep="\t")

#exons with rMATS events
dexseq_res_all_filt_before <- dexseq_res_all[!(is.na(dexseq_res_all$rMATS_ID) | 
                                                 dexseq_res_all$rMATS_ID==""), ]

#SE exons before filtering
dexseq_res_all_filt_before_SE <- dexseq_res_all_filt_before %>% filter(str_detect(rMATS_ID, "SE"))
dexseq_res_all_filt_before_SE_sig <- dexseq_res_all_filt_before_SE %>% filter(padj<=0.05)
#1125/165064 = 0.006815538 significant exons

#RI exons before filtering
dexseq_res_all_filt_before_RI <- dexseq_res_all_filt_before %>% filter(str_detect(rMATS_ID, "RI"))
dexseq_res_all_filt_before_RI_sig <- dexseq_res_all_filt_before_RI %>% filter(padj<=0.05)
#436/46385 = 0.00939959 significant exons

#A3SS exons before filtering
dexseq_res_all_filt_before_A3SS <- dexseq_res_all_filt_before %>% filter(str_detect(rMATS_ID, "A3SS"))
dexseq_res_all_filt_before_A3SS_sig <- dexseq_res_all_filt_before_A3SS %>% filter(padj<=0.05)
#163/19369 = 0.008415509 significant exons

#A5SS exons before filtering
dexseq_res_all_filt_before_A5SS <- dexseq_res_all_filt_before %>% filter(str_detect(rMATS_ID, "A5SS"))
dexseq_res_all_filt_before_A5SS_sig <- dexseq_res_all_filt_before_A5SS %>% filter(padj<=0.05)
#161/14279 = 0.0112753 significant exons

#DEXSeq results after filtering
dexseq_res_after <- read.table("/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/bonnal_dexseq_results_138997.txt", header=TRUE, sep="\t")
dexseq_res_after_sig_merge <- merge(dexseq_res_all, dexseq_res_after, by=c("groupID","featureID"))

#SE exons after filtering
dexseq_res_all_filt_after_SE <- dexseq_res_after_sig_merge %>% filter(str_detect(rMATS_ID, "SE"))
dexseq_res_all_filt_after_SE_sig <- dexseq_res_all_filt_after_SE %>% filter(padj.y<=0.05)
#772/159592 = 0.004837335 significant exons

#RI exons after filtering
dexseq_res_all_filt_after_RI <- dexseq_res_after_sig_merge %>% filter(str_detect(rMATS_ID, "RI"))
dexseq_res_all_filt_after_RI_sig <- dexseq_res_all_filt_after_RI %>% filter(padj.y<=0.05)
#345/44089 = 0.007825081 significant exons

#A3SS exons after filtering
dexseq_res_all_filt_after_A3SS <- dexseq_res_after_sig_merge %>% filter(str_detect(rMATS_ID, "A3SS"))
dexseq_res_all_filt_after_A3SS_sig <- dexseq_res_all_filt_after_A3SS %>% filter(padj.y<=0.05)
#128/18464 = 0.006932409 significant exons

#A5SS exons after filtering
dexseq_res_all_filt_after_A5SS <- dexseq_res_after_sig_merge %>% filter(str_detect(rMATS_ID, "A5SS"))
dexseq_res_all_filt_after_A5SS_sig <- dexseq_res_all_filt_after_A5SS %>% filter(padj.y<=0.05)
#121/13755 = 0.008796801 significant exons

library(reshape2)
df_barPlot <- data.frame(`AS Event` = c("A3SS", "A5SS", "RI", "SE"),
                         Before = c(0.008415509, 0.0112753, 0.00939959, 0.006815538),
                         After = c(0.006932409, 0.008796801, 0.007825081, 0.004837335))
dfm <- melt(df_barPlot[,c('AS.Event','Before','After')],id.vars = 1)

pdf("DEXSeq_before_after_filt_barplot.pdf", width=5.99, height=5.99)
ggplot(dfm,aes(x = AS.Event,y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") + 
  scale_fill_manual(values=c('lightpink1', 'darkseagreen2')) + 
  ylab("Proportion DEXSeq Successful") +
  xlab("AS Event") +
  theme_classic() +
  labs(fill=NULL)
dev.off()