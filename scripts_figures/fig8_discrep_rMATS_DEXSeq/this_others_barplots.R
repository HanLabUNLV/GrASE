library(dplyr)
dexseq_res_all <- read.table("Mapped.ExonsToEvents.txt",
                             header=TRUE, sep="\t")
dexseq_res_all$transcripts <- lapply(dexseq_res_all$transcripts, function(x) unlist(strsplit(x, " "))) #unlist the transcript column
dexseq_res_all$tx_count <- sapply(dexseq_res_all$transcripts, length) #obtain how many transcripts the exon is part of

#dexseq this
tx_count <- sort(unique(dexseq_res_all$tx_count))

sig_results <- c()
norm_tx_count <- c()
for(i in tx_count) { 
  dexseq_res_all_filt <- dexseq_res_all %>% filter(tx_count == i)
  dexseq_sig <- dexseq_res_all_filt %>% filter(padj <= 0.05)
  sig_results <- c(sig_results, dim(dexseq_sig)[1])
  norm_tx_count <- c(norm_tx_count, dim(dexseq_res_all_filt)[1])
}

tx_count_dexseq_sig_this <- data.frame(tx.cnt = tx_count,
                                       sig.cnt = sig_results,
                                       norm.sig.cnt = sig_results/norm_tx_count,
                                       software = "DEXSeq")

#rmats this
rmats_res_all <- read.table("Mapped.EventsToExons.txt",
                            header=TRUE, sep="\t")
rmats_res_all$transcripts <- strsplit(rmats_res_all$transcripts, " ")
rmats_res_all$tx_count <- sapply(rmats_res_all$transcripts, length)

tx_count <- unique(rmats_res_all$tx_count)
sig_results <- c()
norm_tx_count <- c()
for(i in tx_count) { #bin_count
  rmats_res_all_filt <- rmats_res_all %>% filter(tx_count == i) # change to bin if using bin
  rmats_sig <- rmats_res_all_filt %>% filter(FDR <= 0.05)
  sig_results <- c(sig_results, dim(rmats_sig)[1])
  norm_tx_count <- c(norm_tx_count, dim(rmats_res_all_filt)[1])
}

tx_count_rmats_sig_this <- data.frame(tx.cnt = tx_count, # bin_count
                                      sig.cnt = sig_results,
                                      norm.sig.cnt = sig_results/norm_tx_count,
                                      software = "rMATS")

tx_count_both_sig_this <- rbind(tx_count_dexseq_sig_this, tx_count_rmats_sig_this)

#barchart
pdf("this_barchart_both.pdf", width=10, height=7) 
ggplot(tx_count_both_sig_this, aes(x=tx.cnt, y=norm.sig.cnt, fill=software)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.2), alpha = 0.7)+
  theme_classic(base_size = 22) +
  labs(title= "Number of transcripts this exon is part of vs\nNumber of significant events (this)",
       x = "transcript counts this exon is part of", 
       y = "normalized no. significant events") +
  #scale_x_continuous(seq(0, 100, by = 10), limits = c(0, 100)) +
  xlim(0,60) +
  scale_fill_manual(values = c("darkgreen", "lightpink1")) +
  geom_vline(xintercept = 30.5, linetype=3, size=1.5)
dev.off()

#obtain differences between rMATS and DEXSeq on the proportions of significant event for each transcript 

diff_rMATS_minus_DEXSeq <- c() 
for(i in sort(unique(tx_count_both_sig_this$tx.cnt))) { #for each unique transcript counts of each gene
  filt_df <- filter(tx_count_both_sig_this, tx_count_both_sig_this$tx.cnt == i) #obtain only that tx count
  if(dim(filt_df)[1]==2){ #there will only be two rows in the table (one for rMATS another for DEXSeq)
    diff <- filt_df$norm.sig.cnt[2]-filt_df$norm.sig.cnt[1] #subtract rMATS minus DEXSeq normalized counts
    diff_rMATS_minus_DEXSeq <- c(diff_rMATS_minus_DEXSeq, diff) # add the difference to the list
  }
  else { # there will be cases where there is only one transcript count observed in DEXSeq or rMATS
    if(filt_df$software == "DEXSeq") {
      diff <- 0 - filt_df$norm.sig.cnt # there will be 0 for rMATS when you only have a value for DEXSeq
      diff_rMATS_minus_DEXSeq <- c(diff_rMATS_minus_DEXSeq, diff)
    }
    else {
      diff <- filt_df$norm.sig.cnt - 0 # there will be 0 for DEXSeq when you only have a value for rMATS
      diff_rMATS_minus_DEXSeq <- c(diff_rMATS_minus_DEXSeq, diff) 
    }
  }
}
diff_df_this <- data.frame(tx.cnt = sort(unique(tx_count_both_sig_this$tx.cnt)),
                           diff = diff_rMATS_minus_DEXSeq) 
#diff_df_this <- diff_df_this[diff_df_this$diff != 0, ] # remove zero rows
diff_30_or_less_tx <- diff_df_this[1:30,] 
diff_30_or_more_tx <- diff_df_this[31:55,] #diff_df_this[31:length(diff_df_this$tx.cnt),]
t.test(diff_30_or_less_tx$diff, diff_30_or_more_tx$diff, alternative = 'less') #one sided t-test
wilcox.test(diff_30_or_less_tx$diff, diff_30_or_more_tx$diff) #perform wilcoxon ranksum t-test


##################################################################################################################
#dexseq others
#others
dexseq_res_all_filt <- dexseq_res_all %>% dplyr::select(groupID, featureID, padj, transcripts)
geneIDs <- unique(dexseq_res_all$groupID)
genes_tx_cnts <- c()
for(i in geneIDs){
  dexseq_res_all_filt_2 <- dexseq_res_all_filt %>% filter(groupID == i)
  gene_tx <- unique(unlist(dexseq_res_all_filt_2$transcripts))
  genes_tx_cnts <- c(genes_tx_cnts, length(gene_tx))
}
genes_tx_cnt <- data.frame(geneID = geneIDs,
                           tx.cnt = genes_tx_cnts)

dexseq_res_all_filt <- merge(dexseq_res_all_filt, genes_tx_cnt, by.x="groupID", by.y="geneID")
tx_count <- unique(dexseq_res_all_filt$tx.cnt)

sig_results <- c()
norm_tx_count <- c()
dexseq_total <- c()
for(i in tx_count) {
  dexseq_res_all_filt_2 <- dexseq_res_all_filt %>% filter(tx.cnt == i)
  dexseq_sig <- dexseq_res_all_filt_2 %>% filter(padj <= 0.05)
  sig_results <- c(sig_results, dim(dexseq_sig)[1])
  norm_tx_count <- c(norm_tx_count, dim(dexseq_res_all_filt_2)[1])
}

tx_count_dexseq_sig_others <- data.frame(tx.cnt = tx_count, 
                                         sig.cnt = sig_results,
                                         total.cnt = norm_tx_count,
                                         norm.sig.cnt = sig_results/norm_tx_count, 
                                         software = "DEXSeq")

#rmats others
rmats_res_all_filt <- merge(rmats_res_all, genes_tx_cnt, by.x = "GeneID", by.y = "geneID")
tx_count <- unique(rmats_res_all_filt$tx.cnt)

sig_results <- c()
norm_tx_count <- c()
for(i in tx_count) {
  rmats_res_all_filt_2 <- rmats_res_all_filt %>% filter(tx.cnt == i)
  rmats_sig <- rmats_res_all_filt_2 %>% filter(FDR <= 0.05)
  sig_results <- c(sig_results, dim(rmats_sig)[1])
  norm_tx_count <- c(norm_tx_count, dim(rmats_res_all_filt_2)[1])
}

tx_count_rmats_sig_others <- data.frame(tx.cnt = tx_count, 
                                        sig.cnt = sig_results,
                                        total.cnt = norm_tx_count,
                                        norm.sig.cnt = sig_results/norm_tx_count,
                                        software = "rMATS")

tx_count_both_sig_others <- rbind(tx_count_dexseq_sig_others, tx_count_rmats_sig_others)
pdf("others_barchart_both.pdf", width=10, height=7) 
ggplot(tx_count_both_sig_others, aes(x=tx.cnt, y=norm.sig.cnt, fill=software)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.2), alpha = 0.7)+
  theme_classic(base_size = 22) +
  #scale_x_continuous(seq(0, 105, by = 5), limits = c(0, 105)) +
  xlim(0,105) +
  scale_fill_manual(values = c("darkgreen", "lightpink1")) +
  labs(title= "Number of transcripts in a gene vs\nNumber of significant events (Others)",
       x = "transcript counts of genes", 
       y = "normalized no. significant events")
dev.off()

diff_rMATS_minus_DEXSeq <- c() 
for(i in sort(unique(tx_count_both_sig_others$tx.cnt))) { #for each unique transcript counts of each gene
  filt_df <- filter(tx_count_both_sig_others, tx_count_both_sig_others$tx.cnt == i) #obtain only that tx count
  if(dim(filt_df)[1]==2){ #there will only be two rows in the table (one for rMATS another for DEXSeq)
    diff <- filt_df$norm.sig.cnt[2]-filt_df$norm.sig.cnt[1] #subtract rMATS minus DEXSeq normalized counts
    diff_rMATS_minus_DEXSeq <- c(diff_rMATS_minus_DEXSeq, diff) # add the difference to the list
  }
  else { # there will be cases where there is only one transcript count observed in DEXSeq or rMATS
    if(filt_df$software == "DEXSeq") {
      diff <- 0 - filt_df$norm.sig.cnt # there will be 0 for rMATS when you only have a value for DEXSeq
      diff_rMATS_minus_DEXSeq <- c(diff_rMATS_minus_DEXSeq, diff)
    }
    else {
      diff <- filt_df$norm.sig.cnt - 0 # there will be 0 for DEXSeq when you only have a value for rMATS
      diff_rMATS_minus_DEXSeq <- c(diff_rMATS_minus_DEXSeq, diff) 
    }
  }
}

diff_df_others <- data.frame(tx.cnt = sort(unique(tx_count_both_sig_others$tx.cnt)),
                             diff = diff_rMATS_minus_DEXSeq) 
diff_df_others <- diff_df_others[diff_df_others$diff != 0, ] #remove 0 difference
diff_30_or_less_tx <- diff_df_others[1:50,] #genes with 30 or less transcripts
diff_30_or_more_tx <- diff_df_others[51:length(diff_df_others$diff),] #genes with more than 30 transcripts
diff_30_or_more_tx <- diff_df_others[51:71,]
t.test(diff_30_or_less_tx$diff, diff_30_or_more_tx$diff, alternative = 'less') 

