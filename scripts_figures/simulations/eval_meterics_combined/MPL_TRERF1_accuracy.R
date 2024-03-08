#raw pvalue
sens_raw <- read.table("/mnt/storage/jaquino/polyester_sim_2/A3SS_MPL_TRERF1_res/Sens_pvalue_A3SS_MPL_TRERF1_disp3.txt", header=TRUE, sep="\t")
spec_raw <- read.table("/mnt/storage/jaquino/polyester_sim_2/A3SS_MPL_TRERF1_res/Spec_pvalue_A3SS_MPL_TRERF1_disp3.txt", header=TRUE, sep="\t")

MPL_sens <- sens_raw[1:112,1:8]
MPL_spec <- spec_raw %>%filter(gene == "MPL")
MPL_spec <- MPL_spec[1:14,1:8]

DEXSeq_acc <- c()
rMATS_jcec_acc <- c()
rMATS_jc_acc <- c()
for(FC in unique(MPL_sens$FoldChange)){
  MPL_sens_filt <- MPL_sens %>% filter(FoldChange == FC)
  #DEXSeq accuracy by fold change
  DEXSeq_acc_FC <- (MPL_sens_filt$DEXSeqTP + MPL_spec$DEXSeqTN)/
    (MPL_sens_filt$DEXSeqTP + MPL_spec$DEXSeqTN +
    MPL_sens_filt$DEXSeqFN + MPL_spec$DEXSeqFP)
  DEXSeq_acc <- c(DEXSeq_acc, DEXSeq_acc_FC)
  
  #rMATS JCEC accuracy by fold change
  rMATS_jcec_acc_FC <- (MPL_sens_filt$rMATS_jcecTP + MPL_spec$rMATS_jcecTN)/
    (MPL_sens_filt$rMATS_jcecTP + MPL_spec$rMATS_jcecTN +
       MPL_sens_filt$rMATS_jcecFN + MPL_spec$rMATS_jcecFP)
  rMATS_jcec_acc <- c(rMATS_jcec_acc, rMATS_jcec_acc_FC)
  
  #rMATS JC accuracy by fold change
  rMATS_jc_acc_FC <- (MPL_sens_filt$rMATS_jcTP + MPL_spec$rMATS_jcTN)/
    (MPL_sens_filt$rMATS_jcTP + MPL_spec$rMATS_jcTN +
       MPL_sens_filt$rMATS_jcFN + MPL_spec$rMATS_jcFP)
  rMATS_jc_acc <- c(rMATS_jc_acc, rMATS_jc_acc_FC)
}

MPL_accuracy_raw <- data.frame(ReadDepth = MPL_sens$ReadDepth,
                               FoldChange = MPL_sens$FoldChange,
                               DEXSeq_acc = DEXSeq_acc,
                               rMATS_jcec_acc = rMATS_jcec_acc,
                               rMATS_jc_acc = rMATS_jc_acc,
                               gene = "MPL")


TRERF1_sens <- sens_raw %>% filter(gene == "TRERF1")
TRERF1_sens <- TRERF1_sens[1:112, 1:8]
TRERF1_spec <- spec_raw %>%filter(gene == "TRERF1")
TRERF1_spec <- TRERF1_spec[1:14,1:8]

DEXSeq_acc <- c()
rMATS_jcec_acc <- c()
rMATS_jc_acc <- c()
for(FC in unique(TRERF1_sens$FoldChange)){
  TRERF1_sens_filt <- TRERF1_sens %>% filter(FoldChange == FC)
  #DEXSeq accuracy by fold change
  DEXSeq_acc_FC <- (TRERF1_sens_filt$DEXSeqTP + TRERF1_spec$DEXSeqTN)/
    (TRERF1_sens_filt$DEXSeqTP + TRERF1_spec$DEXSeqTN +
       TRERF1_sens_filt$DEXSeqFN + TRERF1_spec$DEXSeqFP)
  DEXSeq_acc <- c(DEXSeq_acc, DEXSeq_acc_FC)
  
  #rMATS JCEC accuracy by fold change
  rMATS_jcec_acc_FC <- (TRERF1_sens_filt$rMATS_jcecTP + TRERF1_spec$rMATS_jcecTN)/
    (TRERF1_sens_filt$rMATS_jcecTP + TRERF1_spec$rMATS_jcecTN +
       TRERF1_sens_filt$rMATS_jcecFN + TRERF1_spec$rMATS_jcecFP)
  rMATS_jcec_acc <- c(rMATS_jcec_acc, rMATS_jcec_acc_FC)
  
  #rMATS JC accuracy by fold change
  rMATS_jc_acc_FC <- (TRERF1_sens_filt$rMATS_jcTP + TRERF1_spec$rMATS_jcTN)/
    (TRERF1_sens_filt$rMATS_jcTP + TRERF1_spec$rMATS_jcTN +
       TRERF1_sens_filt$rMATS_jcFN + TRERF1_spec$rMATS_jcFP)
  rMATS_jc_acc <- c(rMATS_jc_acc, rMATS_jc_acc_FC)
}

TRERF1_accuracy_raw <- data.frame(ReadDepth = TRERF1_sens$ReadDepth,
                               FoldChange = TRERF1_sens$FoldChange,
                               DEXSeq_acc = DEXSeq_acc,
                               rMATS_jcec_acc = rMATS_jcec_acc,
                               rMATS_jc_acc = rMATS_jc_acc,
                               gene = "TRERF1")

accuracy_combine <- rbind(MPL_accuracy_raw, TRERF1_accuracy_raw)
accuracy_combine <- accuracy_combine %>% gather(var, val, 3:5)
accuracy_combine$var2 <- paste0(accuracy_combine$var, accuracy_combine$gene)

write.table(accuracy_combine, "/mnt/storage/jaquino/polyester_sim_2/A3SS_MPL_TRERF1_res/Acc_pvalue_A3SS_MPL_TRERF1_disp3.txt", 
            sep="\t", quote=FALSE, row.names = FALSE)

accuracy_combine$FoldChange <- factor(accuracy_combine$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))
ggplot(accuracy_combine, aes(group=var2, y=val, x=factor(ReadDepth))) +
  geom_line(aes(linetype=gene, color=var), alpha = 0.55, size =1, position = position_dodge(width = 0.2)) +
  labs(title="A3SS Accuracy raw p-value", x="Read Depth", y="Accuracy") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1) +
  theme_bw(base_size = 20) +
  facet_wrap(~FoldChange, nrow=2)
ggsave("/mnt/storage/jaquino/polyester_sim_2/A3SS_MPL_TRERF1_res/Acc_pvalue_lineplot_A3SS_MPL_TRERF1.pdf", width = 28, height = 12, dpi = 100, units = "in", device='pdf')

#adjusted p-value
sens_adj <- read.table("/mnt/storage/jaquino/polyester_sim_2/A3SS_MPL_TRERF1_res/Sens_pvalue_A3SS_MPL_TRERF1_disp3_adj.txt", header=TRUE, sep="\t")
spec_adj <- read.table("/mnt/storage/jaquino/polyester_sim_2/A3SS_MPL_TRERF1_res/Spec_pvalue_A3SS_MPL_TRERF1_disp3_adj.txt", header=TRUE, sep="\t")

MPL_sens <- sens_adj[1:112,1:8]
MPL_spec <- spec_adj %>%filter(gene == "MPL")
MPL_spec <- MPL_spec[1:14,1:8]

DEXSeq_acc <- c()
rMATS_jcec_acc <- c()
rMATS_jc_acc <- c()
for(FC in unique(MPL_sens$FoldChange)){
  MPL_sens_filt <- MPL_sens %>% filter(FoldChange == FC)
  #DEXSeq accuracy by fold change
  DEXSeq_acc_FC <- (MPL_sens_filt$DEXSeqTP + MPL_spec$DEXSeqTN)/
    (MPL_sens_filt$DEXSeqTP + MPL_spec$DEXSeqTN +
       MPL_sens_filt$DEXSeqFN + MPL_spec$DEXSeqFP)
  DEXSeq_acc <- c(DEXSeq_acc, DEXSeq_acc_FC)
  
  #rMATS JCEC accuracy by fold change
  rMATS_jcec_acc_FC <- (MPL_sens_filt$rMATS_jcecTP + MPL_spec$rMATS_jcecTN)/
    (MPL_sens_filt$rMATS_jcecTP + MPL_spec$rMATS_jcecTN +
       MPL_sens_filt$rMATS_jcecFN + MPL_spec$rMATS_jcecFP)
  rMATS_jcec_acc <- c(rMATS_jcec_acc, rMATS_jcec_acc_FC)
  
  #rMATS JC accuracy by fold change
  rMATS_jc_acc_FC <- (MPL_sens_filt$rMATS_jcTP + MPL_spec$rMATS_jcTN)/
    (MPL_sens_filt$rMATS_jcTP + MPL_spec$rMATS_jcTN +
       MPL_sens_filt$rMATS_jcFN + MPL_spec$rMATS_jcFP)
  rMATS_jc_acc <- c(rMATS_jc_acc, rMATS_jc_acc_FC)
}

MPL_accuracy_adj <- data.frame(ReadDepth = MPL_sens$ReadDepth,
                               FoldChange = MPL_sens$FoldChange,
                               DEXSeq_acc = DEXSeq_acc,
                               rMATS_jcec_acc = rMATS_jcec_acc,
                               rMATS_jc_acc = rMATS_jc_acc,
                               gene = "MPL")


TRERF1_sens <- sens_adj %>% filter(gene == "TRERF1")
TRERF1_sens <- TRERF1_sens[1:112, 1:8]
TRERF1_spec <- spec_adj %>%filter(gene == "TRERF1")
TRERF1_spec <- TRERF1_spec[1:14,1:8]

DEXSeq_acc <- c()
rMATS_jcec_acc <- c()
rMATS_jc_acc <- c()
for(FC in unique(TRERF1_sens$FoldChange)){
  TRERF1_sens_filt <- TRERF1_sens %>% filter(FoldChange == FC)
  #DEXSeq accuracy by fold change
  DEXSeq_acc_FC <- (TRERF1_sens_filt$DEXSeqTP + TRERF1_spec$DEXSeqTN)/
    (TRERF1_sens_filt$DEXSeqTP + TRERF1_spec$DEXSeqTN +
       TRERF1_sens_filt$DEXSeqFN + TRERF1_spec$DEXSeqFP)
  DEXSeq_acc <- c(DEXSeq_acc, DEXSeq_acc_FC)
  
  #rMATS JCEC accuracy by fold change
  rMATS_jcec_acc_FC <- (TRERF1_sens_filt$rMATS_jcecTP + TRERF1_spec$rMATS_jcecTN)/
    (TRERF1_sens_filt$rMATS_jcecTP + TRERF1_spec$rMATS_jcecTN +
       TRERF1_sens_filt$rMATS_jcecFN + TRERF1_spec$rMATS_jcecFP)
  rMATS_jcec_acc <- c(rMATS_jcec_acc, rMATS_jcec_acc_FC)
  
  #rMATS JC accuracy by fold change
  rMATS_jc_acc_FC <- (TRERF1_sens_filt$rMATS_jcTP + TRERF1_spec$rMATS_jcTN)/
    (TRERF1_sens_filt$rMATS_jcTP + TRERF1_spec$rMATS_jcTN +
       TRERF1_sens_filt$rMATS_jcFN + TRERF1_spec$rMATS_jcFP)
  rMATS_jc_acc <- c(rMATS_jc_acc, rMATS_jc_acc_FC)
}

TRERF1_accuracy_adj <- data.frame(ReadDepth = TRERF1_sens$ReadDepth,
                                  FoldChange = TRERF1_sens$FoldChange,
                                  DEXSeq_acc = DEXSeq_acc,
                                  rMATS_jcec_acc = rMATS_jcec_acc,
                                  rMATS_jc_acc = rMATS_jc_acc,
                                  gene = "TRERF1")

accuracy_combine <- rbind(MPL_accuracy_adj, TRERF1_accuracy_adj)
accuracy_combine <- accuracy_combine %>% gather(var, val, 3:5)
accuracy_combine$var2 <- paste0(accuracy_combine$var, accuracy_combine$gene)

write.table(accuracy_combine, "/mnt/storage/jaquino/polyester_sim_2/A3SS_MPL_TRERF1_res/Acc_pvalue_A3SS_MPL_TRERF1_disp3_adj.txt", 
            sep="\t", quote=FALSE, row.names = FALSE)

accuracy_combine$FoldChange <- factor(accuracy_combine$FoldChange, levels=c(-8, -4, -2, -1, 8, 4, 2, 1))
ggplot(accuracy_combine, aes(group=var2, y=val, x=factor(ReadDepth))) +
  geom_line(aes(linetype=gene, color=var), alpha = 0.55, size =1, position = position_dodge(width = 0.2)) +
  labs(title="A3SS Accuracy adjusted p-value", x="Read Depth", y="Accuracy") +
  scale_color_discrete(name = "Software", labels=c("DEXSeq", "rMATS JC", "rMATS JCEC")) +
  ylim(0,1) +
  theme_bw(base_size = 20) +
  facet_wrap(~FoldChange, nrow=2)
ggsave("/mnt/storage/jaquino/polyester_sim_2/A3SS_MPL_TRERF1_res/Acc_pvalue_lineplot_A3SS_MPL_TRERF1_adj.pdf", width = 28, height = 12, dpi = 100, units = "in", device='pdf')
