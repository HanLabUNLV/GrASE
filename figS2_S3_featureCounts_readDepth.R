library(DESeq2)
library(Rsubread)

setwd("/mnt/storage/jaquino/bonnal_2015_bulkRNAseq")

#create sample table
samples <- read.table("samples.txt", header=FALSE, stringsAsFactors = FALSE)
colnames(samples) <- c("sample", "condition")

#fastq file path
fastq1_CD8 <- list.files(path="CD8_naive", pattern="*1.fastq.gz", full.names=TRUE)
fastq2_CD8 <- list.files(path="CD8_naive", pattern="*2.fastq.gz", full.names=TRUE) 

fastq1_B <- list.files(path="B_naive", pattern="*1.fastq.gz", full.names=TRUE)
fastq2_B <- list.files(path="B_naive", pattern="*2.fastq.gz", full.names=TRUE) 

#Build index for human reference genome hg38
Rsubread::buildindex(basename="GRCh38",
                      reference="/mnt/storage/jaquino/ref/GRCh38.p13.genome.fa",
                      memory=10000)

#align B naive fastq; generate BAM files
Rsubread::align(index = "GRCh38",
                readfile1 = fastq1_B,
                readfile2 = fastq2_B,
                type = "rna",
                input_format = "gzFASTQ",
                output_format = "BAM",
                PE_orientation = "rf",
                nthreads = 24 )

#align CD8 naive fastq; generate BAM files
Rsubread::align(index = "GRCh38",
                readfile1 = fastq1_CD8,
                readfile2 = fastq2_CD8,
                type = "rna",
                input_format = "gzFASTQ",
                output_format = "BAM",
                PE_orientation = "rf",
                nthreads = 36 )

#obtain bam files path
bamfiles <- list.files(path="bamfiles", pattern="*.subread.BAM$", full.names = TRUE)
#run featurecounts
fc <- featureCounts(bamfiles, annot.ext = "/mnt/storage/jaquino/ref/gencode.v34.annotation.gtf",
                    isGTFAnnotationFile = TRUE, GTF.featureType = "gene", GTF.attrType = "gene_id", 
                    isPairedEnd = TRUE, nthreads=36) 

#reorder featurecounts counts table to match sample table
bam_file_colnames <- str_split(colnames(fc$counts), "_")
bam_file_colnames <- sapply(bam_file_colnames, "[[", 1)
reorder_idx <- match(samples$sample,bam_file_colnames)
fc$counts <- fc$counts[,reorder_idx]
samples$samplename <- colnames(fc$counts)

fc_cnts <- as.data.frame(fc$counts)

ERR431592_fc_cnts <- as.data.frame(fc_cnts$ERR431592_1.fastq.gz.subread.BAM)
colnames(ERR431592_fc_cnts) <- "NumReads"
ERR431592_fc_cnts$sampleID <- "ERR431592"
ERR431592_fc_cnts$celltype <- "T"

ERR431605_fc_cnts <- as.data.frame(fc_cnts$ERR431605_1.fastq.gz.subread.BAM)
colnames(ERR431605_fc_cnts) <- "NumReads"
ERR431605_fc_cnts$sampleID <- "ERR431605"
ERR431605_fc_cnts$celltype <- "T"

ERR431612_fc_cnts <- as.data.frame(fc_cnts$ERR431612_1.fastq.gz.subread.BAM)
colnames(ERR431612_fc_cnts) <- "NumReads"
ERR431612_fc_cnts$sampleID <- "ERR431612"
ERR431612_fc_cnts$celltype <- "T"

ERR431623_fc_cnts <- as.data.frame(fc_cnts$ERR431623_1.fastq.gz.subread.BAM)
colnames(ERR431623_fc_cnts) <- "NumReads"
ERR431623_fc_cnts$sampleID <- "ERR431623"
ERR431623_fc_cnts$celltype <- "T"

ERR431627_fc_cnts <- as.data.frame(fc_cnts$ERR431627_1.fastq.gz.subread.BAM)
colnames(ERR431627_fc_cnts) <- "NumReads"
ERR431627_fc_cnts$sampleID <- "ERR431627"
ERR431627_fc_cnts$celltype <- "T"

ERR431572_fc_cnts <- as.data.frame(fc_cnts$ERR431572_1.fastq.gz.subread.BAM)
colnames(ERR431572_fc_cnts) <- "NumReads"
ERR431572_fc_cnts$sampleID <- "ERR431572"
ERR431572_fc_cnts$celltype <- "B"

ERR431586_fc_cnts <- as.data.frame(fc_cnts$ERR431586_1.fastq.gz.subread.BAM)
colnames(ERR431586_fc_cnts) <- "NumReads"
ERR431586_fc_cnts$sampleID <- "ERR431586"
ERR431586_fc_cnts$celltype <- "B"

ERR431611_fc_cnts <- as.data.frame(fc_cnts$ERR431611_1.fastq.gz.subread.BAM)
colnames(ERR431611_fc_cnts) <- "NumReads"
ERR431611_fc_cnts$sampleID <- "ERR431611"
ERR431611_fc_cnts$celltype <- "B"

ERR431619_fc_cnts <- as.data.frame(fc_cnts$ERR431619_1.fastq.gz.subread.BAM)
colnames(ERR431619_fc_cnts) <- "NumReads"
ERR431619_fc_cnts$sampleID <- "ERR431619"
ERR431619_fc_cnts$celltype <- "B"

ERR431624_fc_cnts <- as.data.frame(fc_cnts$ERR431624_1.fastq.gz.subread.BAM)
colnames(ERR431624_fc_cnts) <- "NumReads"
ERR431624_fc_cnts$sampleID <- "ERR431624"
ERR431624_fc_cnts$celltype <- "B"

master_B_fc <- rbind(ERR431572_fc_cnts, ERR431586_fc_cnts, ERR431611_fc_cnts, ERR431619_fc_cnts, ERR431624_fc_cnts)
master_T_fc <- rbind(ERR431592_fc_cnts, ERR431605_fc_cnts, ERR431612_fc_cnts, ERR431623_fc_cnts, ERR431627_fc_cnts)

p1 <- ggplot(master_B_fc, aes(x = NumReads)) +
  geom_histogram(aes(y = ..density..),
                 fill = "white", colour = "black") +
  xlim(c(0, 500)) +
  ylim(c(0, 0.010)) +
  labs(title = "Gene Level Read Depth (featureCounts)")+
  geom_density(col = "red") +
  facet_grid(celltype~sampleID)
p2 <- ggplot(master_T_fc, aes(x = NumReads)) +
  geom_histogram(aes(y = ..density..),
                 fill = "white", colour = "black") +
  xlim(c(0, 500)) +
  ylim(c(0, 0.010)) +
  geom_density(col = "red") +
  facet_grid(celltype~sampleID)

grid.arrange(p1, p2, nrow = 2)
ggsave("B_T_hist_fc_tx.pdf", width = 7, height = 5, dpi = 100, units = "in", device='pdf')

fitERR431592 <- fitdistr(fc_cnts$ERR431592_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431592_1.fastq.gz.subread.BAM, freq = FALSE, breaks=500,
     xlim = c(0, quantile(fc_cnts$ERR431592_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431592 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431592$estimate), from = 0, col = "red", add = TRUE)

fitERR431605 <- fitdistr(fc_cnts$ERR431605_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431605_1.fastq.gz.subread.BAM, freq = FALSE, breaks=700,
     xlim = c(0, quantile(fc_cnts$ERR431605_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431605 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431605$estimate), from = 0, col = "red", add = TRUE)

fitERR431612 <- fitdistr(fc_cnts$ERR431612_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431612_1.fastq.gz.subread.BAM, freq = FALSE, breaks=700,
     xlim = c(0, quantile(fc_cnts$ERR431612_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431612 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431612$estimate), from = 0, col = "red", add = TRUE)

fitERR431623 <- fitdistr(fc_cnts$ERR431623_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431623_1.fastq.gz.subread.BAM, freq = FALSE, breaks=700,
     xlim = c(0, quantile(fc_cnts$ERR431623_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431623 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431623$estimate), from = 0, col = "red", add = TRUE)

fitERR431627 <- fitdistr(fc_cnts$ERR431627_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431627_1.fastq.gz.subread.BAM, freq = FALSE, breaks=700,
     xlim = c(0, quantile(fc_cnts$ERR431627_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431627 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431627$estimate), from = 0, col = "red", add = TRUE)

fitERR431572 <- fitdistr(fc_cnts$ERR431572_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431572_1.fastq.gz.subread.BAM, freq = FALSE, breaks=400,
     xlim = c(0, quantile(fc_cnts$ERR431572_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.0055), main="ERR431572 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431572$estimate), from = 0, col = "red", add = TRUE)

fitERR431586 <- fitdistr(fc_cnts$ERR431586_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431586_1.fastq.gz.subread.BAM, freq = FALSE, breaks=400,
     xlim = c(0, quantile(fc_cnts$ERR431586_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431586 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431586$estimate), from = 0, col = "red", add = TRUE)

fitERR431611 <- fitdistr(fc_cnts$ERR431611_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431611_1.fastq.gz.subread.BAM, freq = FALSE, breaks=400,
     xlim = c(0, quantile(fc_cnts$ERR431611_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431611 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431611$estimate), from = 0, col = "red", add = TRUE)

fitERR431619 <- fitdistr(fc_cnts$ERR431619_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431619_1.fastq.gz.subread.BAM, freq = FALSE, breaks=400,
     xlim = c(0, quantile(fc_cnts$ERR431619_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431619 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431619$estimate), from = 0, col = "red", add = TRUE)

fitERR431624 <- fitdistr(fc_cnts$ERR431624_1.fastq.gz.subread.BAM, "exponential") 
hist(fc_cnts$ERR431624_1.fastq.gz.subread.BAM, freq = FALSE, breaks=400,
     xlim = c(0, quantile(fc_cnts$ERR431624_1.fastq.gz.subread.BAM, 0.99)), 
     ylim = c(0, 0.002), main="ERR431624 Read Depth per Gene", xlab="Read Depth per Gene")
curve(dexp(x, rate = fitERR431624$estimate), from = 0, col = "red", add = TRUE)

mean(fitERR431624$estimate, fitERR431619$estimate, fitERR431611$estimate, fitERR431586$estimate, fitERR431572$estimate,
     fitERR431627$estimate, fitERR431623$estimate, fitERR431612$estimate, fitERR431605$estimate, fitERR431605$estimate)
#rate = 0.005958647 for gene level read depth

#run DESeq2 on gene level counts to get log2FC between B and T cells
fc_cnts$ensgeneID <- rownames(fc_cnts)
fc_cnts <- fc_cnts[,c(11,1:10)]
rownames(fc_cnts) <- seq(1,nrow(fc_cnts))
dds <- DESeqDataSetFromMatrix(countData=fc_cnts, 
                              colData=samples, 
                              design=~condition, tidy = TRUE)
dds <- DESeq(dds)
DESEq2_res <- as.data.frame(results(dds))
head(results(dds, tidy=TRUE))
head(DESEq2_res)
hist(DESEq2_res$log2FoldChange)
log2FCDESeq2 <- na.omit(DESEq2_res$log2FoldChange)
fit <- fitdistr(log2FCDESeq2, "normal")
para <- fit$estimate

ggplot(data = DESEq2_res) +
  geom_histogram(mapping = aes(x= log2FoldChange, y = ..density..), col="white") +
  labs(x="log 2 Fold Change", title="DESeq2 Gene Level log 2 Fold Change") +
  stat_function(fun = dnorm, 
                args = list(mean = fit$estimate[1], sd = fit$estimate[2], log = F), 
                color="red", lwd=1)
ggsave("DESeq2_log2FC_hist.pdf", width = 6, height = 5, dpi = 100, units = "in", device='pdf')

FC_norm_sample <- rnorm(50, mean=0, sd=1.500618117) 
bins <- c(-100, -5, -2, -1, -0.5, 0.5, 1, 2, 5, 100)
lfc <- c(-10, -5, -2, -1, 0, 1, 2, 5, 10)
xcut <- cut(FC_norm_sample, bins, lfc)
xcut <- as.numeric(as.character(xcut))
hist(xcut, breaks = 20, xlim=c(-6, 6))


#using rate = 0.005958647 for read depth distribution
RD_exp_sample <- rexp(50, rate=0.005958647)
bins <- c(0,  0.5, 1, 2, 5, 10, 50, 100, 500, 1000)
rd <- c(0, 0.5, 1, 2, 5, 10, 50, 100, 500)
RDcut <- cut(RD_exp_sample, bins, rd)
RDcut <- as.numeric(as.character(RDcut))
hist(RDcut, breaks=20)

#salmon read depth
#get transcripts fasta
#rsem-prepare-reference --gtf gencode.v34.annotation.gtf --star -p 36 GRCh38.p13.genome.fa RSEM_ref
#running salmon to get transcript counts
#salmon index -t ../ref/RSEM_ref/ -i salmon_tx_idx
#cat B_naive/B_naive_samples.txt | xargs -n 1 -P 5 -I {} salmon quant -i salmon_tx_idx -l A -1 B_naive/{}_1.fastq.gz -2 B_naive/{}_2.fastq.gz -o {}_tx_quant
#cat CD8_naive/CD8_naive_samples.txt | xargs -n 1 -P 5 -I {} salmon quant -i salmon_tx_idx -l A -1 CD8_naive/{}_1.fastq.gz -2 CD8_naive/{}_2.fastq.gz -o {}_tx_quant

for(i in samples$sample){
  assign(paste0(i, "_salmon_cnt"), read.table(paste0(i,"_tx_quant/quant.sf"), header=TRUE, sep="\t")) 
}

#Bcells
ERR431572_salmon_cnt$celltype <- "B"
ERR431586_salmon_cnt$celltype <- "B"
ERR431611_salmon_cnt$celltype <- "B"
ERR431619_salmon_cnt$celltype <- "B"
ERR431624_salmon_cnt$celltype <- "B"

ERR431572_salmon_cnt$sampleID <- "ERR431572"
ERR431586_salmon_cnt$sampleID <- "ERR431586"
ERR431611_salmon_cnt$sampleID <- "ERR431611"
ERR431619_salmon_cnt$sampleID <- "ERR431619"
ERR431624_salmon_cnt$sampleID <- "ERR431624"

#TCells
ERR431592_salmon_cnt$celltype <- "T"
ERR431605_salmon_cnt$celltype <- "T"
ERR431612_salmon_cnt$celltype <- "T"
ERR431623_salmon_cnt$celltype <- "T"
ERR431627_salmon_cnt$celltype <- "T"

ERR431592_salmon_cnt$sampleID <- "ERR431592"
ERR431605_salmon_cnt$sampleID <- "ERR431605"
ERR431612_salmon_cnt$sampleID <- "ERR431612"
ERR431623_salmon_cnt$sampleID <- "ERR431623"
ERR431627_salmon_cnt$sampleID <- "ERR431627"

master_salmon_cnt <- rbind(ERR431572_salmon_cnt, ERR431586_salmon_cnt, ERR431611_salmon_cnt, ERR431619_salmon_cnt, ERR431624_salmon_cnt,
                           ERR431592_salmon_cnt, ERR431605_salmon_cnt, ERR431612_salmon_cnt, ERR431623_salmon_cnt, ERR431627_salmon_cnt)

master_salmon_cnt_B <- rbind(ERR431572_salmon_cnt, ERR431586_salmon_cnt, ERR431611_salmon_cnt, ERR431619_salmon_cnt, ERR431624_salmon_cnt)
master_salmon_cnt_T <- rbind(ERR431592_salmon_cnt, ERR431605_salmon_cnt, ERR431612_salmon_cnt, ERR431623_salmon_cnt, ERR431627_salmon_cnt)

library(gridExtra)

p1 <- ggplot(master_salmon_cnt_B, aes(x = NumReads)) +
  geom_histogram(aes(y = ..density..),
                 fill = "white", colour = "black") +
  xlim(c(0, 500)) +
  ylim(c(0, 0.010)) +
  labs(title = "Transcript Level Read Depth (Salmon)")+
  geom_density(col = "red") +
  facet_grid(celltype~sampleID)
p2 <- ggplot(master_salmon_cnt_T, aes(x = NumReads)) +
  geom_histogram(aes(y = ..density..),
                 fill = "white", colour = "black") +
  xlim(c(0, 500)) +
  ylim(c(0, 0.010)) +
  geom_density(col = "red") +
  facet_grid(celltype~sampleID)

grid.arrange(p1, p2, nrow = 2)
ggsave("B_T_hist_salmon_tx.pdf", width = 7, height = 5, dpi = 100, units = "in", device='pdf')

fitERR431592_salmon <- fitdistr(ERR431592_salmon_cnt$NumReads, "exponential") 
hist(ERR431592_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431592_salmon_cnt$NumReads, 0.99)), 
     main="ERR431592 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431592_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431605_salmon <- fitdistr(ERR431605_salmon_cnt$NumReads, "exponential") 
hist(ERR431605_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431605_salmon_cnt$NumReads, 0.99)), 
     main="ERR431605 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431605_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431612_salmon <- fitdistr(ERR431612_salmon_cnt$NumReads, "exponential") 
hist(ERR431612_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431612_salmon_cnt$NumReads, 0.99)), 
     main="ERR431612 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431612_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431623_salmon <- fitdistr(ERR431623_salmon_cnt$NumReads, "exponential") 
hist(ERR431623_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431623_salmon_cnt$NumReads, 0.99)), 
     main="ERR431623 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431623_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431627_salmon <- fitdistr(ERR431627_salmon_cnt$NumReads, "exponential") 
hist(ERR431627_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431627_salmon_cnt$NumReads, 0.99)), 
     main="ERR431627 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431627_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431572_salmon <- fitdistr(ERR431572_salmon_cnt$NumReads, "exponential") 
hist(ERR431572_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431572_salmon_cnt$NumReads, 0.99)), 
    main="ERR431572 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431572_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431586_salmon <- fitdistr(ERR431586_salmon_cnt$NumReads, "exponential") 
hist(ERR431586_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431586_salmon_cnt$NumReads, 0.99)), 
     main="ERR431586 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431586_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431611_salmon <- fitdistr(ERR431611_salmon_cnt$NumReads, "exponential") 
hist(ERR431611_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431611_salmon_cnt$NumReads, 0.99)), 
     main="ERR431611 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431611_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431619_salmon <- fitdistr(ERR431619_salmon_cnt$NumReads, "exponential") 
hist(ERR431619_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431619_salmon_cnt$NumReads, 0.99)), 
     main="ERR431619 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431619_salmon$estimate), from = 0, col = "red", add = TRUE)

fitERR431624_salmon <- fitdistr(ERR431624_salmon_cnt$NumReads, "exponential") 
hist(ERR431624_salmon_cnt$NumReads, freq = FALSE, breaks=1500,
     xlim = c(0, quantile(ERR431624_salmon_cnt$NumReads, 0.99)), 
     main="ERR431624 Read Depth per Transcript", xlab="Read Depth per Transcript")
curve(dexp(x, rate = fitERR431624_salmon$estimate), from = 0, col = "red", add = TRUE)

mean(fitERR431624_salmon$estimate, fitERR431619_salmon$estimate, fitERR431611_salmon$estimate, fitERR431586_salmon$estimate, fitERR431572_salmon$estimate,
     fitERR431627_salmon$estimate, fitERR431623_salmon$estimate, fitERR431612_salmon$estimate, fitERR431605_salmon$estimate, fitERR431592_salmon$estimate)
# rate = 0.02433752 for salmon

#using rate = 0.02433752 for read depth distribution
RD_exp_sample <- rexp(50, rate=0.02433752)
bins <- c(0,  0.5, 1, 2, 5, 10, 50, 100, 500 )
rd <- c(0, 0.5, 1, 2, 5, 10, 50, 100)
RDcut_salmon <- cut(RD_exp_sample, bins, rd)
RDcut_salmon <- as.numeric(as.character(RDcut_salmon))
hist(RDcut_salmon, breaks=20)