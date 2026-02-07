library(DEXSeq)
library(tidyverse)
library(BiocParallel)
library(optparse)
library(satuRn)

# Set working directory
dicedir = "/mnt/data1/home/mirahan/DICE/"
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--dir1"), type="character",  
              help="count_files directory 1", metavar="character"),
  make_option(c("--dir2"), type="character", 
              help="count_files directory 2", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

# Flattened annotation file
GFFfile <- "gencode.v34.dexseq.bygene.gff"

celltype1 = 'CD4'
celltype2 = 'CD4_N_STIM'
if (! is.null(opt$dir1)) {
  celltype1 = opt$dir1
  celltype2 = opt$dir2
}

# Get count files from both B and CD8 folders
countFiles1 <- list.files(paste0(dicedir,"/count_files/",celltype1),
                           pattern = "_counts.txt.clean$",
                           full.names = TRUE)

countFiles2 <- list.files(paste0(dicedir,"/count_files/",celltype2),
                             pattern = "_counts.txt.clean$",
                             full.names = TRUE)

# samplesIDs
samples1 <- sub("_.*", "", basename(countFiles1))
samples2 <- sub("_.*", "", basename(countFiles2))
filenames <- c(samples1, samples2)
sampleNames = c(paste0(celltype1, '_', samples1), paste0(celltype2, '_', samples2))
conditions = c(rep(celltype1,length(samples1)), rep(celltype2, length(samples2)))
 
# Create sample table
sampleData <- data.frame(
  row.names = filenames,
  condition = factor(conditions),
  sample_id = factor(sampleNames)  # Add sample column for design
)


# Combine count file paths in the correct sample order
countFiles <- c(countFiles1, countFiles2)
names(countFiles) = sub("_.*", "", basename(countFiles))

dxdfile = paste0(dicedir,"/saturn_CD4_CD_N_STIM_7/dxd.",celltype1,"_",celltype2,".rds")
if (file.exists(dxdfile)) {
  dxd <- readRDS(dxdfile)
} else {
  # Create DEXSeq dataset
  dxd <- DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData = sampleData,
    design = ~ sample + exon + condition:exon,
    flattenedfile = GFFfile
  )
  saveRDS(dxd, dxdfile)
}
print(dim(dxd))

# create SummarizedExperiment
exonInfo <- rowData(dxd)
colnames(exonInfo)[1:2] <- c("isoform_id", "gene_id")
exonInfo$isoform_id <- rownames(exonInfo)
sumExp <- SummarizedExperiment(assays = list(counts=featureCounts(dxd)),
                               colData = sampleAnnotation(dxd),
                               rowData = exonInfo)

sumExp0 <- sumExp
print(dim(sumExp))
# remove exons with zero expression
#keep <- rowSums(assay(sumExp, "counts")) != 0

# at least 10 reads in at least 3 samples
keep <- rowSums(assay(sumExp, "counts") >= 10) >= 3
sumExp <- sumExp[keep, ]
sumExp1 <- sumExp
print(dim(sumExp))

# countsum > 0 in both groups
counts_mat <- assay(sumExp, "counts")
is_celltype1 <- colData(sumExp)$condition == celltype1 
is_celltype2 <- colData(sumExp)$condition == celltype2
safe_rows <- rowSums(counts_mat[, is_celltype1]) > 0 & rowSums(counts_mat[, is_celltype2]) > 0
sumExp <- sumExp[safe_rows, ]
sumExp2 <- sumExp
print(dim(sumExp))

# proportion > 0.05 in at least 5 samples 
gene_ids <- rowData(sumExp)$gene_id
counts_mat <- assay(sumExp, "counts")
gene_sums <- rowsum(counts_mat, gene_ids)
exon_gene_totals <- gene_sums[match(gene_ids, rownames(gene_sums)), ]
usage_prop <- counts_mat / exon_gene_totals
keep_prop <- rowSums(usage_prop >= 0.05, na.rm = TRUE) >= 5
keep_prop[is.na(keep_prop)] <- FALSE
sumExp <- sumExp[keep_prop, ]
sumExp3 <- sumExp
print(dim(sumExp))

# filter out genes with only one exon
keep_groups <- names(which(table(rowData(sumExp)$gene_id) > 1))
sumExp <- sumExp[rowData(sumExp)$gene_id %in% keep_groups, ]
sumExp4 <- sumExp
print(dim(sumExp))

# Identify rows where variance is practically zero
counts_mat <- assay(sumExp, "counts")
row_vars <- apply(counts_mat, 1, var)
problem_rows <- which(row_vars < 0.01 | is.na(row_vars))
if(length(problem_rows) > 0) {
    message("Removing ", length(problem_rows), " numerically unstable rows.")
    sumExp <- sumExp[-problem_rows, ]
}
sumExp5 <- sumExp
print(dim(sumExp))

Explist = list(sumExp1, sumExp2, sumExp3, sumExp4, sumExp5)
for (expnum in 1:(length(Explist)) ) {
sumexpfile = paste0(dicedir,"/saturn_CD4_CD_N_STIM_7/sumexp.",celltype1,"_",celltype2,".sumExp",expnum,".rds")

# logit(p)=βcondI(cond1)+βcondI(cond2) no intercept.
# One coefficient per conditiona
# Each coefficient is a group-specific log-odds

if (file.exists(sumexpfile)) {
  sumExp <- readRDS(sumexpfile)
} else {
  sumExp <- satuRn::fitDTU(object = Explist[[expnum]],
                         formula = ~ 0 + condition,
                         parallel = FALSE,
                         BPPARAM = MulticoreParam(workers=15),
                         verbose = TRUE)


  saveRDS(sumExp, sumexpfile)
}

print(paste0("sumExp", expnum, ":"))

# This finds the indices where the model is NOT an error
valid_models <- sapply(rowData(sumExp)$fitDTUModels, function(x) {
    # Check if the object is a valid StatModel and has coefficients
    is(x, "StatModel") && !is.null(x@params$coefficients)
})

# 3. Subset the object to ONLY include the successful fits
message("Successful fits: ", sum(valid_models))
message("Failed fits (including those 66 mu errors): ", sum(!valid_models))

# check the fit
table(sapply(rowData(sumExp)$fitDTUModels, function(x) x@type))

#  coefficient names 
print(head(rowData(sumExp)$fitDTUModels)) 

# Re-create the contrast matrix to match the formula ~ 0 + condition
# If the levels are "CD4" and "CD4_N_STIM", fitDTU likely named them 
# "conditionCD4" and "conditionCD4_N_STIM"

# Try defining the contrast names based on the design matrix used in fitDTU
group_names <- levels(colData(sumExp)$condition)
# Assuming the fit created names like "conditionCD4"
coef_names <- paste0("condition", group_names) 

L <- matrix(0, ncol = 1, nrow = length(coef_names))
rownames(L) <- coef_names
colnames(L) <- "C1"
L

# Adjust these indices based on which celltype is 1 and 2
L[paste0("condition", celltype1), 1] <- 1
L[paste0("condition", celltype2), 1] <- -1
L

sumExp <- satuRn::testDTU(object = sumExp, 
                          contrasts = L, 
                          diagplot1 = TRUE, 
                          diagplot2 = TRUE, 
                          sort = FALSE, 
                          forceEmpirical = TRUE)



# differential exons
DEU <- rownames(rowData(sumExp)[["fitDTUResult_C1"]][which(rowData(sumExp)[["fitDTUResult_C1"]]$empirical_FDR < 0.05),])
DEUcount <- counts(dxd)[DEU,] 
group1_names = rownames(sampleData)[sampleData$condition==celltype1]
group2_names = rownames(sampleData)[sampleData$condition==celltype2] 
group1_cnts = DEUcount[,as.character(colData(dxd)$sample)  %in% group1_names]
group2_cnts = DEUcount[,as.character(colData(dxd)$sample)  %in% group2_names]
group1_exon_mean = rowMeans(group1_cnts[,1:length(group1_names)])
group1_total_mean = rowMeans(group1_cnts[,(1+length(group1_names)):(length(group1_names)+length(group1_names))])
group2_exon_mean = rowMeans(group2_cnts[,1:length(group2_names)])
group2_total_mean = rowMeans(group2_cnts[,(1+length(group2_names)):(length(group2_names)+length(group2_names))])
DEUstat <- rowData(sumExp)[DEU,][,c('exonBaseMean', 'exonBaseVar')]

models <- rowData(sumExp)[DEU,"fitDTUModels"]
coef_list <- lapply(DEU, function(id) {
  m <- models[[id]]
  if (is.null(m)) return(NULL)
  m@params$coefficients
})
names(coef_list) <-DEU 
coef_list <- Filter(Negate(is.null), coef_list)
# union of all coefficient names
all_coef_names <- sort(unique(unlist(lapply(coef_list, names))))
# build matrix (genes x coefficients) Log-odds change in isoform usage relative to the reference
coef_mat <- do.call(rbind, lapply(coef_list, function(co) {
  out <- setNames(rep(NA_real_, length(all_coef_names)), all_coef_names)
  out[names(co)] <- unname(co)
  out
}))

fitresult <- rowData(sumExp)[DEU,'fitDTUResult_C1']   # fit results estimate βCD-βCD4_N_STI
OR = exp(fitresult$estimates)
DEUresults <- cbind.data.frame(DEU, group1_exon_mean, group1_total_mean, group2_exon_mean, group2_total_mean, DEUstat,  coef_mat, fitresult, OR)
DEUfile = paste0(dicedir,"/saturn_CD4_CD_N_STIM_7/DEU.",celltype1,"_",celltype2,".sumExp",expnum,".txt")
write.table(DEUresults, DEUfile, sep="\t", quote=FALSE, row.names=FALSE)

 
}
