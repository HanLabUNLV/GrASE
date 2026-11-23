#!/usr/bin/env Rscript
## satuRn on DEXSeq exonic-part counts, one two-group contrast: each part is a
## "transcript" of its gene, so testDTU tests differential exon usage. Parts are kept
## with >= 10 reads in >= 3 samples and some reads in both groups, and the fit is run
## twice: on those parts (sumExp_filteredbyCount) and on those of genes with more than
## one kept part (sumExp_filteredbyCountMultiExon). The estimate is group 1 minus
## group 2 (log odds).
## Output, in --outdir: all.<cell1>_<cell2>.<set>.txt (every fitted part),
## DEU.<cell1>_<cell2>.<set>.txt (empirical FDR < 0.05), and the dxd and fitted sumExp rds.
## usage: Rscript saturn_exon.R --gff <flattened gff> --cntdir <dir with <cell>/<sample>_counts.txt>
##          --outdir <dir> --cell1 <group 1> --cell2 <group 2> [--workers 15]
suppressPackageStartupMessages({
  library(DEXSeq); library(tidyverse); library(BiocParallel); library(optparse); library(satuRn)
})

option_list = list(
  make_option(c("--gff"), type="character", help="DEXSeq flattened gff file", metavar="character"),
  make_option(c("--cntdir"), type="character", help="count directory, one subdirectory per group", metavar="character"),
  make_option(c("--outdir"), type="character", help="output directory", metavar="character"),
  make_option(c("--cell1"), type="character", help="group 1", metavar="character"),
  make_option(c("--cell2"), type="character", help="group 2", metavar="character"),
  make_option(c("--workers"), type="integer", default=15, help="fitDTU workers [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))
required_opts <- c("gff", "cntdir", "outdir", "cell1", "cell2")
missing_opts <- required_opts[sapply(required_opts, function(x) is.null(opt[[x]]))]
if (length(missing_opts) > 0) stop("missing options: ", paste(missing_opts, collapse = ", "))
GFFfile = opt$gff
celltype1 = opt$cell1
celltype2 = opt$cell2
outdir = opt$outdir
countdir = opt$cntdir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Get count files from both celltype1 and celltype2 folders
countFiles1 <- list.files(paste0(countdir,"/",celltype1),
                           pattern = "_counts.txt$",
                           full.names = TRUE)

countFiles2 <- list.files(paste0(countdir,"/",celltype2),
                             pattern = "_counts.txt$",
                             full.names = TRUE)

# samplesIDs
samples1 <- sub("_counts\\.txt$", "", basename(countFiles1))
samples2 <- sub("_counts\\.txt$", "", basename(countFiles2))
sampleNames = c(paste0(celltype1, '_', samples1), paste0(celltype2, '_', samples2))
conditions = c(rep(celltype1,length(samples1)), rep(celltype2, length(samples2)))
 
# Create sample table
sampleData <- data.frame(
  row.names = sampleNames,
  condition = factor(conditions),
  sample_id = factor(sampleNames)  # Add sample column for design
)


# Combine count file paths in the correct sample order
countFiles <- c(countFiles1, countFiles2)
names(countFiles) = sampleNames

## DEXSeqDataSetFromHTSeq, with the count files read as tab-separated and unquoted.
## The counters write part ids as "GENE":"001"; DEXSeqDataSetFromHTSeq reads them with
## read.table's defaults (whitespace and quote characters), which splits such an id
## into two fields, so every file stops at its first summary line (_ambiguous ...).
## Otherwise as DEXSeq 1.52.0: summary lines dropped, ids GENE:E001, rows matched to
## the flattened annotation's exonic parts.
dxd_from_counts <- function(countfiles, sampleData, design, flattenedfile) {
  lf <- lapply(countfiles, function(x) {
    d <- read.table(x, header = FALSE, sep = "\t", quote = "", comment.char = "",
                    stringsAsFactors = FALSE, colClasses = c("character", "integer"))
    d$V1 <- gsub('"', "", d$V1, fixed = TRUE)
    d
  })
  if (!all(sapply(lf[-1], function(x) identical(x$V1, lf[[1]]$V1))))
    stop("Count files have differing gene ID column.")
  dcounts <- sapply(lf, `[[`, "V2")
  rownames(dcounts) <- lf[[1]]$V1
  dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ]
  rownames(dcounts) <- sub(":", ":E", rownames(dcounts))
  colnames(dcounts) <- countfiles
  splitted <- strsplit(rownames(dcounts), ":")
  exons <- sapply(splitted, "[[", 2)
  genesrle <- sapply(splitted, "[[", 1)
  aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, header = FALSE)
  colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
  aggregates$strand <- gsub("\\.", "*", aggregates$strand)
  aggregates <- aggregates[which(aggregates$class == "exonic_part"), ]
  aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
  aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
  transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
  transcripts <- strsplit(transcripts, "\\+")
  exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", aggregates$attr)
  exoninfo <- GRanges(as.character(aggregates$chr),
                      IRanges(start = aggregates$start, end = aggregates$end),
                      strand = aggregates$strand)
  names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")
  names(transcripts) <- names(exoninfo)
  if (!all(rownames(dcounts) %in% names(exoninfo)))
    stop("Count files do not correspond to the flattened annotation file")
  matching <- match(rownames(dcounts), names(exoninfo))
  stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
  DEXSeqDataSet(dcounts, sampleData, design, exons, genesrle,
                exoninfo[matching], transcripts[matching])
}

dxdfile = paste0(outdir,"/dxd.",celltype1,"_",celltype2,".rds")
if (file.exists(dxdfile)) {
  dxd <- readRDS(dxdfile)
} else {
  # Create DEXSeq dataset
  dxd <- dxd_from_counts(
    countFiles,
    sampleData = sampleData,
    design = ~ sample_id + exon + condition:exon,
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
print(dim(sumExp))

# countsum > 0 in both groups
counts_mat <- assay(sumExp, "counts")
is_celltype1 <- colData(sumExp)$condition == celltype1 
is_celltype2 <- colData(sumExp)$condition == celltype2
safe_rows <- rowSums(counts_mat[, is_celltype1]) > 0 & rowSums(counts_mat[, is_celltype2]) > 0
sumExp <- sumExp[safe_rows, ]
sumExp_filteredbyCount <- sumExp
print(dim(sumExp))

## proportion > 0.05 in at least 5 samples 
#gene_ids <- rowData(sumExp)$gene_id
#counts_mat <- assay(sumExp, "counts")
#gene_sums <- rowsum(counts_mat, gene_ids)
#exon_gene_totals <- gene_sums[match(gene_ids, rownames(gene_sums)), ]
#usage_prop <- counts_mat / exon_gene_totals
#keep_prop <- rowSums(usage_prop >= 0.05, na.rm = TRUE) >= 5
#keep_prop[is.na(keep_prop)] <- FALSE
#sumExp <- sumExp[keep_prop, ]
#sumExp3 <- sumExp
#print(dim(sumExp))

# filter out genes with only one exon
keep_groups <- names(which(table(rowData(sumExp)$gene_id) > 1))
sumExp <- sumExp[rowData(sumExp)$gene_id %in% keep_groups, ]
sumExp_filteredbyCountMultiExon <- sumExp
print(dim(sumExp))

# Identify rows where variance is practically zero
#counts_mat <- assay(sumExp, "counts")
#row_vars <- apply(counts_mat, 1, var)
#problem_rows <- which(row_vars < 0.01 | is.na(row_vars))
#if(length(problem_rows) > 0) {
#    message("Removing ", length(problem_rows), " numerically unstable rows.")
#    sumExp <- sumExp[-problem_rows, ]
#}
#print(dim(sumExp))

Explist = list(sumExp_filteredbyCount=sumExp_filteredbyCount, sumExp_filteredbyCountMultiExon = sumExp_filteredbyCountMultiExon)
for (expnum in 1:(length(Explist)) ) {
sumexpfile = paste0(outdir,"/sumexp.",celltype1,"_",celltype2,".",names(Explist)[expnum],".rds")

# logit(p)=βcondI(cond1)+βcondI(cond2) no intercept.
# One coefficient per conditiona
# Each coefficient is a group-specific log-odds

if (file.exists(sumexpfile)) {
  sumExp <- readRDS(sumexpfile)
} else {
  sumExp <- satuRn::fitDTU(object = Explist[[expnum]],
                         formula = ~ 0 + condition,
                         parallel = FALSE,
                         BPPARAM = MulticoreParam(workers=opt$workers),
                         verbose = TRUE)


  saveRDS(sumExp, sumexpfile)
}

print(paste0(names(Explist)[expnum], ":"))

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
#print(head(rowData(sumExp)$fitDTUModels)) 

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
count <- counts(dxd)
exondata <- rowData(sumExp)
ids = rownames(exondata)
fitcount <- counts(dxd)[ids,] 
group1_names = rownames(sampleData)[sampleData$condition==celltype1]
group2_names = rownames(sampleData)[sampleData$condition==celltype2] 
group1_cnts = fitcount[,as.character(colData(dxd)$sample)  %in% group1_names]
group2_cnts = fitcount[,as.character(colData(dxd)$sample)  %in% group2_names]
group1_exon_mean = rowMeans(group1_cnts[,1:length(group1_names)])
group1_total_mean = rowMeans(group1_cnts[,(1+length(group1_names)):(length(group1_names)+length(group1_names))])
group2_exon_mean = rowMeans(group2_cnts[,1:length(group2_names)])
group2_total_mean = rowMeans(group2_cnts[,(1+length(group2_names)):(length(group2_names)+length(group2_names))])
stat <- exondata[,c('exonBaseMean', 'exonBaseVar')]
models <- rowData(sumExp)[,"fitDTUModels"]
coef_list <- lapply(ids, function(id) {
  m <- models[[id]]
  if (is.null(m)) return(NULL)
  m@params$coefficients
})
names(coef_list) <-ids 
coef_list <- Filter(Negate(is.null), coef_list)
# union of all coefficient names
all_coef_names <- sort(unique(unlist(lapply(coef_list, names))))
# build matrix (genes x coefficients) Log-odds change in isoform usage relative to the reference
coef_mat <- do.call(rbind, lapply(coef_list, function(co) {
  out <- setNames(rep(NA_real_, length(all_coef_names)), all_coef_names)
  out[names(co)] <- unname(co)
  out
}))

fitresult <- exondata[,'fitDTUResult_C1']   # fit results estimate βCD-βCD4_N_STI
OR = exp(fitresult$estimates)
results <- cbind.data.frame(ids, group1_exon_mean, group1_total_mean, group2_exon_mean, group2_total_mean, stat,  coef_mat, fitresult, OR)
rownames(results) = ids
totalfile = paste0(outdir,"/all.",celltype1,"_",celltype2,".",names(Explist)[expnum],".txt")
write.table(results, totalfile, sep="\t", quote=FALSE, row.names=FALSE)
# filter DEU at FDR < 0.05
DEU <- rownames(exondata[["fitDTUResult_C1"]][which(exondata[["fitDTUResult_C1"]]$empirical_FDR < 0.05),])
DEUfile = paste0(outdir,"/DEU.",celltype1,"_",celltype2,".",names(Explist)[expnum],".txt")
write.table(results[DEU,], DEUfile, sep="\t", quote=FALSE, row.names=FALSE)

}
