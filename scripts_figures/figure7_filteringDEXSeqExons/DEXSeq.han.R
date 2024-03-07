
.libPaths("/mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2")

library(DEXSeq)
library(tidyverse)


pvalueAdjustment <- function(res, independentFiltering, filter,
                             theta, alpha, pAdjustMethod) {
  # perform independent filtering
  if (independentFiltering) {
    if (missing(filter)) {
      filter <- res$exonBaseMean
    }
    if (missing(theta)) {
      lowerQuantile <- mean(filter == 0)
      if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
      theta <- seq(lowerQuantile, upperQuantile, length=50)
    }

    # do filtering using genefilter
    stopifnot(length(theta) > 1)
    stopifnot(length(filter) == nrow(res))
    filtPadj <- filtered_p(filter=filter, test=res$pvalue,
                           theta=theta, method=pAdjustMethod)
    numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
    # prevent over-aggressive filtering when all genes are null,
    # by requiring the max number of rejections is above a fitted curve.
    # If the max number of rejection is not greater than 10, then don't
    # perform independent filtering at all.
    lo.fit <- lowess(numRej ~ theta, f=1/5)
    if (max(numRej) <= 10) {
      j <- 1
    } else {

      residual <- if (all(numRej==0)) {
        0
      } else {
        numRej[numRej > 0] - lo.fit$y[numRej > 0]
      }

      # this usually works: find the threshold at which num rejections
      # surpasses the root mean squared error around the fitted curve.
      # it may not work if there is a sharp uptick in the curve at
      # the end of the grid, and there is very little variation.
      maxFit <- max(lo.fit$y)
      rmse <- sqrt(mean(residual^2))
      thresh <- maxFit - rmse

      j <- if (any(numRej > thresh)) {
             # backup case: if low variation and uptick at end,
             # pick the first point at which num rejections reaches
             # 90% of the fitted curve, or 80% of the fitted curve
             which(numRej > thresh)[1]
           } else if (any(numRej > .9 * maxFit)) {
             which(numRej > .9 * maxFit)[1]
           } else if (any(numRej > .8 * maxFit)) {
             which(numRej > .8 * maxFit)[1]
           } else {
             1
           }
    }
    # j <- which.max(numRej) # old method, not stable
    padj <- filtPadj[, j, drop=TRUE]
    cutoffs <- quantile(filter, theta)
    filterThreshold <- cutoffs[j]
    filterNumRej <- data.frame(theta=theta, numRej=numRej)
    filterTheta <- theta[j]

    metadata(res)[["filterThreshold"]] <- filterThreshold
    metadata(res)[["filterTheta"]] <- filterTheta
    metadata(res)[["filterNumRej"]] <- filterNumRej
    metadata(res)[["lo.fit"]] <- lo.fit
    metadata(res)[["alpha"]] <- alpha
  } else {
    # regular p-value adjustment
    # does not include those rows which were removed
    # by maximum Cook's distance
    padj <- p.adjust(res$pvalue,method=pAdjustMethod)
  }
  #res$padj <- padj
  res['padj'] = padj

  # add metadata to padj column
  mcols(res)$type[names(res)=="padj"] <- "results"
  mcols(res)$description[names(res)=="padj"] <- paste(pAdjustMethod,"adjusted p-values")

  return(res)
}

# function copied from `genefilter` package to avoid gfortran requirement
filtered_p <- function( filter, test, theta, data, method = "none" ) {
  if ( is.function( filter ) )
    U1 <- filter( data )
  else
    U1 <- filter
  cutoffs <- quantile( U1, theta )
  result <- matrix( NA_real_, length( U1 ), length( cutoffs ) )
  colnames( result ) <- names( cutoffs )
  for ( i in 1:length( cutoffs ) ) {
    use <- U1 >= cutoffs[i]
    if( any( use ) ) {
      if( is.function( test ) )
        U2 <- test( data[use,] )
      else
        U2 <- test[use]
      result[use,i] <- p.adjust( U2, method )
    }
  }
  return( result )
}


#
##filtered unique 112402 DEXSeq exons with mapped AS events that don't have any NA's
#dexseq_with_AS_events <- read.table("/mnt/storage/jaquino/GrASE_results_bonnal/before_after_filtering_DEXSeq/dexseq_res_all_filt_after.txt", header=TRUE, sep="\t")
##rMATS_tested exons with NA's and unique 138997 exon ID 
rMATS_TestedExons <- read.table("/mnt/storage/jaquino/GrASE_results_bonnal/ExonParts/rMATS_TestedExons.txt", header=TRUE, sep="\t")
#rMATS_TestedExons$featureID <- gsub("E", "", rMATS_TestedExons$featureID)
rMATS_TestedExons <- rMATS_TestedExons %>% unite(col = "DEXSeqID", c(groupID, featureID), sep=":")
#
##non aggregated gff
#flattenedFile = list.files("/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis",
#                           pattern="gencode.v34.annotation.dexseq.gff",
#                           full.names=TRUE) #read in dexseq annotation file
#sampleTable = data.frame(
#  row.names = c( "Bnaive1", "Bnaive2", "Bnaive3", "Bnaive4", "Bnaive5",
#                "CD8naive1", "CD8naive2", "CD8naive3", "CD8naive4", "CD8naive5"),
#  condition = c( "Bnaive", "Bnaive", "Bnaive", "Bnaive", "Bnaive",
#                 "CD8naive", "CD8naive", "CD8naive", "CD8naive", "CD8naive"))
#origcountFiles = list.files("/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/original_dexseq_cnts/", pattern=".cnts.txt$", full.names=TRUE)
#countFiles = list.files("/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis", pattern=".cnts.txt$", full.names=TRUE)
#
#listOfFiles <- lapply(origcountFiles, function(x) read.table(x, header=FALSE, sep="\t"))
##listOfFiles <- lapply(listOfFiles, function(z) z[z$V1 %in% rMATS_TestedExons$DEXSeqID,]) #keep dexseq exons that mapped to AS events
##listOfFiles <- lapply(listOfFiles, function(z) z[z$V1 %in% dexseq_with_AS_events$DEXSeqID,])
#
##save new counts files
##sapply(1:length(listOfFiles), function(x) 
##  write.table(listOfFiles[x], 
##            countFiles[x], 
##            row.names = FALSE, quote=FALSE, col.names = FALSE))
#
#dxd = DEXSeqDataSetFromHTSeq(
#  origcountFiles,
#  sampleData = sampleTable,
#  design = ~ sample + exon + condition:exon,
#  flattenedfile = flattenedFile )
#
#this_other <- as.data.frame(counts(dxd))
#write.table(this_other, "this_others_counts.orig.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#
#BPPARAM = MulticoreParam(16)
#dxd = estimateSizeFactors( dxd )
#dxd = estimateDispersions( dxd, BPPARAM=BPPARAM )
##dxd = estimateDispersions( dxd, fitType='mean' )        ((keep commented as it was originally commented out))
#dxd = testForDEU( dxd, BPPARAM=BPPARAM )
#dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
#saveRDS(dxd, file="dxd.orig.rds")

dxd <- readRDS("dxd.orig.rds")
dxr1 = DEXSeqResults( dxd )
origtheta <- metadata(dxr1)$filterNumRej[,"theta"]
dxr1df <- as.data.frame(dxr1)
rownames(dxr1df) <- NULL
dxr1df$transcripts <- vapply(dxr1df$transcripts, paste, collapse = ", ", character(1L))
dxr1df <- unite(dxr1df, col='countData.case_Bnaive', c('countData.Bnaive1','countData.Bnaive2','countData.Bnaive3', 'countData.Bnaive4', 'countData.Bnaive5'), sep=',')
dxr1df <- unite(dxr1df, col='countData.case_CD8naive', c('countData.CD8naive1','countData.CD8naive2', 'countData.CD8naive3', 'countData.CD8naive4', 'countData.CD8naive5'), sep=',')
write.table(dxr1df,"/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/bonnal_dexseq_results.orig.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


dxr1df$DEXSeqID = apply(dxr1df[,c("groupID", "featureID")], 1, paste, collapse=":")
keep <- dxr1df$DEXSeqID %in% rMATS_TestedExons$DEXSeqID
print(paste("DEXseq with rMATS cnt: ",sum(keep)))
dxr1df_filtered <- dxr1df[keep,]
write.table(dxr1df_filtered,"/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/bonnal_dexseq_results.orig.filtered.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
dxd = dxd[keep,]

# rerun after filtering
dxr2 = dxr1_filtered = dxr1[keep,]
dxr2 <- pvalueAdjustment(dxr2, independentFiltering=TRUE, filter=dxr2$exonBaseMean, theta = origtheta, alpha=0.1,  pAdjustMethod="BH")
saveRDS(dxd, file="dxd.pvalueAdjustment.rds")
dxr2df <- as.data.frame(dxr2)
rownames(dxr2df) <- NULL
dxr2df$transcripts <- vapply(dxr2df$transcripts, paste, collapse = ", ", character(1L))
dxr2df <- unite(dxr2df, col='countData.case_Bnaive', c('countData.Bnaive1','countData.Bnaive2','countData.Bnaive3', 'countData.Bnaive4', 'countData.Bnaive5'), sep=',')
dxr2df <- unite(dxr2df, col='countData.case_CD8naive', c('countData.CD8naive1','countData.CD8naive2', 'countData.CD8naive3', 'countData.CD8naive4', 'countData.CD8naive5'), sep=',')
write.table(dxr2df,"/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/bonnal_dexseq_results.pvalueAdjustment.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

print(paste("DEXseq/rMATS tested original: ", sum(!is.na(dxr1_filtered$padj))))
print(paste("DEXseq/rMATS sig (0.05) original: ", sum(dxr1_filtered$padj < 0.05, na.rm=TRUE)))
print(paste("DEXseq/rMATS sig (0.1) original: ", sum(dxr1_filtered$padj < 0.1, na.rm=TRUE)))

print(paste("DEXseq/rMATS tested after filtering: ", sum(!is.na(dxr2$padj))))
print(paste("DEXseq/rMATS sig (0.05) after filtering: ", sum(dxr2$padj < 0.05, na.rm=TRUE)))
print(paste("DEXseq/rMATS sig (0.1) after filtering: ", sum(dxr2$padj < 0.1, na.rm=TRUE)))


# rerun after filtering
dxr3 = DEXSeqResults( dxd )
saveRDS(dxd, file="dxd.newrun.rds")
dxr3df <- as.data.frame(dxr3)
rownames(dxr3df) <- NULL
dxr3df$transcripts <- vapply(dxr3df$transcripts, paste, collapse = ", ", character(1L))
dxr3df <- unite(dxr3df, col='countData.case_Bnaive', c('countData.Bnaive1','countData.Bnaive2','countData.Bnaive3', 'countData.Bnaive4', 'countData.Bnaive5'), sep=',')
dxr3df <- unite(dxr3df, col='countData.case_CD8naive', c('countData.CD8naive1','countData.CD8naive2', 'countData.CD8naive3', 'countData.CD8naive4', 'countData.CD8naive5'), sep=',')
write.table(dxr3df,"/mnt/storage/jaquino/bonnal_dexseq_counts/DEXSeq_analysis/bonnal_dexseq_results.newrun.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

print(paste("DEXseq/rMATS tested original: ", sum(!is.na(dxr1_filtered$padj))))
print(paste("DEXseq/rMATS sig (0.05) original: ", sum(dxr1_filtered$padj < 0.05, na.rm=TRUE)))
print(paste("DEXseq/rMATS sig (0.1) original: ", sum(dxr1_filtered$padj < 0.1, na.rm=TRUE)))

print(paste("DEXseq/rMATS tested after filtering: ", sum(!is.na(dxr3$padj))))
print(paste("DEXseq/rMATS sig (0.05) after filtering: ", sum(dxr3$padj < 0.05, na.rm=TRUE)))
print(paste("DEXseq/rMATS sig (0.1) after filtering: ", sum(dxr3$padj < 0.1, na.rm=TRUE)))



dexseq_successful_pvalues <- dxr1df %>%
  filter(!is.na(pvalue) & is.na(padj))
dexseq_unsuccessful_pvalues <- dxr1df %>%
  filter(is.na(pvalue) & is.na(padj))
dexseq_successful_both <- dxr1df %>%
  filter(!is.na(pvalue) & !is.na(padj))

dexseq_successful <- rbind(dexseq_successful_both, dexseq_successful_pvalues)
dexseq_successful$corrected.padj <- p.adjust(dexseq_successful$pvalue, method="BH")
#write.table(dexseq_successful,"bonnal_dexseq_corrected_padj.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
