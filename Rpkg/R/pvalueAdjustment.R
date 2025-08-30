###########################################################
# unexported functions  copied from DESeq2
###########################################################

#' pvalueAdjustment function 
#' @export
pvalueAdjustment <- function(res, independentFiltering, filter,
                             theta, alpha, pAdjustMethod) {
  # perform independent filtering
  if (independentFiltering) {
    if (missing(filter)) {
      filter <- res$baseMean
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
    numRej  <- MatrixGenerics::colSums(filtPadj < alpha, na.rm = TRUE)
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

    #metadata(res)[["filterThreshold"]] <- filterThreshold
    #metadata(res)[["filterTheta"]] <- filterTheta
    #metadata(res)[["filterNumRej"]] <- filterNumRej
    #metadata(res)[["lo.fit"]] <- lo.fit
    #metadata(res)[["alpha"]] <- alpha
  } else {
    # regular p-value adjustment
    # does not include those rows which were removed
    # by maximum Cook's distance
    padj <- p.adjust(res$pvalue,method=pAdjustMethod)  
  }
  res$padj <- padj
  
  # add metadata to padj column
  #mcols(res)$type[names(res)=="padj"] <- "results"
  #mcols(res)$description[names(res)=="padj"] <- paste(pAdjustMethod,"adjusted p-values")
  
  res
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
