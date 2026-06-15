###########################################################
# unexported functions  copied from DESeq2
###########################################################

#' pvalueAdjustment function
#' @param res A data frame of test results containing at least a \code{pvalue} column and,
#'   when \code{independentFiltering = TRUE}, a \code{baseMean} column.
#' @param independentFiltering Logical. Whether to perform independent filtering to increase
#'   the number of discoveries.
#' @param filter Numeric vector used as the independent filter statistic. If missing and
#'   \code{independentFiltering = TRUE}, defaults to \code{res$baseMean}.
#' @param theta Numeric vector of quantile thresholds for the independent filter. If missing,
#'   a sequence from the fraction of zeros up to 0.95 is used.
#' @param alpha Numeric. The significance level for the independent filtering optimization.
#' @param pAdjustMethod Character string. The p-value adjustment method passed to
#'   \code{p.adjust} (e.g., \code{"BH"}, \code{"bonferroni"}).
#' @export
#' @examples
#' \dontrun{
#' results <- data.frame(
#'   gene = c("GENE1", "GENE1", "GENE2"),
#'   event = c("1", "2", "1"),
#'   pvalue = c(0.01, 0.04, 0.20),
#'   baseMean = c(100, 50, 10)
#' )
#' results_adj <- grase::pvalueAdjustment(results,
#'   independentFiltering = FALSE, alpha = 0.05, pAdjustMethod = "BH"
#' )
#' results_adj[, c("gene", "event", "pvalue", "padj")]
#' }
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

#' Nested BH procedure for hierarchical FDR control
#'
#' Two-stage procedure (Benjamini-Bogomolov 2014):
#' Stage 1 applies BH across genes on Simes gene-level p-values.
#' Stage 2 applies BH within each gene that passes stage 1.
#' Controls global FDR on individual events and removes gene-size detection bias.
#'
#' @param res Data frame with columns \code{gene} and \code{pvalue}.
#' @param alpha Numeric. FDR level for both stages (default 0.05).
#' @return \code{res} with two new columns:
#'   \code{padj_gene} (BH-adjusted Simes p-value, one value per gene, replicated to all rows);
#'   \code{padj} (within-gene BH adjusted p-value for events in selected genes, NA otherwise).
#' @export
#' @examples
#' \dontrun{
#' res <- data.frame(
#'   gene   = c("G1","G1","G1","G2","G2"),
#'   pvalue = c(0.001, 0.01, 0.8, 0.4, 0.5)
#' )
#' nested_bh(res, alpha = 0.05)
#' }
nested_bh <- function(res, alpha = 0.05) {
  genes <- unique(res$gene)

  # Stage 1: Simes p-value per gene
  pv_by_gene <- split(res$pvalue, res$gene)
  simes_p <- vapply(genes, function(g) {
    pv <- sort(pv_by_gene[[g]][!is.na(pv_by_gene[[g]])])
    K  <- length(pv)
    if (K == 0L) return(NA_real_)
    min(K * pv / seq_len(K))
  }, numeric(1L))

  padj_gene_vec          <- p.adjust(simes_p, method = "BH")
  names(padj_gene_vec)   <- genes
  res$padj_gene          <- padj_gene_vec[res$gene]

  # Stage 2: within-gene BH only for genes passing stage 1
  res$padj <- NA_real_
  selected <- names(padj_gene_vec)[!is.na(padj_gene_vec) & padj_gene_vec < alpha]
  for (g in selected) {
    idx          <- which(res$gene == g)
    res$padj[idx] <- p.adjust(res$pvalue[idx], method = "BH")
  }

  res
}

#' Apply p-value adjustment with optional independent filtering
#'
#' Dispatches on \code{method}:
#' \itemize{
#'   \item \code{"nested_BH"}: optionally masks low-count events via independent
#'     filtering, then runs the two-stage nested BH procedure (\code{nested_bh}).
#'     Adds \code{padj_gene} (gene-level) and \code{padj} (event-level) columns.
#'   \item any \code{p.adjust} method (e.g. \code{"BH"}, \code{"bonferroni"}):
#'     delegates directly to \code{pvalueAdjustment}, which handles independent
#'     filtering internally and adds a single \code{padj} column.
#' }
#'
#' @param res Data frame with columns \code{gene}, \code{pvalue}, and (when
#'   \code{independentFiltering = TRUE}) \code{baseMean}.
#' @param independentFiltering Logical. Apply independent filtering before
#'   adjustment.
#' @param alpha Numeric. FDR level (default 0.05).
#' @param method Character. Adjustment method: \code{"nested_BH"} (default) or
#'   any method accepted by \code{p.adjust} (e.g. \code{"BH"}).
#' @return \code{res} with \code{padj} added (and \code{padj_gene} for
#'   \code{"nested_BH"}).
#' @export
adjust_pvalues <- function(res, independentFiltering = FALSE, alpha = 0.05,
                           method = "nested_BH") {
  if (method == "nested_BH") {
    if (independentFiltering && "baseMean" %in% names(res)) {
      filt       <- pvalueAdjustment(res, independentFiltering = TRUE,
                                     alpha = alpha, pAdjustMethod = "BH")
      res$pvalue <- ifelse(is.na(filt$padj), NA_real_, res$pvalue)
    }
    nested_bh(res, alpha = alpha)
  } else {
    pvalueAdjustment(res, independentFiltering = independentFiltering,
                     alpha = alpha, pAdjustMethod = method)
  }
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
