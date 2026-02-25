# Function to sum counts across columns for matching gene:exon keys
# Function to compute column-wise sum for a gene's exons
#' Sum exon counts for a comma-separated list of exon fragments
#' @param count_col Character string. A comma-separated list of exon fragment identifiers
#'   (e.g., \code{"E001,E003"}). An empty, \code{NA}, or \code{"NA"} value returns a row of
#'   \code{NA_real_}.
#' @param counts_df A data frame of exon counts with sample count columns followed by
#'   \code{gene} and \code{exon} columns.
#' @export
#' @examples
#' counts_df <- data.frame(
#'   s1 = c(10, 20, 5, 0), s2 = c(8, 15, 3, 1),
#'   gene = "GENE1", exon = c("001", "002", "003", "004"),
#'   stringsAsFactors = FALSE
#' )
#' grase::sum_exon_counts("E001,E003", counts_df)
sum_exon_counts <- function(count_col, counts_df) {
  if (is.na(count_col) || count_col == "" || count_col == "NA") return(rep(NA_real_, (ncol(counts_df)-2)))

  # Split the exon list
  exon_list <- unlist(strsplit(count_col, ","))
  exon_list <- gsub("^E", "", exon_list)

  rownames(counts_df) = counts_df$exon
  matched_rows <- counts_df[exon_list,1:(ncol(counts_df)-2)]
  # Sum across rows for each column (column-wise sum)
  return_row = matrix(data=NA, ncol=ncol(matched_rows), dimnames=list(c(), colnames(matched_rows)))
  if (nrow(matched_rows) > 0) {
    return_row[1,] = colSums(matched_rows, na.rm = TRUE)
  } 
  return (return_row)
}


#' Write exon counts in long format for downstream statistical testing
#' @param ref_counts_df A data frame of reference exonic part read counts with sample columns
#'   and an \code{event} column.
#' @param diff_counts_df A data frame of differential (alternative) exonic part read counts with
#'   the same structure as \code{ref_counts_df}.
#' @param events Character vector of event identifiers to include from the count data frames.
#' @param gene_name Character string. Gene identifier used to construct the output file name.
#' @param sampleinfo A named character vector mapping sample names to group labels.
#' @param outfilesuffix Character string. Suffix appended to output file names.
#' @param outdir Character string. Path to output directory.
#' @export
#' @examples
#' \dontrun{
#' sampleinfo <- c(sampleA = "group1", sampleB = "group1",
#'                 sampleC = "group2", sampleD = "group2")
#' ref_counts_df <- data.frame(
#'   sampleA = c(100, 80), sampleB = c(90, 70),
#'   sampleC = c(50, 40), sampleD = c(55, 45),
#'   gene = "GENE1", event = c("1", "2"),
#'   source = "3", sink = "7",
#'   ref_ex_part = c("E001", "E002"), setdiff = c("E003", "E004")
#' )
#' grase::write_exoncnt_long(ref_counts_df, ref_counts_df,
#'   events = c("1", "2"), gene_name = "GENE1",
#'   sampleinfo = sampleinfo, outfilesuffix = "bipartition", outdir = tempdir()
#' )
#' }
write_exoncnt_long <- function(ref_counts_df, diff_counts_df, events, gene_name, sampleinfo, outfilesuffix, outdir) {

  ref_long <- ref_counts_df[ref_counts_df$event %in% events,] %>% 
      pivot_longer(
      cols = all_of(names(sampleinfo)),
      names_to = "sample",
      values_to = "ref")
  diff_long <- diff_counts_df[diff_counts_df$event %in% events,] %>% 
      pivot_longer(
      cols = all_of(names(sampleinfo)),
      names_to = "sample",
      values_to = "diff")
  
  exoncnts <- cbind.data.frame(ref_long, diff = diff_long$diff)  
  exoncnts <- exoncnts %>%
      mutate(groups = sampleinfo[sample])
  exoncnts$ref <- as.numeric(exoncnts$ref)
  exoncnts$diff <- as.numeric(exoncnts$diff)
  exoncnts$n <- exoncnts$ref + exoncnts$diff
#  exoncnts <- exoncnts[!is.na(exoncnts$diff),]
#  exoncnts <- exoncnts[exoncnts$n > 10,]
  if (nrow(exoncnts) > 1) {
    write.table(exoncnts, file=paste0(outdir,'/', gene_name, '.', outfilesuffix, '.exoncnt.txt'), quote=FALSE, sep="\t")
  }

}

#' Count exonic part reads for bipartition events
#' @param bipartition_file Character string. Path to a bipartition results file (tab-separated,
#'   as written by \code{bipartition_paths}).
#' @param countmat A data frame or matrix of exon fragment read counts with sample columns
#'   followed by \code{gene} and \code{exon} columns.
#' @param sampleinfo A named character vector mapping sample names to group labels.
#' @param outdir Character string. Path to output directory.
#' @param outfilesuffix Character string. Suffix appended to output file names.
#' @export
#' @examples
#' \dontrun{
#' countmat <- read.table("GENE1.dexseq.counts.txt", header = TRUE)
#' sampleinfo <- c(sampleA = "group1", sampleB = "group1",
#'                 sampleC = "group2", sampleD = "group2")
#' grase::count_bipartitions(
#'   bipartition_file = "GENE1.bipartition.txt",
#'   countmat = countmat, sampleinfo = sampleinfo,
#'   outdir = tempdir(), outfilesuffix = "bipartition"
#' )
#' }
count_bipartitions <- function(bipartition_file, countmat, sampleinfo, outdir, outfilesuffix) {
  bipartitions <- as.data.frame(read_tsv(bipartition_file, col_types = cols(.default = "c")))  # Reads columns as characters
  bipartitions$event = rownames(bipartitions)
  bipartitions <- bipartitions[!(is.na(bipartitions$setdiff1) & is.na(bipartitions$setdiff2)),]
  if (length(bipartitions) == 0) { return (0) }
  gene <- bipartitions$gene
  n_samples = ncol(countmat)-2

  bipartitions$ref_ex_part[is.na(bipartitions$ref_ex_part)] <- 'NA'
  bipartitions$setdiff1[is.na(bipartitions$setdiff1)] <- 'NA'
  bipartitions$setdiff2[is.na(bipartitions$setdiff2)] <- 'NA'
  ref_counts_dict <- t(sapply(unique(bipartitions$ref_ex_part), sum_exon_counts, counts_df=countmat))
  diff1_counts_dict <- t(sapply(unique(bipartitions$setdiff1), sum_exon_counts, counts_df=countmat))
  diff2_counts_dict <- t(sapply(unique(bipartitions$setdiff2), sum_exon_counts, counts_df=countmat))

  ref_counts_df = matrix(0, nrow(bipartitions), n_samples)
  diff1_counts_df = matrix(0, nrow(bipartitions), n_samples)
  diff2_counts_df = matrix(0, nrow(bipartitions), n_samples)
  if (sum(!is.na(bipartitions$ref_ex_part)) > 0) ref_counts_df = as.data.frame(matrix(ref_counts_dict[bipartitions$ref_ex_part,], nrow = nrow(bipartitions), ncol=n_samples))
  if (sum(!is.na(bipartitions$setdiff1)) > 0) diff1_counts_df = as.data.frame(matrix(diff1_counts_dict[bipartitions$setdiff1,], nrow = nrow(bipartitions), ncol=n_samples))
  if (sum(!is.na(bipartitions$setdiff2)) > 0) diff2_counts_df = as.data.frame(matrix(diff2_counts_dict[bipartitions$setdiff2,], nrow = nrow(bipartitions), ncol=n_samples))

  rownames(ref_counts_df) = rownames(bipartitions) 
  rownames(diff1_counts_df) = rownames(bipartitions) 
  rownames(diff2_counts_df) = rownames(bipartitions) 
 
  colnames(ref_counts_df) = colnames(countmat[,1:n_samples]) 
  colnames(diff1_counts_df) = colnames(countmat[,1:n_samples]) 
  colnames(diff2_counts_df) = colnames(countmat[,1:n_samples]) 

  bipartitions$ref_mean = rowMeans(ref_counts_df)
  bipartitions$diff1_mean = rowMeans(diff1_counts_df)
  bipartitions$diff2_mean = rowMeans(diff2_counts_df)

  ref_counts_df$event = bipartitions$event
  diff1_counts_df$event = bipartitions$event
  diff2_counts_df$event = bipartitions$event


  bipartitions$diff_mean = do.call(pmax, c(bipartitions[, c("diff1_mean", "diff2_mean")], na.rm = TRUE))
  bipartitions <- bipartitions %>%
    mutate(
      which   = case_when(
        is.na(diff1_mean) & is.na(diff2_mean) ~ NA_character_,
        is.na(diff1_mean)                     ~ "diff2",
        is.na(diff2_mean)                     ~ "diff1",
        diff1_mean > diff2_mean               ~ "diff1",
        TRUE                                  ~ "diff2"
      ),
      setdiff = if_else(which == "diff1", setdiff1, setdiff2)
    ) %>%
    filter(!is.na(setdiff))

  ref_counts_df <- inner_join(bipartitions[,c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], ref_counts_df, by = "event")

  diff1_counts_df <- inner_join(bipartitions[bipartitions$which=='diff1',c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], diff1_counts_df, by = "event")

  diff2_counts_df <- inner_join(bipartitions[bipartitions$which=='diff2',c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')], diff2_counts_df, by = "event")

  diff_counts_df <- bind_rows(diff1_counts_df, diff2_counts_df) %>%
       arrange(as.numeric(event))


  if (nrow(bipartitions[bipartitions$diff_mean > 0,])) {
    print(paste0("writing ",gene[1]))
    write.table(bipartitions, file=paste0(outdir,'/', gene[1], '.', outfilesuffix, '.txt'), quote=FALSE, sep="\t")
    write_exoncnt_long(ref_counts_df, diff_counts_df, events=bipartitions$event, gene_name = gene[1], sampleinfo, outfilesuffix, outdir) 
  } else {
    print(paste0("zero read counts at setdiff exons. skipping ",gene[1]))
  }
 
  return (0)
}


# Like count_bipartitions() but emits BOTH setdiff1 and setdiff2 as separate
# testable events (suffixed _s1 / _s2) instead of picking the higher-count side.
# This lets the statistical test see both sides of each bipartition independently.
#' Count exonic part reads for both setdiff1 and setdiff2 of bipartition events
#' @param bipartition_file Character string. Path to a bipartition results file (tab-separated,
#'   as written by \code{bipartition_paths}).
#' @param countmat A data frame or matrix of exon fragment read counts with sample columns
#'   followed by \code{gene} and \code{exon} columns.
#' @param sampleinfo A named character vector mapping sample names to group labels.
#' @param outdir Character string. Path to output directory.
#' @param outfilesuffix Character string. Suffix appended to output file names.
#' @export
#' @examples
#' \dontrun{
#' countmat <- read.table("GENE1.dexseq.counts.txt", header = TRUE)
#' sampleinfo <- c(sampleA = "group1", sampleB = "group1",
#'                 sampleC = "group2", sampleD = "group2")
#' grase::count_bipartitions_both(
#'   bipartition_file = "GENE1.bipartition.txt",
#'   countmat = countmat, sampleinfo = sampleinfo,
#'   outdir = tempdir(), outfilesuffix = "bipartition"
#' )
#' }
count_bipartitions_both <- function(bipartition_file, countmat, sampleinfo, outdir, outfilesuffix) {
  bipartitions <- as.data.frame(read_tsv(bipartition_file, col_types = cols(.default = "c")))
  bipartitions$event = rownames(bipartitions)
  bipartitions <- bipartitions[!(is.na(bipartitions$setdiff1) & is.na(bipartitions$setdiff2)), ]
  if (nrow(bipartitions) == 0) { return(0) }
  gene <- bipartitions$gene
  n_samples <- ncol(countmat) - 2
  sample_cols <- colnames(countmat)[1:n_samples]

  bipartitions$ref_ex_part[is.na(bipartitions$ref_ex_part)] <- 'NA'
  bipartitions$setdiff1[is.na(bipartitions$setdiff1)]       <- 'NA'
  bipartitions$setdiff2[is.na(bipartitions$setdiff2)]       <- 'NA'

  ref_counts_dict   <- t(sapply(unique(bipartitions$ref_ex_part), sum_exon_counts, counts_df = countmat))
  diff1_counts_dict <- t(sapply(unique(bipartitions$setdiff1),    sum_exon_counts, counts_df = countmat))
  diff2_counts_dict <- t(sapply(unique(bipartitions$setdiff2),    sum_exon_counts, counts_df = countmat))

  ref_counts_df   <- as.data.frame(matrix(ref_counts_dict[bipartitions$ref_ex_part, ],  nrow = nrow(bipartitions), ncol = n_samples))
  diff1_counts_df <- as.data.frame(matrix(diff1_counts_dict[bipartitions$setdiff1, ],   nrow = nrow(bipartitions), ncol = n_samples))
  diff2_counts_df <- as.data.frame(matrix(diff2_counts_dict[bipartitions$setdiff2, ],   nrow = nrow(bipartitions), ncol = n_samples))

  colnames(ref_counts_df) <- colnames(diff1_counts_df) <- colnames(diff2_counts_df) <- sample_cols

  bipartitions$ref_mean   <- rowMeans(ref_counts_df,   na.rm = TRUE)
  bipartitions$diff1_mean <- rowMeans(diff1_counts_df, na.rm = TRUE)
  bipartitions$diff2_mean <- rowMeans(diff2_counts_df, na.rm = TRUE)

  # Helper: build one side's bipartition table + count dfs with suffixed event IDs
  build_side <- function(mask, setdiff_col, diff_counts, suffix) {
    bp <- bipartitions[mask, , drop = FALSE]
    if (nrow(bp) == 0) return(NULL)
    bp$event   <- paste0(bp$event, suffix)
    bp$setdiff <- bp[[setdiff_col]]
    bp$side    <- setdiff_col
    bp$diff_mean <- rowMeans(diff_counts[mask, , drop = FALSE], na.rm = TRUE)
    rc <- ref_counts_df[mask, , drop = FALSE]
    dc <- diff_counts[mask,   , drop = FALSE]
    rc$event <- bp$event
    dc$event <- bp$event
    list(bp = bp, ref = rc, diff = dc)
  }

  has_s1 <- bipartitions$setdiff1 != 'NA'
  has_s2 <- bipartitions$setdiff2 != 'NA'

  s1 <- build_side(has_s1, "setdiff1", diff1_counts_df, "_s1")
  s2 <- build_side(has_s2, "setdiff2", diff2_counts_df, "_s2")

  sides <- Filter(Negate(is.null), list(s1, s2))
  if (length(sides) == 0) return(0)

  bp_both   <- bind_rows(lapply(sides, `[[`, "bp"))
  ref_both  <- bind_rows(lapply(sides, `[[`, "ref"))
  diff_both <- bind_rows(lapply(sides, `[[`, "diff"))

  # sort by original numeric event (strip suffix first)
  ord <- order(as.numeric(sub("_s[12]$", "", bp_both$event)))
  bp_both   <- bp_both[ord, ]
  ref_both  <- ref_both[ord, ]
  diff_both <- diff_both[ord, ]

  # only write if at least one side has non-zero counts
  if (nrow(bp_both[bp_both$diff_mean > 0, ]) > 0) {
    print(paste0("writing (both sides) ", gene[1]))
    write.table(bp_both, file = paste0(outdir, '/', gene[1], '.', outfilesuffix, '.txt'),
                quote = FALSE, sep = "\t")
    # build joined ref/diff for write_exoncnt_long
    meta_cols <- c('gene', 'event', 'source', 'sink', 'ref_ex_part', 'setdiff')
    ref_joined  <- inner_join(bp_both[, meta_cols], ref_both,  by = "event")
    diff_joined <- inner_join(bp_both[, meta_cols], diff_both, by = "event")
    write_exoncnt_long(ref_joined, diff_joined,
                       events    = bp_both$event,
                       gene_name = gene[1],
                       sampleinfo, outfilesuffix, outdir)
  } else {
    print(paste0("zero read counts at setdiff exons. skipping (both sides) ", gene[1]))
  }

  return(0)
}


#' Count exonic part reads for multinomial events
#' @param multinomial_file Character string. Path to a multinomial results file (tab-separated,
#'   as written by \code{multinomial_paths}).
#' @param countmat A data frame or matrix of exon fragment read counts with sample columns
#'   followed by \code{gene} and \code{exon} columns.
#' @param sampleinfo A named character vector mapping sample names to group labels.
#' @param outdir Character string. Path to output directory.
#' @export
#' @examples
#' \dontrun{
#' countmat <- read.table("GENE1.dexseq.counts.txt", header = TRUE)
#' sampleinfo <- c(sampleA = "group1", sampleB = "group1",
#'                 sampleC = "group2", sampleD = "group2")
#' grase::count_multinomial(
#'   multinomial_file = "GENE1.multinomial.txt",
#'   countmat = countmat, sampleinfo = sampleinfo, outdir = tempdir()
#' )
#' }
count_multinomial <- function(multinomial_file, countmat, sampleinfo, outdir) {
  # 1. Load the definitions file
  multi_df <- as.data.frame(read_tsv(multinomial_file, col_types = cols(.default = "c")))
  if (nrow(multi_df) == 0) return(0)
  
  # Ensure events are uniquely identified by their row names/indices
  multi_df$event <- rownames(multi_df)
  gene_id <- multi_df$gene[1]
  n_samples <- ncol(countmat) - 2
  
  # Identify all potential alternative option columns (setdiff1, setdiff2, etc.)
  setdiff_cols <- grep("^setdiff[0-9]+$", colnames(multi_df), value = TRUE)
  
  # 2. Robust Helper to retrieve counts for a column of exon strings
  get_counts_df <- function(col_name) {
    vals <- multi_df[[col_name]]
    vals[is.na(vals)] <- 'NA' # Standardize missing values
    unique_vals <- unique(vals)
    
    # Efficiency: Calculate sums only for unique combinations found in this column
    counts_dict <- t(sapply(unique_vals, sum_exon_counts, counts_df=countmat))
    
    # Handle edge case where sapply doesn't return a matrix (single unique value)
    if (length(unique_vals) == 1) {
       counts_dict <- matrix(counts_dict, nrow=1)
       rownames(counts_dict) <- unique_vals
    }
    
    # Map counts back to every event in the gene and force into a data frame
    # This ensures mutate() and pivot_longer() downstream will not error.
    counts_out <- as.data.frame(matrix(counts_dict[vals,], nrow = nrow(multi_df), ncol=n_samples))
    colnames(counts_out) <- colnames(countmat)[1:n_samples]
    counts_out$exon_part <- vals
    return(counts_out)
  }
  
  long_dfs <- list()
  
  # 3. Process the Reference Column (The baseline for the Wald Test)
  multi_df$ref_ex_part[is.na(multi_df$ref_ex_part)] <- 'NA'
  ref_counts_df <- get_counts_df("ref_ex_part")
  
  long_dfs[["ref"]] <- ref_counts_df %>%
    mutate(event = multi_df$event) %>%
    pivot_longer(cols = all_of(names(sampleinfo)), names_to = "sample", values_to = "count") %>%
    mutate(type = "ref")
    
  # 4. Process Setdiff Columns
  for (col in setdiff_cols) {
    # FIX: Skip columns that are entirely NA (like your setdiff4 example)
    if (all(is.na(multi_df[[col]]))) {
      next
    }
    
    counts_df <- get_counts_df(col)
    
    # Add mean counts to the summary data frame
    multi_df[[paste0(col, "_mean")]] <- rowMeans(counts_df[,1:n_samples], na.rm=TRUE)
    
    # Convert to long format and store
    long_dfs[[col]] <- counts_df %>%
      mutate(event = multi_df$event) %>%
      pivot_longer(cols = all_of(names(sampleinfo)), names_to = "sample", values_to = "count") %>%
      mutate(type = col)
  }
  
  # 5. Combine everything and cleanup
  # filter(!is.na(count)) removes rows where a specific event didn't have that setdiff option
  all_long <- bind_rows(long_dfs) %>%
    filter(!is.na(count)) %>% 
    mutate(
      gene = gene_id,
      groups = sampleinfo[sample]
    ) %>%
    relocate(exon_part, .after = type) %>%
    arrange(as.numeric(event))
    
  # Save the updated definitions (with means) and the long-format counts
  write.table(multi_df, file=paste0(outdir,'/', gene_id, '.multinomial.txt'), 
              quote=FALSE, sep="\t", row.names=FALSE)
  
  if (nrow(all_long) > 0) {
    write.table(all_long, file=paste0(outdir,'/', gene_id, '.multinomial.exoncnt.txt'), 
                quote=FALSE, sep="\t", row.names=FALSE)
  }
  
  return(0)
}
