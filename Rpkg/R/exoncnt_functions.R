# Function to sum counts across columns for matching gene:exon keys
# Function to compute column-wise sum for a gene's exons
#' @export
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


#' @export
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

#' @export
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


#' @export
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
