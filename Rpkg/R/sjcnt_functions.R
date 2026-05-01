
#' Read one STAR SJ.out.tab file
#'
#' @param path Path to a STAR SJ.out.tab file.
#' @return Data frame with columns chr, start, end, strand, motif, annotated,
#'   n_uniq, n_multi, max_overhang, and key (chr:start:end).
#' @export
read_sj_file <- function(path) {
  df <- read.table(path, header=FALSE, sep="\t", stringsAsFactors=FALSE,
    col.names=c("chr","start","end","strand","motif","annotated",
                "n_uniq","n_multi","max_overhang"))
  df$key <- paste0(df$chr, ":", df$start, ":", df$end)
  df
}


#' Build a junction count matrix from a directory of STAR SJ.out.tab files
#'
#' @param sj_dir Directory containing *.SJ.out.tab files (searched recursively).
#' @param use_multi Logical. If TRUE, add n_multi to n_uniq (default FALSE).
#' @return Named list: mat (junction x sample integer matrix), samples (sample names).
#' @export
build_sj_matrix <- function(sj_dir, use_multi=FALSE) {
  files <- list.files(sj_dir, pattern="\\.SJ\\.out\\.tab$", full.names=TRUE, recursive=TRUE)
  if (length(files) == 0L) stop("no SJ.out.tab files found in: ", sj_dir)
  samples <- sub("\\.SJ\\.out\\.tab$", "", basename(files))

  sj_list <- lapply(files, read_sj_file)

  all_keys <- unique(unlist(lapply(sj_list, function(x) x$key)))
  mat <- matrix(0L, nrow=length(all_keys), ncol=length(samples),
                dimnames=list(all_keys, samples))

  for (j in seq_along(samples)) {
    df  <- sj_list[[j]]
    cnt <- if (use_multi) df$n_uniq + df$n_multi else df$n_uniq
    idx <- match(df$key, all_keys)
    mat[idx, j] <- cnt
  }
  list(mat=mat, samples=samples)
}


#' Sum SJ counts for a comma-separated string of junction identifiers
#'
#' Analogous to \code{sum_exon_counts} in exoncnt_functions.R.
#'
#' @param junction_str Comma-separated "chr:start:end" junction identifiers,
#'   or NA / "" / "NA" to return a row of NA.
#' @param sj_mat Junction count matrix (rows = junction keys, cols = samples).
#' @return Numeric vector of length ncol(sj_mat).
#' @export
sum_sj_counts <- function(junction_str, sj_mat) {
  if (is.na(junction_str) || junction_str == "" || junction_str == "NA") {
    return(rep(NA_real_, ncol(sj_mat)))
  }
  keys <- trimws(unlist(strsplit(junction_str, ",")))
  keys <- keys[nchar(keys) > 0L]
  # The graph stores to_pos as either intron_end or next_exon_start depending on
  # edge type; try both end and end-1 to handle either convention without
  # double-counting (each physical intron has exactly one entry in the SJ file).
  parts     <- strsplit(keys, ":")
  keys_alt  <- vapply(parts, function(p)
    paste0(p[1L], ":", p[2L], ":", as.integer(p[3L]) - 1L), character(1L))
  all_keys  <- unique(c(keys, keys_alt))
  idx <- match(all_keys, rownames(sj_mat))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) return(rep(0, ncol(sj_mat)))
  if (length(idx) == 1L) return(as.numeric(sj_mat[idx, ]))
  colSums(sj_mat[idx, , drop=FALSE])
}


#' Write SJ counts in wide long format for downstream statistical testing
#'
#' Output format mirrors \code{write_exoncnt_wide} and is directly compatible
#' with exontest.R.
#'
#' @param ref_counts_df   Data frame of shared intron counts with sample columns.
#' @param diff1_counts_df Data frame of distinct1 intron counts.
#' @param diff2_counts_df Data frame of distinct2 intron counts.
#' @param events          Character vector of event identifiers to include.
#' @param gene_name       Gene identifier string.
#' @param sampleinfo      Named character vector mapping sample names to groups.
#' @param outfilesuffix   Suffix for the output filename.
#' @param outdir          Output directory.
#' @export
write_sjcnt_wide <- function(ref_counts_df, diff1_counts_df, diff2_counts_df,
                             events, gene_name, sampleinfo, outfilesuffix, outdir) {
  ref_long <- tidyr::pivot_longer(ref_counts_df[ref_counts_df$event %in% events, ],
    cols=dplyr::all_of(names(sampleinfo)), names_to="sample", values_to="ref")
  diff1_long <- tidyr::pivot_longer(diff1_counts_df[diff1_counts_df$event %in% events, ],
    cols=dplyr::all_of(names(sampleinfo)), names_to="sample", values_to="diff1")
  diff2_long <- tidyr::pivot_longer(diff2_counts_df[diff2_counts_df$event %in% events, ],
    cols=dplyr::all_of(names(sampleinfo)), names_to="sample", values_to="diff2")

  sjcnts <- cbind.data.frame(ref_long, diff1=diff1_long$diff1, diff2=diff2_long$diff2)
  sjcnts <- dplyr::mutate(sjcnts, groups=sampleinfo[sample])
  sjcnts$ref   <- as.numeric(sjcnts$ref)
  sjcnts$diff1 <- as.numeric(sjcnts$diff1)
  sjcnts$diff2 <- as.numeric(sjcnts$diff2)

  if (nrow(sjcnts) > 1L) {
    write.table(sjcnts,
      file=paste0(outdir, '/', gene_name, '.', outfilesuffix, '.sjcnt.txt'),
      quote=FALSE, sep="\t", row.names=FALSE)
  }
}


#' Count SJ reads for bipartition events in one gene
#'
#' For each bipartition row, sums junction read counts from a prebuilt SJ
#' matrix for the distinct1, distinct2, and shared intron sets, then writes
#' long-format output compatible with exontest.R.
#'
#' @param splits_df   Data frame of labeled bipartition rows for one gene
#'   (must contain intron_distinct1, intron_distinct2, intron_shared columns
#'   from \code{label_all_bipartition_introns}).
#' @param sj_mat      Combined junction count matrix (rows = chr:start:end keys,
#'   cols = samples with condition-prefixed names).
#' @param sampleinfo  Named character vector mapping sample names to groups.
#' @param outdir      Output directory.
#' @param outfilesuffix Suffix for output filenames.
#' @export
count_bipartition_sj <- function(splits_df, sj_mat, sampleinfo, outdir, outfilesuffix) {
  splits_df$event <- rownames(splits_df)
  splits_df <- splits_df[!(is.na(splits_df$intron_distinct1) &
                            is.na(splits_df$intron_distinct2)), ]
  if (nrow(splits_df) == 0L) return(invisible(NULL))

  gene       <- splits_df$gene[1L]
  n_samples  <- ncol(sj_mat)
  sample_names <- colnames(sj_mat)

  # replace NA with string "NA" so dict lookup works (same pattern as count_bipartitions)
  splits_df$intron_distinct1[is.na(splits_df$intron_distinct1)] <- "NA"
  splits_df$intron_distinct2[is.na(splits_df$intron_distinct2)] <- "NA"
  splits_df$intron_shared[is.na(splits_df$intron_shared)]       <- "NA"

  d1_dict <- t(sapply(unique(splits_df$intron_distinct1), sum_sj_counts, sj_mat=sj_mat))
  d2_dict <- t(sapply(unique(splits_df$intron_distinct2), sum_sj_counts, sj_mat=sj_mat))
  sh_dict <- t(sapply(unique(splits_df$intron_shared),    sum_sj_counts, sj_mat=sj_mat))

  diff1_mat <- as.data.frame(matrix(d1_dict[splits_df$intron_distinct1, ], nrow=nrow(splits_df), ncol=n_samples))
  diff2_mat <- as.data.frame(matrix(d2_dict[splits_df$intron_distinct2, ], nrow=nrow(splits_df), ncol=n_samples))
  ref_mat   <- as.data.frame(matrix(sh_dict[splits_df$intron_shared, ],    nrow=nrow(splits_df), ncol=n_samples))
  # NA shared means no shared introns in bubble -- treat as 0 counts, not unknown
  ref_mat[is.na(ref_mat)] <- 0L

  colnames(diff1_mat) <- sample_names
  colnames(diff2_mat) <- sample_names
  colnames(ref_mat)   <- sample_names
  rownames(diff1_mat) <- splits_df$event
  rownames(diff2_mat) <- splits_df$event
  rownames(ref_mat)   <- splits_df$event

  always_cols   <- c("gene","event","source","sink",
                     "intron_distinct1","intron_distinct2","intron_shared")
  optional_cols <- intersect(c("ref_ex_part","setdiff1","setdiff2",
                                "transcripts1","transcripts2","path1","path2"),
                             names(splits_df))
  meta_cols <- c(always_cols, optional_cols)
  diff1_mat$event <- splits_df$event
  diff2_mat$event <- splits_df$event
  ref_mat$event   <- splits_df$event

  diff1_mat <- dplyr::inner_join(splits_df[, meta_cols], diff1_mat, by="event")
  diff2_mat <- dplyr::inner_join(splits_df[, meta_cols], diff2_mat, by="event")
  ref_mat   <- dplyr::inner_join(splits_df[, meta_cols], ref_mat,   by="event")

  diff1_mean <- rowMeans(diff1_mat[, names(sampleinfo), drop=FALSE], na.rm=TRUE)
  diff2_mean <- rowMeans(diff2_mat[, names(sampleinfo), drop=FALSE], na.rm=TRUE)
  if (any(diff1_mean > 0, na.rm=TRUE) || any(diff2_mean > 0, na.rm=TRUE)) {
    write_sjcnt_wide(ref_mat, diff1_mat, diff2_mat,
      events=splits_df$event, gene_name=gene,
      sampleinfo=sampleinfo, outfilesuffix=outfilesuffix, outdir=outdir)
  }
  invisible(NULL)
}
