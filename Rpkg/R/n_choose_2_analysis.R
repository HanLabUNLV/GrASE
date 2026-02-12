# N-choose-2 pairwise analysis functions for alternative splicing


#' Find differential and reference exonic parts for pairwise combinations
#' @export
find_diff_and_ref_exparts_pairwise <- function(
    g, source, sink, path_pairs, 
    parsed_partitions, parsed_paths, tx_ex_parts, 
    pairwise_list, gene=gene) {
  
  # Setup graph components
  ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
  g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part'])
  g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex'])
  ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
  names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")
  dexseq_frag <- igraph::edge_attr(g_expart, "dexseq_fragment")
  names(dexseq_frag) <- igraph::edge_attr(g_expart, "name")

  format_edges <- function(edges) {
    if (length(edges) > 0) paste0("E", paste(dexseq_frag[edges], collapse=",E")) else ""
  }

  # Process each pairwise combination
  for (pair in path_pairs) {
    path1_idx <- pair[1]
    path2_idx <- pair[2]
    
    # Get transcript groups for this pair
    tx1 <- unique(unlist(parsed_partitions[path1_idx]))
    tx2 <- unique(unlist(parsed_partitions[path2_idx]))
    
    # Create vertex paths
    vpath1 <- list(c(source, parsed_paths[[path1_idx]], sink))
    vpath2 <- list(c(source, parsed_paths[[path2_idx]], sink))
    
    # Convert to exon paths
    epath1 <- lapply(vpath1, from_vpath_to_exon_path_simple, g=g_exon)
    epath2 <- lapply(vpath2, from_vpath_to_exon_path_simple, g=g_exon)
    
    # Convert to exonic part paths
    ex_part1 <- lapply(epath1, from_epath_to_expart_path_simple, g=g_expart)
    ex_part1 <- lapply(ex_part1, function(x) x[ex_or_in_gexpart[x] == 'ex_part'])
    ex_part1_intersect <- intersect_stream(ex_part1)
    ex_part1_union <- unique(unlist(ex_part1))

    ex_part2 <- lapply(epath2, from_epath_to_expart_path_simple, g=g_expart)
    ex_part2 <- lapply(ex_part2, function(x) x[ex_or_in_gexpart[x] == 'ex_part'])
    ex_part2_intersect <- intersect_stream(ex_part2)
    ex_part2_union <- unique(unlist(ex_part2))
    
    # Find differential parts (using binary logic)
    setdiff1 <- lapply(ex_part1, setdiff, ex_part2_union)
    setdiff2 <- lapply(ex_part2, setdiff, ex_part1_union)
    setdiff1 <- intersect_stream(setdiff1)
    setdiff2 <- intersect_stream(setdiff2)

    # Get transcript exonic parts
    tx_ex_part1 <- tx_ex_parts[tx1[1]]
    tx_ex_part1 <- unique(unlist(tx_ex_part1))
    tx_ex_part2 <- tx_ex_parts[tx2[1]] 
    tx_ex_part2 <- unique(unlist(tx_ex_part2))

    # Find reference exonic parts
    ref_ex_part <- find_reference_exonic_part_simple(
      source = source, 
      sink = sink, 
      ex_part1_set = ex_part1_intersect, 
      ex_part2_set = ex_part2_intersect, 
      tx_ex_part1_set = tx_ex_part1, 
      tx_ex_part2_set = tx_ex_part2, 
      g = g_expart 
    )
    
    ref_ex_part_set <- unique(unlist(c(ref_ex_part$common, ref_ex_part$source, ref_ex_part$sink)))

    # Build new row as a list
    new_row <- list(
      gene = gene,
      source = source,
      sink = sink,
      path1_idx = path1_idx,
      path2_idx = path2_idx,
      ref_ex_part = format_edges(ref_ex_part_set),
      setdiff1 = format_edges(setdiff1),
      setdiff2 = format_edges(setdiff2),
      transcripts1 = paste(tx1, collapse=", "),
      transcripts2 = paste(tx2, collapse=", "),
      path1 = paste(sapply(vpath1, function(vec) paste(vec, collapse = "-")), collapse=", "),
      path2 = paste(sapply(vpath2, function(vec) paste(vec, collapse = "-")), collapse=", ")
    )

    # Append to list
    pairwise_list <- append(pairwise_list, list(new_row))
  }

  return(list(
    pairwise_list = pairwise_list
  ))
}

#' Generate all n-choose-2 combinations of path indices
generate_path_pairs <- function(n_paths) {
  if (n_paths < 2) return(list())
  
  pairs <- list()
  for (i in 1:(n_paths-1)) {
    for (j in (i+1):n_paths) {
      pairs <- append(pairs, list(c(i, j)))
    }
  }
  pairs
}

#' Convert pairwise analysis results to standardized data frame
#' @export
standardize_pairwise_columns <- function(pairwise_list) {
  if (length(pairwise_list) == 0) return(data.frame())
  
  # Define expected columns (fixed structure for pairwise analysis)
  expected_cols <- c("gene", "source", "sink", "path1_idx", "path2_idx", 
                    "ref_ex_part", "setdiff1", "setdiff2", 
                    "transcripts1", "transcripts2", "path1", "path2")
  
  # Convert each row to have all expected columns
  standardized_rows <- lapply(pairwise_list, function(row) {
    # Create a list with all expected columns, defaulting to empty string
    std_row <- setNames(rep("", length(expected_cols)), expected_cols)
    
    # Fill in values that exist in the original row
    for (col_name in names(row)) {
      if (col_name %in% expected_cols) {
        std_row[[col_name]] <- row[[col_name]]
      }
    }
    
    std_row
  })
  
  # Convert list of standardized rows to data frame
  df <- do.call(rbind, lapply(standardized_rows, function(x) as.data.frame(as.list(x), stringsAsFactors = FALSE)))
  
  # Ensure column order
  df[expected_cols]
}

#' Compute diff exons from graph using n-choose-2 pairwise analysis
#' @export
n_choose_2_paths <- function(gene, g, sg, outdir, max_path = 20, max_span = Inf, collapse_bubbles=FALSE) {

  pairwise_list <- list()
  g <- grase::set_edge_names(g)

  txpaths <- grase::txpath_from_edgeattr(g)
  trans <- grase::transcripts_from_igraph(g)
  bubbles_df <- grase::detect_bubbles_igraph(g)     # high_mem
  if ( nrow(bubbles_df) == 0) {
    message(paste("single transcript, no bubbles ", gene))
    return (NULL) 
  }

  bubbles_ordered = grase::bubble_ordering4(g, bubbles_df)
  bubbles_orig = bubbles_ordered
  #bubbles_ordered = bubbles_ordered[,1:ncol(bubbles_df)]

  for (bubble_idx in 1:nrow(bubbles_ordered)) {
  
    ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
    g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part'])
    g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex'])

    ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
    names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")

    if (bubble_idx > nrow(bubbles_ordered)) {
        log_debug("collapsed all bubbles. exiting the loop 1")
        break
    }
    source <- bubbles_ordered$source[[bubble_idx]]
    sink <- bubbles_ordered$sink[[bubble_idx]]
    span <- bubbles_ordered$span[[bubble_idx]]
    num_paths <- bubbles_ordered$num_paths[[bubble_idx]]
    if (source == "R" && sink == "L") next

    if (num_paths > max_path || span > max_span) {
        message(paste("bubble with ", num_paths, " alternative paths and span ", span, " exonic parts: skipping ", gene, source, sink))
        next;
    }
    log_debug(paste("source: ", source, " sink: ", sink))

    txpaths <- grase::txpath_from_edgeattr(g)
    tx_exonpaths = lapply(txpaths, grase::from_vpath_to_exon_path_simple, g=g_exon)
    tx_ex_parts = lapply(tx_exonpaths, grase::from_epath_to_expart_path_simple, g=g_expart)

    parsed_paths <- lapply(bubbles_ordered$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_ordered$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])

    vpaths = lapply(parsed_paths, function(vec) c(source, vec, sink))
    epaths = lapply(vpaths, grase::from_vpath_to_exon_path_simple, g=g_exon) 

    ex_part_paths <- lapply(epaths, grase::from_epath_to_expart_path_simple, g=g_expart)
    ex_part_paths <- lapply(ex_part_paths, function(x) {return (x[ex_or_in_gexpart[x] == 'ex_part'])})

    names(ex_part_paths)  = 1:length(ex_part_paths) 

    log_debug(paste("ex_part_paths len:", length(ex_part_paths)))
    

    # For n-choose-2: generate all pairwise combinations
    n_paths <- length(parsed_partitions)
    if (n_paths < 2) {
        log_debug("Need at least 2 paths for pairwise analysis")
        next
    }
    
    # Generate all n-choose-2 combinations
    path_pairs <- generate_path_pairs(n_paths)
    log_debug(paste("Generated", length(path_pairs), "pairwise combinations for", n_paths, "paths"))
    
    # Call the pairwise analysis function
    retval = find_diff_and_ref_exparts_pairwise(
      g = g, 
      source = source, 
      sink = sink, 
      path_pairs = path_pairs,
      parsed_partitions = parsed_partitions, 
      parsed_paths = parsed_paths, 
      tx_ex_parts = tx_ex_parts, 
      pairwise_list = pairwise_list, 
      gene = gene
    )
    pairwise_list = retval$pairwise_list

    if (collapse_bubbles) {
      if (source == "R" || sink == "L") next
      n_partitions = length(parsed_partitions)
      tx_counts = unlist(lapply(parsed_partitions, length))
      ordered_idx = order(tx_counts)
      txs_with_epath = grase::find_tx_with_epath(g, trans, epaths)
      idx_to_remove = c()
      idx_to_keep = c() 
      for (idx in ordered_idx) {
        tx_list = txs_with_epath[[idx]]
        all_vec <- sort(unique(unlist(tx_list)))
        parsed_part <- sort(parsed_partitions[[idx]])
        if (identical(all_vec, parsed_part)) {
          idx_to_remove = c(idx_to_remove, idx)
        }
        else {
          idx_to_keep = c(idx_to_keep, idx)
        }
      }
      if(length(idx_to_keep) == 0) {
        idx_to_keep = idx_to_remove[length(idx_to_remove)]
        idx_to_remove = idx_to_remove[-length(idx_to_remove)]
      }
      if ((length(idx_to_remove) == 0) || (length(idx_to_keep) == 0)) {
        log_debug("no edges to remove: next")
        next;
      }
      tx_to_update = txs_with_epath[idx_to_remove] 
      edges_to_remove = unique(unlist(epaths[idx_to_remove]))
      g <- igraph::delete_edges(g, edges_to_remove)
      epaths = lapply(vpaths, grase::from_vpath_to_exon_path, g=g) 
      epaths_to_keep = epaths[idx_to_keep]
      g <- grase::update_txpaths_after_bubble_collapse2(g, tx_to_update, epaths_to_keep) 

      # update bubbles_df
      g <- grase::set_txpath_to_vertex_attr(g)
      bubbles_updated <- grase::detect_bubbles_igraph(g)
      if (nrow(bubbles_updated) == 0) {
        log_debug("collapsed all bubbles. exiting the loop 2")
        break
      }
      bubbles_updated_ordered = grase::bubble_ordering4(g, bubbles_updated)     
      finished_sources = bubbles_ordered$source[(1:bubble_idx)]
      finished_sinks = bubbles_ordered$sink[(1:bubble_idx)]
      finished_bubbles_idx = (bubbles_updated_ordered$source %in% finished_sources) & (bubbles_updated_ordered$sink %in% finished_sinks)
      bubbles_updated_ordered = bubbles_updated_ordered[!finished_bubbles_idx,]
      bubbles_ordered = S4Vectors::rbind(bubbles_ordered[1:bubble_idx,],bubbles_updated_ordered)
    }
  }

  # Convert list to standardized data frame for final processing
  if (length(pairwise_list) > 0) {
    pairwise_df <- standardize_pairwise_columns(pairwise_list)
    
    pairwise_df <- pairwise_df %>% dplyr::rowwise() %>% dplyr::mutate(ref_part_cnt = count_items(ref_ex_part)) %>% dplyr::ungroup()
    
    # Remove rows where both setdiff1 and setdiff2 are empty
    pairwise_filtered <- pairwise_df[!((pairwise_df$setdiff1=="") & (pairwise_df$setdiff2=="")),]
    
    # For pairwise analysis, use the same grouping logic as binary case
    if (nrow(pairwise_filtered) > 0) {
      pairwise_filtered <- pairwise_filtered %>%
        dplyr::group_by(setdiff1, setdiff2) %>%
        dplyr::filter(ref_part_cnt == min(ref_part_cnt)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-ref_part_cnt)
    }
  } else {
    # Handle empty list case
    pairwise_df <- data.frame()
    pairwise_filtered <- data.frame()
  }

  filename_all = file.path(outdir, paste0(gene, ".n_choose_2.all.txt"))
  filename = file.path(outdir, paste0(gene, ".n_choose_2.txt"))
  log_debug(paste0("filename :", filename))

  write.table(pairwise_df, filename_all , sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(pairwise_filtered, filename, sep = "\t", quote = FALSE, row.names = FALSE)

  pairwise_filtered
}
