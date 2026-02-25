# Multinomial partition analysis functions for alternative splicing


#' Find reference exonic parts for multiple partitions (generalization of simple binary case)
#' @param source Character string. Name of the source vertex of the bubble.
#' @param sink Character string. Name of the sink vertex of the bubble.
#' @param ex_part_sets A list of character vectors, one per partition, each containing exonic part edge names for that partition.
#' @param tx_ex_part_sets A list of character vectors, one per partition, each containing exonic part edge names covered by the representative transcript.
#' @param g_expart An `igraph` graph of exonic parts.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
#' g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == "ex"])
#' ref <- grase::find_reference_exonic_part_general(
#'   source = "3", sink = "7",
#'   ex_part_sets = list(character(0), character(0)),
#'   tx_ex_part_sets = list(character(0), character(0)),
#'   g_expart = g_expart
#' )
#' }
find_reference_exonic_part_general <- function(source, sink, ex_part_sets, tx_ex_part_sets, g_expart) {
  ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
  names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")

  # Common parts: intersection across ALL partitions
  common_edges <- Reduce(intersect, ex_part_sets)
  common_exonic_part <- if(length(common_edges) > 0) common_edges[ex_or_in_gexpart[common_edges] == "ex_part"] else character(0)
  
  # Source adjacent parts: shared across ALL transcript sets
  incident_edges <- igraph::incident(g_expart, source, mode = 'in')
  incident_names = igraph::edge_attr(g_expart, 'name', incident_edges)
  source_adj_expart_name <- incident_names[ex_or_in_gexpart[incident_names] == "ex_part"]
  # Intersect with parts common to all transcript sets
  tx_ex_part_intersection <- Reduce(intersect, tx_ex_part_sets)
  source_adj_expart_name <- intersect(source_adj_expart_name, tx_ex_part_intersection)

  # Sink adjacent parts: shared across ALL transcript sets  
  incident_edges <- igraph::incident(g_expart, sink, mode = 'out')
  incident_names = igraph::edge_attr(g_expart, 'name', incident_edges)
  sink_adj_expart_name <- incident_names[ex_or_in_gexpart[incident_names] == "ex_part"]
  sink_adj_expart_name <- intersect(sink_adj_expart_name, tx_ex_part_intersection)
  
  return(list(common = common_exonic_part, source = source_adj_expart_name, sink = sink_adj_expart_name))
}

#' Find differential and reference exonic parts for multinomial partitions
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @param source Character string. Name of the source vertex of the bubble.
#' @param sink Character string. Name of the sink vertex of the bubble.
#' @param groups A list of integer vectors, each giving the path indices belonging to one group.
#' @param parsed_partitions A list of character vectors mapping each path index to its transcript names.
#' @param parsed_paths A list of character vectors of internal vertex names for each path.
#' @param tx_ex_parts A named list of exonic part edge names covered by each transcript.
#' @param multipaths_list A list accumulating result rows; new rows are appended and returned.
#' @param gene Character string. Gene identifier.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' g <- grase::set_edge_names(g)
#' result <- grase::multinomial_paths(
#'   gene = "ENSG00000000003.14", g = g, sg = NULL, outdir = tempdir()
#' )
#' }
find_diff_and_ref_exparts_general <- function(
    g, source, sink, groups, 
    parsed_partitions, parsed_paths, tx_ex_parts, 
    multipaths_list, gene=gene) {
  
  r <- length(groups)

  # Collect exonic part paths for each partition
  ex_parts <- list()
  vpaths_list <- list()
  tx_list <- list()

  ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
  g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part']) # have to delete all edges in one command
  g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex']) # have to delete all edges in one command
  ex_or_in_gexon <- igraph::edge_attr(g_exon, "ex_or_in")
  names(ex_or_in_gexon) <- igraph::edge_attr(g_exon, "name")
  ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
  names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")
  dexseq_frag <- igraph::edge_attr(g_expart, "dexseq_fragment")
  names(dexseq_frag) <- igraph::edge_attr(g_expart, "name")

  for (k in seq_len(r)) {
    txk <- unique(unlist(parsed_partitions[groups[[k]]]))
    vpaths <- lapply(parsed_paths[groups[[k]]], function(vec) c(source, vec, sink))
    epaths <- lapply(vpaths, from_vpath_to_exon_path_simple, g=g_exon)
    exparts_k <- lapply(epaths, from_epath_to_expart_path_simple, g=g_expart)
    exparts_k <- lapply(exparts_k, function(x) { return (x[ex_or_in_gexpart[x] == 'ex_part']) })

    ex_parts[[k]] <- exparts_k
    vpaths_list[[k]] <- vpaths
    tx_list[[k]] <- txk
  }

  # Compute I_k (intersection within each partition) and U_k (union)
  I <- lapply(ex_parts, intersect_stream) 
  U <- lapply(ex_parts, function(lst) unique(unlist(lst)))

  # Distinct per partition - generalize the binary logic:
  # For each group k, find what's common within that group but not in ANY other group
  D <- lapply(seq_len(r), function(k) {
    # Get union of all OTHER groups (equivalent to ex_part2_union in binary case)
    other_groups_union <- Reduce(union, U[-k])
    if (all(sapply(ex_parts[[k]], length) == 0)) {
      return("no_exon_part")
    }
    # Apply setdiff to each path in group k, then intersect results
    setdiff_results <- lapply(ex_parts[[k]], setdiff, other_groups_union)
    intersect_stream(setdiff_results)
  })

  # Get transcript exonic parts for each group
  tx_ex_part_sets <- lapply(seq_len(r), function(k) {
    txk <- tx_list[[k]]
    tx_ex_partk <- tx_ex_parts[txk[1]]  # Take first transcript as representative
    unique(unlist(tx_ex_partk))
  })

  # Find reference exonic parts using the generalized function
  ref_ex_part <- find_reference_exonic_part_general(
    source = source,
    sink = sink, 
    ex_part_sets = I,  # Use intersections within each group
    tx_ex_part_sets = tx_ex_part_sets,
    g_expart = g_expart
  )

  ref_ex_part_set <- unique(unlist(c(ref_ex_part$common, ref_ex_part$source, ref_ex_part$sink)))

  format_edges <- function(edges) {
    frags <- sort(as.numeric(dexseq_frag[edges]))
    if (length(frags) > 0) paste0("E", paste(sprintf("%03d", frags), collapse=",E")) else ""
  }

  # Build new row as a list (avoid rbind issues with variable columns)
  new_row <- list(
    gene = gene,
    source = source,
    sink = sink,
    ref_ex_part = format_edges(ref_ex_part_set)
  )

  # Add variable number of columns
  for (k in seq_len(r)) {
    if (identical(D[[k]], "no_exon_part")) {
      new_row[[paste0("setdiff", k)]] <- "no_exon_part"
    } else if (length(D[[k]]) == 0) {
      new_row[[paste0("setdiff", k)]] <- "{}"
    } else {
      new_row[[paste0("setdiff", k)]] <- format_edges(if (length(D[[k]]) > 0) D[[k]] else "empty")
    }
    new_row[[paste0("transcripts", k)]] <- paste(tx_list[[k]], collapse=", ")
    new_row[[paste0("path", k)]] <- paste(
      sapply(vpaths_list[[k]], function(vec) paste(vec, collapse="-")),
      collapse=", "
    )
    new_row[[paste0("exparts", k)]] <- paste(
      sapply(ex_parts[[k]], function(vec) if(length(vec) > 0) paste(vec, collapse="-") else "empty"),
      collapse=", "
    )
  }

  # Append to list instead of rbind
  multipaths_list <- append(multipaths_list, list(new_row))

  return(list(
    multipaths_list = multipaths_list,
    shared = ref_ex_part_set,
    distinct = D
  ))
}

#' Compute diff exons from graph using multinomial partitions (each path as its own group)
#' @param gene Character string. Gene identifier (e.g., Ensembl gene ID).
#' @param g An `igraph` directed acyclic graph representing a gene splicing graph.
#' @param sg Unused; retained for API compatibility. Pass \code{NULL}.
#' @param outdir Character string. Path to output directory.
#' @param max_path Integer. Maximum number of alternative paths in a bubble to process; bubbles
#'   with more paths are skipped. Default is \code{20}.
#' @param max_span Numeric. Maximum topological span (number of nodes between source and sink) of
#'   a bubble to process; bubbles exceeding this are skipped. Default is \code{Inf}.
#' @param collapse_bubbles Logical. Whether to collapse processed bubbles by removing redundant
#'   paths from the graph. Default is \code{FALSE}.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' result <- grase::multinomial_paths(
#'   gene = "ENSG00000000003.14", g = g, sg = NULL, outdir = tempdir(), max_path = 20
#' )
#' head(result[, c("gene", "source", "sink", "setdiff1", "setdiff2", "ref_ex_part")])
#' }
multinomial_paths <- function(gene, g, sg, outdir, max_path = 20, max_span = Inf, collapse_bubbles=FALSE) {

  multipaths_list <- list()
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
    g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part']) # have to delete all edges in one command
    g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex']) # have to delete all edges in one command

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
    tx_exonpaths = lapply(txpaths, grase::from_vpath_to_exon_path_simple, g=g_exon)         # high_mem
    tx_ex_parts = lapply(tx_exonpaths, grase::from_epath_to_expart_path_simple, g=g_expart) #high_mem

    parsed_paths <- lapply(bubbles_ordered$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_ordered$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])

    vpaths = lapply(parsed_paths, function(vec) c(source, vec, sink))
    epaths = lapply(vpaths, grase::from_vpath_to_exon_path_simple, g=g_exon) 

    ex_part_paths <- lapply(epaths, grase::from_epath_to_expart_path_simple, g=g_expart)
    ex_part_paths <- lapply(ex_part_paths, function(x) {return (x[ex_or_in_gexpart[x] == 'ex_part'])})

    names(ex_part_paths)  = 1:length(ex_part_paths) 

    log_debug(paste("ex_part_paths len:", length(ex_part_paths)))
    
    # For multinomial: each path is its own group
    n_paths <- length(parsed_partitions)
    if (n_paths < 2) {
        log_debug("Need at least 2 paths for multinomial analysis")
        next
    }
    
    # Create groups where each path index forms its own group
    groups <- lapply(1:n_paths, function(i) i)
    
    # Call the general function with multinomial groups
    retval = grase::find_diff_and_ref_exparts_general(
      g = g, 
      source = source, 
      sink = sink, 
      groups = groups, 
      parsed_partitions = parsed_partitions, 
      parsed_paths = parsed_paths, 
      tx_ex_parts = tx_ex_parts, 
      multipaths_list = multipaths_list, 
      gene = gene
    )
    multipaths_list = retval$multipaths_list

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
        # have to check if there are any other transcripts using these edges outside of the bubble. 
        # if so cannot remove the path, if not can remove.
        all_vec <- sort(unique(unlist(tx_list)))
        parsed_part <- sort(parsed_partitions[[idx]])
        if (identical(all_vec, parsed_part)) {
          # all the edges are covered by same transcripts. no partial coverage. can remove the edge.
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
      g <- igraph::delete_edges(g, edges_to_remove) # have to delete all edges in one command
      epaths = lapply(vpaths, grase::from_vpath_to_exon_path, g=g) 
      epaths_to_keep = epaths[idx_to_keep]
      g <- grase::update_txpaths_after_bubble_collapse2(g, tx_to_update, epaths_to_keep) 

      # update bubbles_df
      g <- grase::set_txpath_to_vertex_attr(g)
      bubbles_updated <- grase::detect_bubbles_igraph(g)            # high_mem
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
  if (length(multipaths_list) > 0) {
    multipaths_df <- standardize_multinomial_columns(multipaths_list, max_path)
    
    multipaths_df <- multipaths_df %>% dplyr::rowwise() %>% dplyr::mutate(ref_part_cnt = count_items(ref_ex_part)) %>% dplyr::ungroup()
    # remove rows where all setdiff columns are empty
    setdiff_cols <- grep("^setdiff[0-9]+$", names(multipaths_df), value = TRUE)
    path_cols <- grep("^path[0-9]+$", names(multipaths_df), value = TRUE)

    if (length(setdiff_cols) > 0) {
      # Pass the entire dataframe to apply so we can access both path and setdiff columns
      keep_rows <- apply(multipaths_df, 1, function(row) {
        # Get values for paths that exist (non-empty string in standardized df)
        p_vals <- row[path_cols]
        active_mask <- p_vals != "" & !is.na(p_vals)
        if (!any(active_mask)) return(FALSE)

        # Get corresponding setdiff columns
        active_path_cols <- names(p_vals)[active_mask]
        active_setdiff_cols <- gsub("path", "setdiff", active_path_cols)

        # Check if corresponding setdiff values are valid
        # We want to filter out "{}", but keep "no_exon_part" and valid edge IDs (starting with "E")
        sd_vals <- row[active_setdiff_cols]
        all(sd_vals != "{}" & !is.na(sd_vals))
      })
      multipaths_filtered <- multipaths_df[keep_rows, ]
    } else {
      multipaths_filtered <- multipaths_df
    }
     
    # For multinomial, we can't use the same grouping logic as binary case
    # Just remove duplicates based on all setdiff columns and ref_ex_part
    if (nrow(multipaths_filtered) > 0 && length(setdiff_cols) > 0) {
      multipaths_filtered <- multipaths_filtered %>%
        dplyr::distinct(dplyr::across(all_of(c("ref_ex_part", setdiff_cols))), .keep_all = TRUE) %>%
        dplyr::group_by(dplyr::across(all_of(setdiff_cols))) %>%
        dplyr::filter(ref_part_cnt == min(ref_part_cnt)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-ref_part_cnt)
    } else if (nrow(multipaths_filtered) > 0) {
      multipaths_filtered <- multipaths_filtered %>% dplyr::select(-ref_part_cnt)
    }
  } else {
    # Handle empty list case
    multipaths_df <- data.frame()
    multipaths_filtered <- data.frame()
  }

  filename_all = file.path(outdir, paste0(gene, ".multinomial.all.txt"))
  filename = file.path(outdir, paste0(gene, ".multinomial.txt"))
  log_debug(paste0("filename :", filename))

  write.table(multipaths_df, filename_all , sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(multipaths_filtered, filename, sep = "\t", quote = FALSE, row.names = FALSE)

  multipaths_filtered
}

#' Standardize multinomial results from list to data frame with consistent column structure
#' @param row_list List of row-lists with variable multinomial columns
#' @param max_paths Maximum number of paths to standardize to
#' @export
#' @examples
#' row_list <- list(
#'   list(gene="GENE1", source="3", sink="7", ref_ex_part="E001",
#'        setdiff1="E002", transcripts1="TX1", path1="3-5-7",
#'        setdiff2="E003", transcripts2="TX2", path2="3-6-7"),
#'   list(gene="GENE1", source="3", sink="9", ref_ex_part="E004",
#'        setdiff1="E005", transcripts1="TX1", path1="3-8-9",
#'        setdiff2="",     transcripts2="TX3", path2="3-9")
#' )
#' df <- grase::standardize_multinomial_columns(row_list, max_paths = 2)
#' colnames(df)
standardize_multinomial_columns <- function(row_list, max_paths) {
  if (length(row_list) == 0) return(data.frame())
  
  # Find the maximum number of paths across all rows
  max_paths_found <- 0
  for (row in row_list) {
    # Extract numeric suffixes from column names
    col_names <- names(row)
    setdiff_nums <- as.numeric(gsub("setdiff", "", col_names[grepl("^setdiff[0-9]+$", col_names)]))
    transcripts_nums <- as.numeric(gsub("transcripts", "", col_names[grepl("^transcripts[0-9]+$", col_names)]))
    path_nums <- as.numeric(gsub("path", "", col_names[grepl("^path[0-9]+$", col_names)]))
    
    row_max <- max(c(setdiff_nums, transcripts_nums, path_nums, 0), na.rm = TRUE)
    max_paths_found <- max(max_paths_found, row_max)
  }
  
  # Use the larger of max_paths or max_paths_found
  target_max <- max(max_paths_found, max_paths)
  
  # Define all expected columns
  base_cols <- c("gene", "source", "sink", "ref_ex_part")
  dynamic_cols <- c()
  for (i in 1:target_max) {
    dynamic_cols <- c(dynamic_cols, 
                     paste0("setdiff", i),
                     paste0("transcripts", i), 
                     paste0("path", i))
  }
  all_expected_cols <- c(base_cols, dynamic_cols)
  
  # Convert each row to have all expected columns
  standardized_rows <- lapply(row_list, function(row) {
    # Create a list with all expected columns, defaulting to empty string
    std_row <- setNames(rep("", length(all_expected_cols)), all_expected_cols)
    
    # Fill in values that exist in the original row
    for (col_name in names(row)) {
      if (col_name %in% all_expected_cols) {
        std_row[[col_name]] <- row[[col_name]]
      }
    }
    
    std_row
  })
  
  # Convert list of standardized rows to data frame
  df <- do.call(rbind, lapply(standardized_rows, function(x) as.data.frame(as.list(x), stringsAsFactors = FALSE)))
  
  # Ensure column order
  df[all_expected_cols]
}
