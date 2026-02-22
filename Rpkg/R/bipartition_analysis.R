# Binary partition analysis functions for alternative splicing


#' Find reference exonic part near bubble region between two partitions
#' @export
find_reference_exonic_part_simple <- function(source, sink, ex_part1_set, ex_part2_set, tx_ex_part1_set, tx_ex_part2_set, g_expart) {
  source_adj_exonic_part <- c()
  sink_adj_exonic_part <- c()
  ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
  names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")

  common_edges <- intersect(ex_part1_set, ex_part2_set)
  common_exonic_part <- if(length(common_edges) > 0) common_edges[ex_or_in_gexpart[common_edges] == "ex_part"] else character(0)
  #node_current_pos <- igraph::vertex_attr(g_expart, 'position', index=source) 
  incident_edges <- igraph::incident(g_expart, source, mode = 'in')
  incident_names = igraph::edge_attr(g_expart, 'name', incident_edges)
  source_adj_expart_name <- incident_names[ex_or_in_gexpart[incident_names] == "ex_part"]
  source_adj_expart_name <- intersect(source_adj_expart_name, intersect(tx_ex_part1_set, tx_ex_part2_set))

  #node_current_pos <- igraph::vertex_attr(g_expart, 'position', index=sink) 
  incident_edges <- igraph::incident(g_expart, sink, mode = 'out')
  incident_names = igraph::edge_attr(g_expart, 'name', incident_edges)
  sink_adj_expart_name <- incident_names[ex_or_in_gexpart[incident_names] == "ex_part"]
  sink_adj_expart_name <- intersect(sink_adj_expart_name, intersect(tx_ex_part1_set, tx_ex_part2_set))
  
  return(list(common = common_exonic_part, source = source_adj_expart_name, sink = sink_adj_expart_name ))
}

#' Find reference exonic part near bubble region (legacy version)
#' @export
find_reference_exonic_part <- function(source, sink, ex_part1_set, ex_part2_set, tx_ex_part1_set, tx_ex_part2_set, g) {
  source_adj_exonic_part <- c()
  sink_adj_exonic_part <- c()
  common_edges <- intersect(ex_part1_set, ex_part2_set)
  common_exonic_part <- common_edges[igraph::edge_attr(g)$ex_or_in[as.numeric(common_edges)] == "ex_part"]

  node_current_pos <- igraph::vertex_attr(g, 'position', index=source) 
  incident_edges <- igraph::incident(g, source, mode = 'in')
  source_adj_exonic_part <- incident_edges[igraph::edge_attr(g)$ex_or_in[as.numeric(incident_edges)] == "ex_part"]
  source_adj_exonic_part <- intersect(source_adj_exonic_part, intersect(tx_ex_part1_set, tx_ex_part2_set))

  node_current_pos <- igraph::vertex_attr(g, 'position', index=sink) 
  incident_edges <- igraph::incident(g, sink, mode = 'out')
  sink_adj_exonic_part <- incident_edges[igraph::edge_attr(g)$ex_or_in[as.numeric(incident_edges)] == "ex_part"]
  sink_adj_exonic_part <- intersect(sink_adj_exonic_part, intersect(tx_ex_part1_set, tx_ex_part2_set))

  return(list(common = common_exonic_part, source = source_adj_exonic_part, sink = sink_adj_exonic_part))
}

#' Find differential and reference exonic parts for binary splits
#' @export
find_diff_and_ref_exparts_for_split <- function(g, source, sink, split, parsed_partitions, parsed_paths, tx_ex_parts, bipartitions_df, gene=gene,
                                                  g_exon=NULL, g_expart=NULL, ex_or_in_gexpart=NULL, dexseq_frag=NULL) {

      group1 = as.numeric(split$group1)
      group2 = as.numeric(split$group2)
      tx1 = unique(unlist( parsed_partitions[group1]))
      tx2 = unique(unlist( parsed_partitions[group2]))
      log_debug(tx1)
      log_debug(tx2)

      # Use precomputed graph derivatives if provided, otherwise compute them
      if (is.null(g_exon) || is.null(g_expart)) {
        ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
        g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part'])
        g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex'])
      }
      if (is.null(ex_or_in_gexpart)) {
        ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
        names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")
      }
      if (is.null(dexseq_frag)) {
        dexseq_frag <- igraph::edge_attr(g_expart, "dexseq_fragment")
        names(dexseq_frag) <- igraph::edge_attr(g_expart, "name")
      }

      vpath1 = lapply(parsed_paths[group1], function(vec) c(source, vec, sink)) 
      vpath2 = lapply(parsed_paths[group2], function(vec) c(source, vec, sink)) 
   
      #epath1.0 = lapply(vpath1, from_vpath_to_exon_path, g=g)     # high_mem 
      #epath2.0 = lapply(vpath2, from_vpath_to_exon_path, g=g)     # high_mem 
      epath1 = lapply(vpath1, from_vpath_to_exon_path_simple, g=g_exon)     # high_mem 
      epath2 = lapply(vpath2, from_vpath_to_exon_path_simple, g=g_exon)     # high_mem 
     
      #ex_part1.0 <- lapply(epath1.0, from_exon_path_to_exonic_part_path, g=g)   # high_mem
      ex_part1 <- lapply(epath1, from_epath_to_expart_path_simple, g=g_expart)   # high_mem
      ex_part1 <- lapply(ex_part1, function(x) {return (x[ex_or_in_gexpart[x] == 'ex_part'])})
      ex_part1_intersect <- intersect_stream(ex_part1)
      ex_part1_union <- unique(unlist(ex_part1))

      #ex_part2.0 <- lapply(epath2.0, from_exon_path_to_exonic_part_path, g=g)   # high_mem 
      ex_part2 <- lapply(epath2, from_epath_to_expart_path_simple, g=g_expart)   # high_mem
      ex_part2 <- lapply(ex_part2, function(x) {return (x[ex_or_in_gexpart[x] == 'ex_part'])})
      ex_part2_intersect <- intersect_stream(ex_part2)
      ex_part2_union <- unique(unlist(ex_part2))
     
      setdiff1 <- lapply(ex_part1, setdiff, ex_part2_union)
      setdiff2 <- lapply(ex_part2, setdiff, ex_part1_union)
      setdiff1 <- intersect_stream(setdiff1)
      setdiff2 <- intersect_stream(setdiff2)

      frag1 <- sort(as.numeric(dexseq_frag[setdiff1]))
      frag2 <- sort(as.numeric(dexseq_frag[setdiff2]))
      setdiff1ID <- if (length(frag1) > 0) paste0("E", paste(frag1, collapse=",E")) else ""
      setdiff2ID <- if (length(frag2) > 0) paste0("E", paste(frag2, collapse=",E")) else ""

      tx_ex_part1 <- tx_ex_parts[tx1[1]]
      tx_ex_part1 <- unique(unlist(tx_ex_part1))
      tx_ex_part2 <- tx_ex_parts[tx2[1]] 
      tx_ex_part2 <- unique(unlist(tx_ex_part2))

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
      ref_frag <- sort(as.numeric(dexseq_frag[ref_ex_part_set]))
      ref_ex_part_set_ID <- if (length(ref_frag) > 0) paste0("E", paste(ref_frag, collapse=",E")) else ""

      new_row <- data.frame(
        gene = gene,
        source = source,
        sink = sink,
        ref_ex_part = ref_ex_part_set_ID,
        setdiff1 = setdiff1ID,
        setdiff2 = setdiff2ID,
        transcripts1 = paste(tx1, collapse=", "),
        transcripts2 = paste(tx2, collapse=", "),
        path1 = paste(sapply(vpath1, function(vec) paste(vec, collapse = "-")), collapse=", "),
        path2 = paste(sapply(vpath2, function(vec) paste(vec, collapse = "-")), collapse=", ")
      )
      log_debug(new_row)
      bipartitions_df <- rbind(bipartitions_df, new_row)
  return(list(bipartitions_df = bipartitions_df, new_row = new_row))
}

#' Compute diff exons from graph and splicing structure for specific transcripts
#' @export
diff_exons_between_tx <- function(gene, g, tx_ids, outdir) {

  bipartitions_df <- data.frame()
  g <- grase::set_edge_names(g)
  bubbles_df <- grase::detect_bubbles_igraph(g)     # high_mem
  if ( nrow(bubbles_df) == 0) {
    message(paste("single transcript, no bubbles ", gene))
    return (NULL) 
  }

  for (bubble_idx in 1:nrow(bubbles_df)) {
    source <- bubbles_df$source[[bubble_idx]]
    sink <- bubbles_df$sink[[bubble_idx]]
    parsed_paths <- lapply(bubbles_df$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_df$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    contains_matrix <- do.call(rbind, lapply(parsed_partitions, function(x) tx_ids %in% x))
    
    if (all(!contains_matrix)) {
      next
    } else {
      rows_with_true <- apply(contains_matrix, 2, function(col) which(col))
      
      for (tx_idx in seq_along(tx_ids)) {
        if (length(rows_with_true[[tx_idx]]) > 0) {
          # This transcript is in some partitions
          for (i in rows_with_true[[tx_idx]]) {
            # Current partition
            current_partition <- parsed_partitions[[i]]
            current_path <- parsed_paths[[i]]
            
            # Find all other partitions
            other_partitions <- seq_along(parsed_partitions)[-i]
            other_tx <- unique(unlist(parsed_partitions[other_partitions]))
            other_paths <- parsed_paths[other_partitions]
            
            vpath1 = list(c(source, current_path, sink))
            vpath_rest = lapply(parsed_paths[seq_along(parsed_partitions)[-i]], function(vec) c(source, vec, sink)) 
            
            new_row <- data.frame(
              gene = gene,
              source = source,
              sink = sink,
              ref_ex_part = "",
              setdiff1 = "",
              setdiff2 = "",
              transcripts1 = tx_ids[tx_idx],
              transcripts2 = paste(other_tx, collapse=", "),
              path1 = paste(sapply(vpath1, function(vec) paste(vec, collapse = "-")), collapse=", "),
              path2 = paste(sapply(vpath_rest, function(vec) paste(vec, collapse = "-")), collapse=", ")
            )
            
            bipartitions_df <- rbind(bipartitions_df, new_row)
          }
        }
      }
    }
  }
  
  filename = file.path(outdir, paste0(gene, ".diff_tx.txt"))
  log_debug(paste0("filename :", filename))
  write.table(bipartitions_df, filename, sep = "\t", quote = FALSE, row.names = FALSE)
  bipartitions_df
}

#' Compute diff exons with one-vs-rest approach for each partition
#' @export
diff_exons_gene_ovr <- function(gene, g, outdir) {

  bipartitions_df <- data.frame()
  g <- grase::set_edge_names(g)
  bubbles_df <- grase::detect_bubbles_igraph(g)     # high_mem
  if ( nrow(bubbles_df) == 0) {
    message(paste("single transcript, no bubbles ", gene))
    return (NULL) 
  }

  for (bubble_idx in 1:nrow(bubbles_df)) {
    source <- bubbles_df$source[[bubble_idx]]
    sink <- bubbles_df$sink[[bubble_idx]]
    parsed_paths <- lapply(bubbles_df$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_df$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    #contains_matrix <- do.call(rbind, lapply(parsed_partitions, function(x) tx_ids %in% x))
    
    for (i in seq_along(parsed_partitions)) {
      # Current partition vs rest
      current_partition <- parsed_partitions[[i]]
      current_path <- parsed_paths[[i]]
      
      # All other partitions combined
      other_partitions <- seq_along(parsed_partitions)[-i]
      other_tx <- unique(unlist(parsed_partitions[other_partitions]))
      
      vpath1 = list(c(source, current_path, sink))
      vpath_rest = lapply(parsed_paths[seq_along(parsed_partitions)[-i]], function(vec) c(source, vec, sink)) 
      
      new_row <- data.frame(
        gene = gene,
        source = source,
        sink = sink,
        ref_ex_part = "",
        setdiff1 = "",
        setdiff2 = "",
        transcripts1 = paste(current_partition, collapse=", "),
        transcripts2 = paste(other_tx, collapse=", "),
        path1 = paste(sapply(vpath1, function(vec) paste(vec, collapse = "-")), collapse=", "),
        path2 = paste(sapply(vpath_rest, function(vec) paste(vec, collapse = "-")), collapse=", ")
      )
      
      bipartitions_df <- rbind(bipartitions_df, new_row)
    }
  }
  
  filename = file.path(outdir, paste0(gene, ".gene_ovr.txt"))
  log_debug(paste0("filename :", filename))
  write.table(bipartitions_df, filename, sep = "\t", quote = FALSE, row.names = FALSE)
  bipartitions_df
}

#' Compute diff exons from graph using binary partitions
#' For each bubble, find all valid binary partitions of paths, then find diff exons for each partition. 
#' Optionally collapse bubbles with redundant paths.
#' @export
bipartition_paths <- function(gene, g, outdir, max_path = 15, max_span = Inf, collapse_bubbles=FALSE) {

  bipartitions_df <- data.frame()
  g <- grase::set_edge_names(g)

  txpaths <- grase::txpath_from_edgeattr(g)
  trans <- grase::transcripts_from_igraph(g)
  bubbles_df <- grase::detect_bubbles_igraph(g)     # high_mem
  if ( nrow(bubbles_df) == 0) {
    message(paste("single transcript, no bubbles ", gene))
    return (NULL)
  }

  bubbles_ordered = grase::bubble_ordering4(g, bubbles_df)
  #bubbles_ordered = bubbles_ordered[,1:ncol(bubbles_df)]

  # Precompute graph derivatives (invariant when collapse_bubbles=FALSE)
  .recompute_graph_derivs <- function(g) {
    ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
    g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part'])
    g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex'])
    ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
    names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")
    dexseq_frag <- igraph::edge_attr(g_expart, "dexseq_fragment")
    names(dexseq_frag) <- igraph::edge_attr(g_expart, "name")
    txpaths <- grase::txpath_from_edgeattr(g)
    tx_exonpaths <- lapply(txpaths, grase::from_vpath_to_exon_path_simple, g=g_exon)
    tx_ex_parts <- lapply(tx_exonpaths, grase::from_epath_to_expart_path_simple, g=g_expart)
    list(g_exon=g_exon, g_expart=g_expart, ex_or_in_gexpart=ex_or_in_gexpart,
         dexseq_frag=dexseq_frag, txpaths=txpaths, tx_exonpaths=tx_exonpaths, tx_ex_parts=tx_ex_parts)
  }
  gd <- .recompute_graph_derivs(g)

  n_bubbles <- nrow(bubbles_ordered)
  message(sprintf("[%s] %s: processing %d bubbles", format(Sys.time(), "%H:%M:%S"), gene, n_bubbles))
  pb <- utils::txtProgressBar(min = 0, max = n_bubbles, style = 3)

  for (bubble_idx in 1:nrow(bubbles_ordered)) {

    utils::setTxtProgressBar(pb, bubble_idx)

    if (bubble_idx > nrow(bubbles_ordered)) {
        log_debug("collapsed all bubbles. exiting the loop 1")
        break
    }
    source <- bubbles_ordered$source[[bubble_idx]]
    sink <- bubbles_ordered$sink[[bubble_idx]]
    num_paths <- bubbles_ordered$num_paths[[bubble_idx]]
    span <- bubbles_ordered$span[[bubble_idx]] 
    if (source == "R" && sink == "L") next
    if (num_paths > max_path || span > max_span) {
        message(paste("bubble with ", num_paths, " alternative paths and span ", span, " exonic parts: skipping ", gene, source, sink))
        next;
    }
    log_debug(paste("source: ", source, " sink: ", sink))

    parsed_paths <- lapply(bubbles_ordered$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_ordered$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])

    vpaths = lapply(parsed_paths, function(vec) c(source, vec, sink))
    epaths = lapply(vpaths, grase::from_vpath_to_exon_path_simple, g=gd$g_exon) 

    ex_part_paths <- lapply(epaths, grase::from_epath_to_expart_path_simple, g=gd$g_expart)
    ex_part_paths <- lapply(ex_part_paths, function(x) {return (x[gd$ex_or_in_gexpart[x] == 'ex_part'])})

    names(ex_part_paths)  = 1:length(ex_part_paths) 

    log_debug(paste("ex_part_paths len:", length(ex_part_paths)))
    contain_dag <- grase::path_subset_relation(ex_part_paths) 
    log_debug(paste("contain_dag len:", length(igraph::E(contain_dag))))

    int_splits <- grase::valid_partitions_idx(contain_dag)
    nodes <- igraph::V(contain_dag)$name
    valid_splits <- lapply(int_splits, function(idx) {
      g1 <- sort(nodes[idx])
      g2 <- sort(nodes[-idx])
      list(group1 = g1, group2 = g2)
    })
    valid_splits <- grase::remove_symmetric_splits(valid_splits)
    rm(contain_dag)
    
    for (split in valid_splits) {
      retval = grase::find_diff_and_ref_exparts_for_split(
        g=g, source=source, sink=sink, split=split,
        parsed_partitions=parsed_partitions, parsed_paths=parsed_paths,
        tx_ex_parts=gd$tx_ex_parts, bipartitions_df=bipartitions_df, gene=gene,
        g_exon=gd$g_exon, g_expart=gd$g_expart,
        ex_or_in_gexpart=gd$ex_or_in_gexpart, dexseq_frag=gd$dexseq_frag
      )
      new_row <- retval$new_row
      # Early filter: skip if both setdiff columns are empty
      if (new_row$setdiff1 == "" && new_row$setdiff2 == "") next
      bipartitions_df = retval$bipartitions_df
    }

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

      # update bubbles_df and recompute graph derivatives after collapse
      g <- grase::set_txpath_to_vertex_attr(g)
      gd <- .recompute_graph_derivs(g)
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

  close(pb)
  log_debug(sprintf("[%s] %s: %d rows after empty-setdiff filter",
                  format(Sys.time(), "%H:%M:%S"), gene, nrow(bipartitions_df)))

  if (nrow(bipartitions_df) > 0) {
    bipartitions_df <- bipartitions_df %>% dplyr::mutate(ref_part_cnt = sapply(ref_ex_part, count_items))
    # deduplicate exact events (same ref_ex_part + setdiff1 + setdiff2 from different graph paths),
    # then among remaining rows with same (setdiff1, setdiff2) keep the one with fewest ref parts
    bipartitions_filtered <- bipartitions_df %>%
      dplyr::distinct(ref_ex_part, setdiff1, setdiff2, .keep_all = TRUE) %>%
      dplyr::group_by(setdiff1, setdiff2) %>%
      dplyr::filter(ref_part_cnt == min(ref_part_cnt)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-ref_part_cnt)
  } else {
    bipartitions_filtered <- bipartitions_df
  }

  filename_all = file.path(outdir, paste0(gene, ".bipartition.all.txt"))
  filename = file.path(outdir, paste0(gene, ".bipartition.txt"))
  log_debug(paste0("filename :", filename))

  write.table(bipartitions_df, filename_all , sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(bipartitions_filtered, filename, sep = "\t", quote = FALSE, row.names = FALSE)

  bipartitions_filtered
}

#' Infer differential exonic parts given two groups of transcripts (e.g. changed vs constant)
#' @export
infer_diff_exons_from_tx_groups <- function(gene, g, tx_changed, tx_constant, outdir) {
  
  bipartitions_df <- data.frame()
  g <- grase::set_edge_names(g)
  
  # Precompute graph derivatives
  ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
  g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part'])
  g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex'])
  ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
  names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")
  dexseq_frag <- igraph::edge_attr(g_expart, "dexseq_fragment")
  names(dexseq_frag) <- igraph::edge_attr(g_expart, "name")
  txpaths <- grase::txpath_from_edgeattr(g)
  tx_exonpaths <- lapply(txpaths, grase::from_vpath_to_exon_path_simple, g=g_exon)
  tx_ex_parts <- lapply(tx_exonpaths, grase::from_epath_to_expart_path_simple, g=g_expart)

  bubbles_df <- grase::detect_bubbles_igraph(g)
  if (nrow(bubbles_df) == 0) {
    message(paste("single transcript, no bubbles ", gene))
    return(NULL)
  }

  for (bubble_idx in 1:nrow(bubbles_df)) {
    source <- bubbles_df$source[[bubble_idx]]
    sink <- bubbles_df$sink[[bubble_idx]]
    parsed_paths <- lapply(bubbles_df$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_df$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    
    # Identify which partitions contain transcripts from each group
    idx_changed <- which(sapply(parsed_partitions, function(p) any(p %in% tx_changed)))
    idx_constant <- which(sapply(parsed_partitions, function(p) any(p %in% tx_constant)))
    
    if (length(idx_changed) == 0 || length(idx_constant) == 0) {
      next
    }
    
    split <- list(group1 = idx_changed, group2 = idx_constant)
    
    retval <- find_diff_and_ref_exparts_for_split(
      g=g, source=source, sink=sink, split=split,
      parsed_partitions=parsed_partitions, parsed_paths=parsed_paths,
      tx_ex_parts=tx_ex_parts, bipartitions_df=bipartitions_df, gene=gene,
      g_exon=g_exon, g_expart=g_expart,
      ex_or_in_gexpart=ex_or_in_gexpart, dexseq_frag=dexseq_frag
    )
    bipartitions_df <- retval$bipartitions_df
  }
  
  filename <- file.path(outdir, paste0(gene, ".inferred_diff.txt"))
  write.table(bipartitions_df, filename, sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(bipartitions_df)
}
