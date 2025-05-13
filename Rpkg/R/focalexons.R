
#' @importFrom magrittr %>%
NULL

#' Count items in comma-separated string
#' @export
count_items <- function(x) {
  if (x == "") return(0)
  length(stringr::str_split(x, ",")[[1]])
}

#' Find reference exonic part near bubble region
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



#' Compute focal exons from graph and splicing structure
#' @export
focal_exons_between_tx <- function(gene, g, sg, tx_ids, outdir) {

  bubbles_df <- SplicingGraphs::bubbles(sg)
  tx_ex_part1 <- from_exon_path_to_exonic_part_path(g, from_vpath_to_exon_path(g, as.character(SplicingGraphs::txpath(sg)[[tx_ids[1]]])))
  tx_ex_part2 <- from_exon_path_to_exonic_part_path(g, from_vpath_to_exon_path(g, as.character(SplicingGraphs::txpath(sg)[[tx_ids[2]]])))

  focalexons_df <- data.frame()
  for (bubble_idx in seq_along(bubbles_df$partitions)) {
    source <- bubbles_df$source[[bubble_idx]]
    sink <- bubbles_df$sink[[bubble_idx]]
    if (source == "R" && sink == "L") next

    parsed_paths <- lapply(bubbles_df$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_df$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    contains_matrix <- do.call(rbind, lapply(parsed_partitions, function(x) tx_ids %in% x))
    colnames(contains_matrix) <- tx_ids

    if (all(colSums(contains_matrix) > 0) && sum(rowSums(contains_matrix) > 0) > 1) {
      rows_with_true <- apply(contains_matrix, 2, function(col) which(col))
      vpath1 <- c(source, parsed_paths[[rows_with_true[[1]]]], sink)
      vpath2 <- c(source, parsed_paths[[rows_with_true[[2]]]], sink)
      epath1 <- from_vpath_to_exon_path(g, vpath1)
      epath2 <- from_vpath_to_exon_path(g, vpath2)
      ex_part1 <- from_exon_path_to_exonic_part_path(g, epath1)
      ex_part2 <- from_exon_path_to_exonic_part_path(g, epath2)

      setdiff1 <- setdiff(ex_part1, ex_part2)
      setdiff2 <- setdiff(ex_part2, ex_part1)
      setdiff1 <- setdiff1[igraph::edge_attr(g)$ex_or_in[setdiff1] == "ex_part"]
      setdiff2 <- setdiff2[igraph::edge_attr(g)$ex_or_in[setdiff2] == "ex_part"]

      setdiff1ID <- if (length(setdiff1) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[setdiff1], collapse=",E")) else ""
      setdiff2ID <- if (length(setdiff2) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[setdiff2], collapse=",E")) else ""

      ref_ex_part <- find_reference_exonic_part(
        source = source, 
        sink = sink, 
        ex_part1_set = ex_part1, 
        ex_part2_set = ex_part2, 
        tx_ex_part1_set = tx_ex_part1, 
        tx_ex_part2_set = tx_ex_part2, 
        g = g, 
      ) 
      ref_ex_part_set <- unique(c(ref_ex_part$common, ref_ex_part$source, ref_ex_part$sink))
      ref_ex_part_set_ID <- if (length(ref_ex_part_set) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[ref_ex_part_set], collapse=",E")) else ""

      new_row <- data.frame(
        gene = gene,
        ref_ex_part = ref_ex_part_set_ID,
        setdiff1 = setdiff1ID,
        setdiff2 = setdiff2ID,
        transcripts1 = tx_ids[1],
        transcripts2 = tx_ids[2],
        path1 = paste(vpath1, collapse=", "),
        path2 = paste(vpath2, collapse=", ")
      )
      focalexons_df <- rbind(focalexons_df, new_row)
    }
  }

  if (nrow(focalexons_df) > 0) {
    focalexons_df <- focalexons_df %>% dplyr::rowwise() %>% dplyr::mutate(adj_cnt = count_items(ref_ex_part)) %>% dplyr::ungroup()
    focalexons_filtered <- focalexons_df %>%
      dplyr::group_by(setdiff1, setdiff2) %>%
      dplyr::filter(adj_cnt == min(adj_cnt)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-adj_cnt)
  } else {
    focalexons_filtered <- focalexons_df
  }

  write.table(focalexons_df, file.path(outdir, "focalexons", paste0(gene, ".focalexons.all.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(focalexons_filtered, file.path(outdir, "focalexons", paste0(gene, ".focalexons.txt")), sep = "\t", quote = FALSE, row.names = FALSE)

  return(focalexons_filtered)
}



#' Compute focal exons from graph and splicing structure
#' @export
focal_exons_gene_ovr <- function(gene, g, sg, outdir) {

  bubbles_df <- SplicingGraphs::bubbles(sg)

  txpaths <- lapply(SplicingGraphs::txpath(sg), as.character)
  tx_exon_paths <- grase::txpath_from_edgeattr(g)
  tx_exonpaths = lapply(txpaths, from_vpath_to_exon_path, g=g)
  tx_ex_parts = lapply(tx_exonpaths, from_exon_path_to_exonic_part_path, g=g)

  focalexons_df <- data.frame()
  for (bubble_idx in seq_along(bubbles_df$partitions)) {
    source <- bubbles_df$source[[bubble_idx]]
    sink <- bubbles_df$sink[[bubble_idx]]
    if (source == "R" && sink == "L") next

    parsed_paths <- lapply(bubbles_df$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_df$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    #contains_matrix <- do.call(rbind, lapply(parsed_partitions, function(x) tx_ids %in% x))
    #colnames(contains_matrix) <- tx_ids

    for (i in 1:length(parsed_partitions)) {

      tx1 = parsed_partitions[i][[1]]
      tx_rest = unique(unlist( parsed_partitions[seq_along(parsed_partitions)[-i]]))
      print (tx1)

      vpath1 = c(source, parsed_paths[[i]], sink)
      vpath_rest = lapply(parsed_paths[seq_along(parsed_partitions)[-i]], function(vec) c(source, vec, sink)) 
   
      epath1 <- from_vpath_to_exon_path(g, vpath1)
      epath_rest = lapply(vpath_rest, from_vpath_to_exon_path, g=g) 
     
      ex_part1 <- from_exon_path_to_exonic_part_path(g, epath1) 
      ex_part_rest <- lapply(epath_rest, from_exon_path_to_exonic_part_path, g=g) 
      ex_part2 <- unique(unlist(ex_part_rest))
  
      setdiff1 <- setdiff(ex_part1, ex_part2)
      setdiff2 <- setdiff(ex_part2, ex_part1)
      setdiff1 <- setdiff1[igraph::edge_attr(g)$ex_or_in[setdiff1] == "ex_part"]
      setdiff2 <- setdiff2[igraph::edge_attr(g)$ex_or_in[setdiff2] == "ex_part"]

      setdiff1ID <- if (length(setdiff1) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[setdiff1], collapse=",E")) else ""
      setdiff2ID <- if (length(setdiff2) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[setdiff2], collapse=",E")) else ""

      tx_ex_part1 <- tx_ex_parts[tx1[1]]
      tx_ex_part_rest <- tx_ex_parts[tx_rest] 
      tx_ex_part2 <- unique(unlist(tx_ex_part_rest))

      ref_ex_part <- find_reference_exonic_part(
        source = source, 
        sink = sink, 
        ex_part1_set = ex_part1, 
        ex_part2_set = ex_part2, 
        tx_ex_part1_set = tx_ex_part1, 
        tx_ex_part2_set = tx_ex_part2, 
        g = g, 
      ) 
 
      ref_ex_part_set <- unique(unlist(c(ref_ex_part$common, ref_ex_part$source, ref_ex_part$sink)))
      ref_ex_part_set_ID <- if (length(ref_ex_part_set) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[ref_ex_part_set], collapse=",E")) else ""

      new_row <- data.frame(
        gene = gene,
        ref_ex_part = ref_ex_part_set_ID,
        setdiff1 = setdiff1ID,
        setdiff2 = setdiff2ID,
        transcripts1 = paste(tx1, collapse=", "),
        transcripts2 = paste(tx_rest, collapse=", "),
        path1 = paste(vpath1, collapse=", "),
        path2 = paste(sapply(vpath_rest, function(vec) paste(vec, collapse = "-")), collapse=", ")
      )
      print(new_row)
      focalexons_df <- rbind(focalexons_df, new_row)
    }
  }

  if (nrow(focalexons_df) > 0) {
    focalexons_df <- focalexons_df %>% dplyr::rowwise() %>% dplyr::mutate(adj_cnt = count_items(ref_ex_part)) %>% dplyr::ungroup()
    focalexons_filtered <- focalexons_df %>%
      dplyr::group_by(setdiff1, setdiff2) %>%
      dplyr::filter(adj_cnt == min(adj_cnt)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-adj_cnt)
  } else {
    focalexons_filtered <- focalexons_df
  }

  write.table(focalexons_df, file.path(outdir, "focalexons", paste0(gene, ".focalexons.all.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(focalexons_filtered, file.path(outdir, "focalexons", paste0(gene, ".focalexons.txt")), sep = "\t", quote = FALSE, row.names = FALSE)

}




get_interval <- function(bubble, topo_idx) {
  vals = c(topo_idx[bubble$source], topo_idx[bubble$sink])
  vals <- setNames(vals, c("start", "end"))
  return (vals)
}

get_depths <- function(intervals) {
  n <- length(intervals)
  depths <- rep(0, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) {
        a <- intervals[[i]]
        b <- intervals[[j]]
        if (b["start"] <= a["start"] && a["end"] <= b["end"]) {
          depths[i] <- depths[i] + 1 
        }   
      }   
    }   
  }
  return(depths)
}


#' hierarchical ordering of bubbles in the DAG grase graph from innermost to outermost based on nesting relationships. 
#' @export
bubble_ordering <- function(g, bubbles_df) {

  topo <- igraph::topo_sort(g, mode = "out")
  topo_idx <- setNames(seq_along(topo), names(topo))

  bubbles_list = split(as.data.frame(bubbles_df), seq(nrow(bubbles_df)))
  bubble_intervals <- lapply(bubbles_list, get_interval, topo_idx = topo_idx)
  bubble_depths <- get_depths(bubble_intervals)
  ordered_bubbles <- bubbles_list[order(-bubble_depths)]
  do.call(rbind.data.frame, ordered_bubbles)

}





# Find subset relationships
is_proper_subset <- function(a, b) {
  all(a %in% b) && length(a) < length(b)
}


#' find a containment DAG that represents the hierarchical relationship between paths 
#' @export
path_subset_relation <-function(sets) {

#ex_index = igraph::edge_attr(g)$ex_or_in == 'ex'
#g_exonic <- igraph::delete_edges(g_exonic, igraph::E(g_exonic)[ex_index])
#g_exonic
  n <- length(sets)
  names <- names(sets)
  edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j && is_proper_subset(sets[[i]], sets[[j]])) {
        # Ensure it's a direct subset (no intermediate sets)
        intermediate <- FALSE
        for (k in seq_len(n)) {
          if (k != i && k != j &&
              is_proper_subset(sets[[i]], sets[[k]]) &&
              is_proper_subset(sets[[k]], sets[[j]])) {
            intermediate <- TRUE
            break
          }
        }
        if (!intermediate) {
          edges <- rbind(edges, data.frame(from = names[j], to = names[i]))
        }
      }
    }
  }

  # Create graph with all vertices
  g <- igraph::graph_from_data_frame(edges, vertices = names(sets), directed = TRUE)
  return(g)
}


#' find valid partitions that represents the hierarchical relationship between paths 
#' @export
valid_partitions <- function(g, max_powerset) {

  nodes <- igraph::V(g)$name
  edges <- igraph::E(g)

  # Generate all non-empty, non-full subsets of nodes
  powerset <- function(x) {
    unlist(lapply(1:(length(x) - 1), function(k) combn(x, k, simplify = FALSE)), recursive = FALSE)
  }

  # Check if a subset is downward-closed
  is_downward_closed <- function(g, subset) {
    for (node in subset) {
      descendants <- igraph::subcomponent(g, node, mode = "out")$name
      if (!all(descendants %in% subset)) {
        return(FALSE)
      }
    }
    return(TRUE)
  }

  
  valid_splits <- list()
  # check size of powerset
  powernodes = powerset(nodes)
  print(paste("powerset: ", length(powernodes)))
  if (length(powernodes) > max_powerset) {
    print (paste("power set size greater than ", max_powerset, "skipping."))
    return (valid_splits)
  }
  # Compute all valid binary splits
  for (subset in powernodes) {
    if (is_downward_closed(g, subset)) {
      complement <- setdiff(nodes, subset)
      valid_splits[[length(valid_splits) + 1]] <- list(group1 = sort(subset), group2 = sort(complement))
    }
  }

  # Print splits
  for (i in seq_along(valid_splits)) {
    cat(sprintf("Split %d:\n", i))
    cat("  Group 1:", paste(valid_splits[[i]]$group1, collapse = ", "), "\n")
    cat("  Group 2:", paste(valid_splits[[i]]$group2, collapse = ", "), "\n\n")
  }
  return (valid_splits)

}


#' Compute focal exons from graph and splicing structure
#' @export
focal_exons_gene_powerset <- function(gene, g, sg, outdir, max_powerset = 10000) {


  focalexons_df <- data.frame()

  txpaths <- grase::txpath_from_edgeattr(g)
  txmat <- grase::make_matrix_from_txpath_igraph(txpaths)
  bubble_df <- grase::detect_bubbles_from_mat(txmat)

  if ( nrow(bubbles_df) == 0) {
    print("single transcript, no bubbles")
    return (NULL) 
  }

  #txpaths <- lapply(SplicingGraphs::txpath(sg), as.character)
  tx_exonpaths = lapply(txpaths, from_vpath_to_exon_path, g=g)
  tx_ex_parts = lapply(tx_exonpaths, from_exon_path_to_exonic_part_path, g=g)

  bubbles_ordered = bubble_ordering(g, bubbles_df)
  for (bubble_idx in seq_along(bubbles_ordered$partitions)) {
    source <- bubbles_ordered$source[[bubble_idx]]
    sink <- bubbles_ordered$sink[[bubble_idx]]
    if (source == "R" && sink == "L") next

    parsed_paths <- lapply(bubbles_ordered$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_ordered$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    #contains_matrix <- do.call(rbind, lapply(parsed_partitions, function(x) tx_ids %in% x))
    #colnames(contains_matrix) <- tx_ids


    #vpath1 = c(source, parsed_paths[[i]], sink)
    #vpath_rest = lapply(parsed_paths[seq_along(parsed_partitions)[-i]], function(vec) c(source, vec, sink)) 
    vpaths = lapply(parsed_paths, function(vec) c(source, vec, sink))


    #epath1 <- from_vpath_to_exon_path(g, vpath1)
    #epath_rest = lapply(vpath_rest, from_vpath_to_exon_path, g=g) 
    epaths = lapply(vpaths, from_vpath_to_exon_path, g=g) 

   
    #ex_part1 <- from_exon_path_to_exonic_part_path(g, epath1) 
    #ex_part_rest <- lapply(epath_rest, from_exon_path_to_exonic_part_path, g=g) 
    ex_part_paths <- lapply(epaths, from_exon_path_to_exonic_part_path, g=g)
    ex_index = igraph::edge_attr(g)$ex_or_in == 'ex_part'
    ex_part_paths <- lapply(ex_part_paths, function(x) {return (x[igraph::edge_attr(g)$ex_or_in[x] == 'ex_part'])})
    names(ex_part_paths)  = sapply(parsed_partitions, function(vec) paste(vec, collapse=","))
    names(ex_part_paths)  = 1:length(ex_part_paths) 
    #ex_part2 <- unique(unlist(ex_part_rest))

    print (paste("ex_part_paths len:", length(ex_part_paths)))
    contain_dag <- path_subset_relation(ex_part_paths) 
    print (paste("contain_dag len:", length(igraph::E(g))))
    valid_splits <- valid_partitions(contain_dag, max_powerset)
    rm(contain_dag)
    
    for (split in valid_splits) {

      group1 = as.numeric(split$group1)
      group2 = as.numeric(split$group2)
      tx1 = unique(unlist( parsed_partitions[seq_along(parsed_partitions)[group1]]))
      tx2 = unique(unlist( parsed_partitions[seq_along(parsed_partitions)[group2]]))
      print (tx1)
      print (tx2)

      vpath1 = lapply(parsed_paths[seq_along(parsed_partitions)[group1]], function(vec) c(source, vec, sink)) 
      vpath2 = lapply(parsed_paths[seq_along(parsed_partitions)[group2]], function(vec) c(source, vec, sink)) 
   
      epath1 = lapply(vpath1, from_vpath_to_exon_path, g=g) 
      epath2 = lapply(vpath2, from_vpath_to_exon_path, g=g) 
     
      ex_part1 <- lapply(epath1, from_exon_path_to_exonic_part_path, g=g) 
      ex_part1 <- lapply(ex_part1, function(x) {return (x[igraph::edge_attr(g)$ex_or_in[x] == 'ex_part'])})
      ex_part1_set <- Reduce(intersect, ex_part1)
      ex_part2 <- lapply(epath2, from_exon_path_to_exonic_part_path, g=g) 
      ex_part2 <- lapply(ex_part2, function(x) {return (x[igraph::edge_attr(g)$ex_or_in[x] == 'ex_part'])})
      ex_part2_set <- Reduce(intersect, ex_part2)
  
      setdiff1 <- setdiff(ex_part1_set, ex_part2_set)
      setdiff2 <- setdiff(ex_part2_set, ex_part1_set)

      setdiff1ID <- if (length(setdiff1) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[setdiff1], collapse=",E")) else ""
      setdiff2ID <- if (length(setdiff2) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[setdiff2], collapse=",E")) else ""

      tx_ex_part1 <- tx_ex_parts[tx1[1]]
      tx_ex_part1 <- unique(unlist(tx_ex_part1))
      tx_ex_part2 <- tx_ex_parts[tx2[1]] 
      tx_ex_part2 <- unique(unlist(tx_ex_part2))

      ref_ex_part <- find_reference_exonic_part(
        source = source, 
        sink = sink, 
        ex_part1_set = ex_part1_set, 
        ex_part2_set = ex_part2_set, 
        tx_ex_part1_set = tx_ex_part1, 
        tx_ex_part2_set = tx_ex_part2, 
        g = g 
      ) 
 
      ref_ex_part_set <- unique(unlist(c(ref_ex_part$common, ref_ex_part$source, ref_ex_part$sink)))
      ref_ex_part_set_ID <- if (length(ref_ex_part_set) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[ref_ex_part_set], collapse=",E")) else ""

      new_row <- data.frame(
        gene = gene,
        ref_ex_part = ref_ex_part_set_ID,
        setdiff1 = setdiff1ID,
        setdiff2 = setdiff2ID,
        transcripts1 = paste(tx1, collapse=", "),
        transcripts2 = paste(tx2, collapse=", "),
        path1 = paste(sapply(vpath1, function(vec) paste(vec, collapse = "-")), collapse=", "),
        path2 = paste(sapply(vpath2, function(vec) paste(vec, collapse = "-")), collapse=", ")
      )
      print(new_row)
      focalexons_df <- rbind(focalexons_df, new_row)
    }
    rm(valid_splits)
    # collapse source to sink on graph
    # collapse source to sink on txmat
    # update bubble_df

  }

  if (nrow(focalexons_df) > 0) {
    focalexons_df <- focalexons_df %>% dplyr::rowwise() %>% dplyr::mutate(adj_cnt = count_items(ref_ex_part)) %>% dplyr::ungroup()
    focalexons_filtered <- focalexons_df %>%
      dplyr::group_by(setdiff1, setdiff2) %>%
      dplyr::filter(adj_cnt == min(adj_cnt)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-adj_cnt)
  } else {
    focalexons_filtered <- focalexons_df
  }

  write.table(focalexons_df, file.path(outdir, "focalexons", paste0(gene, ".focalexons.all.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(focalexons_filtered, file.path(outdir, "focalexons", paste0(gene, ".focalexons.txt")), sep = "\t", quote = FALSE, row.names = FALSE)

}


