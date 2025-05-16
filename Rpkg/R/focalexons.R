
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
        g = g 
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
  half_powerset <-  function(x) {
    unlist(lapply(1:floor(length(x)/2), function(k) combn(x, k, simplify = FALSE)), recursive = FALSE)
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


#' @export
remove_symmetric_splits <- function(splits){
  # function to produce a symmetric-invariant key
  canonical_key <- function(split) {
    g1 <- paste(sort(strsplit(split$group1, ",")[[1]]), collapse = ",")
    g2 <- paste(sort(strsplit(split$group2, ",")[[1]]), collapse = ",")
    # sort the two side‑strings so that key(A|B)=key(B|A)
    paste(sort(c(g1, g2)), collapse = "|")
  }

  # compute keys
  keys <- vapply(splits, canonical_key, "")
  # keep only the first occurrence of each key
  unique_splits <- splits[!duplicated(keys)]
  unique_splits
}

#' @export
find_focal_and_ref_exparts_for_split <- function(g, source, sink, a_split, a_parsed_partitions, a_parsed_paths, a_tx_ex_parts, a_exonpath_to_remove, a_expartpath_to_replace, a_focalexons_df) {

      group1 = as.numeric(a_split$group1)
      group2 = as.numeric(a_split$group2)
      tx1 = unique(unlist( a_parsed_partitions[group1]))
      tx2 = unique(unlist( a_parsed_partitions[group2]))
      print (tx1)
      print (tx2)

      vpath1 = lapply(a_parsed_paths[group1], function(vec) c(source, vec, sink)) 
      vpath2 = lapply(a_parsed_paths[group2], function(vec) c(source, vec, sink)) 
   
      epath1 = lapply(vpath1, from_vpath_to_exon_path, g=g) 
      epath2 = lapply(vpath2, from_vpath_to_exon_path, g=g) 
      a_exonpath_to_remove = c(a_exonpath_to_remove, epath1)
      a_exonpath_to_remove = c(a_exonpath_to_remove, epath2)
      print("here1")
     
      ex_part1 <- lapply(epath1, from_exon_path_to_exonic_part_path, g=g) 
      a_expartpath_to_replace <- c(a_expartpath_to_replace, ex_part1)
      ex_part1 <- lapply(ex_part1, function(x) {return (x[igraph::edge_attr(g)$ex_or_in[x] == 'ex_part'])})
      ex_part1_set <- Reduce(intersect, ex_part1)
      ex_part2 <- lapply(epath2, from_exon_path_to_exonic_part_path, g=g) 
      a_expartpath_to_replace <- c(a_expartpath_to_replace, ex_part2)
      ex_part2 <- lapply(ex_part2, function(x) {return (x[igraph::edge_attr(g)$ex_or_in[x] == 'ex_part'])})
      ex_part2_set <- Reduce(intersect, ex_part2)
      print("here2")
     
  
      setdiff1 <- setdiff(ex_part1_set, ex_part2_set)
      setdiff2 <- setdiff(ex_part2_set, ex_part1_set)

      setdiff1ID <- if (length(setdiff1) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[setdiff1], collapse=",E")) else ""
      setdiff2ID <- if (length(setdiff2) > 0) paste0("E", paste(igraph::edge_attr(g)$dexseq_fragment[setdiff2], collapse=",E")) else ""

      tx_ex_part1 <- a_tx_ex_parts[tx1[1]]
      tx_ex_part1 <- unique(unlist(tx_ex_part1))
      tx_ex_part2 <- a_tx_ex_parts[tx2[1]] 
      tx_ex_part2 <- unique(unlist(tx_ex_part2))
      print("here3")

      ref_ex_part <- find_reference_exonic_part(
        source = source, 
        sink = sink, 
        ex_part1_set = ex_part1_set, 
        ex_part2_set = ex_part2_set, 
        tx_ex_part1_set = tx_ex_part1, 
        tx_ex_part2_set = tx_ex_part2, 
        g = g 
      ) 
      print("here4")
 
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
      a_focalexons_df <- rbind(a_focalexons_df, new_row)
  return(list(exonpath_to_remove = a_exonpath_to_remove, expartpath_to_replace = a_expartpath_to_replace, focalexons_df = a_focalexons_df))
}



#' @export
find_edges_to_replace_after_bubble_collapse <- function(g, source, sink) {  

    edges_to_replace = list()
    sg_id_list = igraph::vertex_attr(g)$sg_id
    idx_start <- which(sg_id_list == igraph::V(g)[source]$sg_id )
    idx_end <- which(sg_id_list == igraph::V(g)[sink]$sg_id )
    sublist = sg_id_list[idx_start:idx_end]
    subnamelist = igraph::V(g)[igraph::vertex_attr(g, 'sg_id') %in% sublist]$name
    # now go node by node and find the paths.
    # if no path along the way replace with new edge
    # this only works because the nodes are ordered by their genomic positions
    for (x in 1:(length(sublist)-1)) {
      sg_start = sublist[x]
      sg_end = sublist[x+1]
      current_start <- igraph::V(g)[sg_id == sg_start]
      current_end <- igraph::V(g)[sg_id == sg_end]
      remaining_edges <- igraph::E(g)[current_start %--% current_end]
      if (length(remaining_edges) == 0) {
        edgetype = NA
        if ((current_start$name == 'R') || (current_end$name == 'R')) {
          edgetype = 'R'
        } 
        else if ((current_start$name == 'L') || (current_end$name == 'L')) {
          edgetype = 'L'
        }
        else {
          edgetype = 'in'
        }
        edges_to_replace = c(edges_to_replace, list(list(start=current_start$name, end=current_end$name, type = edgetype, add=TRUE, txupdate = TRUE)))
      }
      else if (all(igraph::edge_attr(g, 'ex_or_in', index=remaining_edges) == 'ex_part')) {
        edges_to_replace = c(edges_to_replace, list(list(start=current_start$name, end=current_end$name, type = 'ex', add=TRUE, txupdate = TRUE)))
      }
      else {
        print ("comes here when there is both exon and expart remaining")
        print ("we don't have to add since there is an exon that covers the range")
        print(remaining_edges)
        #print(igraph::edge_attr(g, index=remaining_edges))
        edges_to_replace = c(edges_to_replace, list(list(start=current_start$name, end=current_end$name, type = 'ex', add=FALSE, txupdate = TRUE)))
      }
    }
  return (list(edges_to_replace=edges_to_replace, subnamelist=subnamelist))
}


#' @export
replace_edges_after_bubble_collapse <- function(g, edges_to_replace) {

    # add missing edges 
    for (e in edges_to_replace) {
      print (paste("replacing edges ", e$start,"-", e$end, ":", e$type)) 
      if (e$add) {
        g <- igraph::add_edges(g, c(e$start, e$end), attr = c(list('ex_or_in' = e$type)))
      } 
    }
  return (g)
}


#' @export
update_txpaths_after_bubble_collapse <- function(g, tx_to_update, subnamelist, txmat) {
  for (x in 1:(length(subnamelist)-1)) {
    e = igraph::E(g)[ .from(subnamelist[x]) & .to(subnamelist[x+1]) & (igraph::E(g)$ex_or_in != 'ex_part') ]
    for (tx in rownames(txmat)) {
      igraph::edge_attr(g, tx, index = e) <- (txmat[tx,subnamelist[x]] && txmat[tx,subnamelist[x+1]])
      if (tx %in% tx_to_update) {
        igraph::edge_attr(g, tx, index = e) <- TRUE
      }
    }
  } 
  return (g)
}


#' @export
edge_in_txpath <- function(g, trans, e){
   sapply(igraph::edge_attr(g)[trans], `[`, e)
}


#' @export
find_tx_to_update <- function(g, trans, edges_to_remove) {
  tx_to_update = c()
  for ( e in edges_to_remove) {
    tx_to_update = c(tx_to_update, trans[edge_in_txpath(g, trans, e)])
  }
  unique(tx_to_update)
}


#' Compute focal exons from graph and splicing structure
#' @export
focal_exons_gene_powerset <- function(gene, g, sg, outdir, max_powerset = 10000, collapse_bubbles=TRUE) {


  focalexons_df <- data.frame()

  g_orig <- g
  txpaths <- grase::txpath_from_edgeattr(g)
  txpaths_orig <- txpaths
  txmat <- grase::make_matrix_from_txpath_igraph(txpaths)
  txmat_orig = txmat 
  trans <- rownames(txmat)
  bubbles_df <- grase::detect_bubbles_from_mat(txmat)
  bubbles_ordered = grase::bubble_ordering(g, bubbles_df)
  bubbles_orig = bubbles_ordered

  if ( nrow(bubbles_df) == 0) {
    print("single transcript, no bubbles")
    return (NULL) 
  }

  for (bubble_idx in seq_along(bubbles_ordered$partitions)) {
    source <- bubbles_ordered$source[[bubble_idx]]
    sink <- bubbles_ordered$sink[[bubble_idx]]
    if (source == "R" && sink == "L") next

    #txpaths <- lapply(SplicingGraphs::txpath(sg), as.character)
    txpaths <- grase::txpath_from_edgeattr(g)
    tx_exonpaths = lapply(txpaths, grase::from_vpath_to_exon_path, g=g)
    tx_ex_parts = lapply(tx_exonpaths, grase::from_exon_path_to_exonic_part_path, g=g)

    parsed_paths <- lapply(bubbles_ordered$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_ordered$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    n_partitions = length(parsed_partitions)
    vpaths = lapply(parsed_paths, function(vec) c(source, vec, sink))
    epaths = lapply(vpaths, grase::from_vpath_to_exon_path, g=g) 
    ex_part_paths <- lapply(epaths, grase::from_exon_path_to_exonic_part_path, g=g)
    ex_index = igraph::edge_attr(g)$ex_or_in == 'ex_part'
    ex_part_paths <- lapply(ex_part_paths, function(x) {return (x[igraph::edge_attr(g)$ex_or_in[x] == 'ex_part'])})
    names(ex_part_paths)  = sapply(parsed_partitions, function(vec) paste(vec, collapse=","))
    names(ex_part_paths)  = 1:length(ex_part_paths) 

    print (paste("ex_part_paths len:", length(ex_part_paths)))
    contain_dag <- grase::path_subset_relation(ex_part_paths) 
    print (paste("contain_dag len:", length(igraph::E(g))))
    valid_splits <- grase::valid_partitions(contain_dag, max_powerset)
    valid_splits <- grase::remove_symmetric_splits(valid_splits)
    rm(contain_dag)
    
    exonpath_to_remove <- integer(0)
    expartpath_to_replace <- integer(0)
    for (split in valid_splits) {

      retval = grase::find_focal_and_ref_exparts_for_split(g=g, source=source, sink=sink,  a_split=split, a_parsed_partitions=parsed_partitions, a_parsed_paths=parsed_paths, a_tx_ex_parts=tx_ex_parts, a_exonpath_to_remove=exonpath_to_remove, a_expartpath_to_replace=expartpath_to_replace, a_focalexons_df=focalexons_df) 
      exonpath_to_remove = retval$exonpath_to_remove 
      expartpath_to_replace = retval$expartpath_to_replace 
      focalexons_df = retval$focalexons_df
    }

    if (collapse_bubbles) {
      # collapse source to sink on graph
      # remove edges that are part of the bubble and replace with exon edges
      edges_to_remove = unique(unlist(exonpath_to_remove))
      print('edges_to_remove')
      print(igraph::E(g)[edges_to_remove])
      tx_to_update <- grase::find_tx_to_update(g, trans, edges_to_remove)
      g <- igraph::delete_edges(g, edges_to_remove)
      # replace missing edges
      retval = grase::find_edges_to_replace_after_bubble_collapse(g, source, sink)
      edges_to_replace = retval$edges_to_replace
      subnamelist = retval$subnamelist
      #print('edges_to_replace')
      #print(matrix(unlist(edges_to_replace), ncol = 5, byrow = TRUE))
      g <- grase::replace_edges_after_bubble_collapse(g, edges_to_replace) 
      g <- grase::update_txpaths_after_bubble_collapse(g, tx_to_update, subnamelist, txmat) 

      # collapse source to sink on txmat
      txmat[tx_to_update, subnamelist] = TRUE
      # update bubbles_df
      bubbles_updated <- grase::detect_bubbles_from_mat(txmat)
      if (nrow(bubbles_updated) == 0) {
        print ("collapsed all bubbles. exiting the loop")
        break
      }
      bubbles_updated_ordered = grase::bubble_ordering(g, bubbles_updated)     
      bubbles_ordered = rbind.data.frame(bubbles_ordered[(1:bubble_idx),],bubbles_updated_ordered)
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


