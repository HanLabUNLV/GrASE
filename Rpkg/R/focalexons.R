
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

      if (length(setdiff1) == 0 && length(setdiff2) == 0) { next }    # when both paths have no exonfrags ENSG00000067445.21
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




#' @export
get_bubble_topo_interval <- function(bubble, topo_idx) {
  vals = c(topo_idx[bubble$source], topo_idx[bubble$sink])
  vals <- setNames(vals, c("start", "end"))
  return (vals)
}

#' @export
get_bubble_depths <- function(intervals) {
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
  bubble_topo_intervals <- lapply(bubbles_list, grase::get_bubble_topo_interval, topo_idx = topo_idx)
  bubble_depths <- get_bubble_depths(bubble_topo_intervals)
  ordered_bubbles <- bubbles_list[order(-bubble_depths)]
  do.call(rbind.data.frame, ordered_bubbles)

}


# maybe should order by the distance in nodes rather then depth.
#' @export
bubble_ordering2 <- function(g, bubbles_df) {
# suppose `bubble_topo_intervals` is a list of named numeric vectors
#   bubble_topo_intervals[[i]] = c(start=<topo_idx>, end=<topo_idx>)

  topo <- igraph::topo_sort(g, mode = "out")
  topo_idx <- setNames(seq_along(topo), names(topo))

  bubbles_list = split(as.data.frame(bubbles_df), seq(nrow(bubbles_df)))
  bubble_topo_intervals <- lapply(bubbles_list, get_bubble_topo_interval, topo_idx = topo_idx)

  # 1. build a data.frame
  df <- do.call(rbind, lapply(seq_along(bubble_topo_intervals), function(i) {
    iv <- bubble_topo_intervals[[i]]
    data.frame(i = i, start = iv["start"], end = iv["end"])
  }))

  # 2. sort by start increasing, end decreasing
  df <- df[order(df$start, -df$end), ]

  # 3. scan with stack
  stack  <- integer(0)
  parent <- integer(nrow(df))
  depth  <- integer(nrow(df))

  for (k in seq_len(nrow(df))) {
    # Pop any intervals that have already ended before this one starts
    while (length(stack) > 0 && df$end[stack[length(stack)]] < df$end[k]) {
      stack <- stack[-length(stack)]
    }
    # The current top of the stack is the parent (if any)
    parent[k] <- if (length(stack)>0) stack[length(stack)] else NA
    # Depth is how many intervals currently open
    depth[k]  <- length(stack)
    # Push this interval onto the stack
    stack     <- c(stack, k)
  }

  # 4. retrieve order from most inner → outer:
  inner_to_outer <- order(-depth)    # highest depth first
  ordered_bubbles_df <- bubbles_df[inner_to_outer, ]
  ordered_bubbles_df

}


#' @export
bubble_ordering3 <- function(g, bubbles_df) {

  topo <- igraph::topo_sort(g, mode = "out")
  topo_idx <- setNames(seq_along(topo), names(topo))

  bubbles_list = split(as.data.frame(bubbles_df), seq(nrow(bubbles_df)))
  bubble_topo_intervals <- lapply(bubbles_list, get_bubble_topo_interval, topo_idx = topo_idx)


# build df of intervals
  df <- do.call(rbind, lapply(seq_along(bubble_topo_intervals), function(i) {
        iv <- bubble_topo_intervals[[i]]
        data.frame(i      = i,
          topo_start  = iv["start"],
          topo_end    = iv["end"],
          stringsAsFactors = FALSE)
        }))
  rownames(df) = NULL

  bubbles_df = cbind.data.frame(bubbles_df, df)
# compute length
  bubbles_df$length <- igraph::vertex_attr(g, index=bubbles_df$sink)$sg_id - igraph::vertex_attr(g, index=bubbles_df$source)$sg_id

# sort by start ↑, end ↓ (for nesting scan)
  bubbles_df <- bubbles_df[order(bubbles_df$topo_start, -bubbles_df$topo_end), ]

# stack‐scan to get parent & depth
  stack  <- integer(0)
  bubbles_df$parent <- NA_integer_
  bubbles_df$depth  <- 0L

  for (k in seq_len(nrow(bubbles_df))) {
    while (length(stack) > 0 && bubbles_df$topo_end[stack[length(stack)]] < bubbles_df$topo_end[k]) {
      stack <- stack[-length(stack)]
    }
    bubbles_df$parent[k] <- if (length(stack)>0) stack[length(stack)] else NA
      bubbles_df$depth[k]  <- length(stack)
      stack        <- c(stack, k)
  }

# final ordering: length ↓, then depth ↓
  order_idx      <- order(bubbles_df$length, -bubbles_df$depth)
  ordered_bubbles <- bubbles_df[order_idx, ]
  ordered_bubbles
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




#' @export
valid_partitions <- function(dag) {
  nodes <- igraph::V(dag)$name
  n     <- length(nodes)

  # reverse the graph so that "downward-closed" in dag becomes "lower set" in dag_rev
  dag_rev <- igraph::reverse_edges(dag)

  # get a topological ordering of dag_rev
  topo <- igraph::topo_sort(dag_rev, mode="out")$name

  # for each node, precompute its direct predecessors in dag_rev
  preds <- lapply(topo, function(v) {
    igraph::neighbors(dag_rev, v, mode="in")$name
  })
  names(preds) <- topo

  valid_splits <- vector("list", 0)
  split_idx    <- 1L

  # recursively build all lower‐sets (order ideals) in dag_rev
  recurse <- function(i, current_set) {
    if (i > n) {
      # record only non‐empty, non‐full sets
      k <- length(current_set)
      if (k > 0 && k < n) {
        # inside your leaf handler
        valid_splits[[split_idx]] <<- list(
          group1 = sort(current_set),
          group2 = sort(setdiff(nodes, current_set))
        )
        split_idx <<- split_idx + 1L
      }
      return()
    }

    v <- topo[i]
    # if all predecessors of v are already in current_set, we *can* include v
    if (all(preds[[v]] %in% current_set)) {
      recurse(i + 1L, c(current_set, v))
    }
    # always also branch on *not* including v
    recurse(i + 1L, current_set)
  }

  recurse(1L, character(0))
  valid_splits
}




#' @export
valid_partitions_idx <- function(dag) {
  # 1) Map node names → 1:n integer indices
  nodes   <- igraph::V(dag)$name
  n       <- length(nodes)
  idx_map <- setNames(seq_len(n), nodes)

  # 2) Reverse & topo-sort once
  dag_rev  <- igraph::reverse_edges(dag)
  topo_vec <- igraph::topo_sort(dag_rev, mode = "out")$name

  # 3) Precompute predecessors as integer vectors
  preds <- lapply(topo_vec, function(v) {
    idx_map[ igraph::neighbors(dag_rev, v, mode = "in")$name ]
  })
  names(preds) <- topo_vec

  # 4) Prepare a growing list buffer
  capacity <- 256L
  out      <- vector("list", capacity)
  out_i    <- 0L

  # 5) Recursive builder (on integer sets)
  recurse <- function(i, cur) {
    if (i > n) {
      k <- length(cur)
      if (k > 0L && k < n) {
        out_i <<- out_i + 1L
        # grow buffer if needed
        if (out_i > capacity) {
          capacity <<- capacity * 2L
          length(out) <<- capacity
        }
        out[[out_i]] <<- cur
      }
      return()
    }
    v_name <- topo_vec[i]
    v_idx  <- idx_map[v_name]

    # include v if all its preds are already in cur
    if (all(preds[[v_name]] %in% cur)) {
      recurse(i + 1L, c(cur, v_idx))
    }
    # also the branch that skips v
    recurse(i + 1L, cur)
  }

  recurse(1L, integer(0))
  # only return the filled portion
  out[seq_len(out_i)]
}




#' find valid partitions that represents the hierarchical relationship between paths 
#' @export
valid_partitions_slow <- function(g, max_powerset) {

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


intersect_stream <- function(lst) {
  if (length(lst) == 0L) return(vector(mode = "integer", length = 0))
  # sort by length so we intersect the smallest sets first
  ord <- order(vapply(lst, length, integer(1)))
  curr <- unique(lst[[ord[1]]])
  for (i in ord[-1]) {
    # filter in place
    curr <- curr[curr %in% lst[[i]]]
    if (length(curr) == 0L) break
  }
  curr
}

#' @export
find_focal_and_ref_exparts_for_split <- function(g, source, sink, split, parsed_partitions, parsed_paths, tx_ex_parts, focalexons_df, gene=gene) {

      group1 = as.numeric(split$group1)
      group2 = as.numeric(split$group2)
      tx1 = unique(unlist( parsed_partitions[group1]))
      tx2 = unique(unlist( parsed_partitions[group2]))
      print (tx1)
      print (tx2)

      ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
      g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part']) # have to delete all edges in one command
      g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex']) # have to delete all edges in one command
      ex_or_in_gexon <- igraph::edge_attr(g_exon, "ex_or_in")
      names(ex_or_in_gexon) <- igraph::edge_attr(g_exon, "name")
      ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
      names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")
      dexseq_frag <- igraph::edge_attr(g_expart, "dexseq_fragment")
      names(dexseq_frag) <- igraph::edge_attr(g_expart, "name")

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

      setdiff1ID <- if (length(setdiff1) > 0) paste0("E", paste(dexseq_frag[setdiff1], collapse=",E")) else ""
      setdiff2ID <- if (length(setdiff2) > 0) paste0("E", paste(dexseq_frag[setdiff2], collapse=",E")) else ""

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
      ref_ex_part_set_ID <- if (length(ref_ex_part_set) > 0) paste0("E", paste(dexseq_frag[ref_ex_part_set], collapse=",E")) else ""

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
      print(new_row)
      focalexons_df <- rbind(focalexons_df, new_row)
  return(list(focalexons_df = focalexons_df))
}






#' @export
edge_in_txpath <- function(g, trans, e){
   sapply(igraph::edge_attr(g)[trans], `[`, e)
}


#' @export
find_tx_with_epath <- function(g, trans, edges_list) {
  eid = as.integer(igraph::E(g))
  names(eid) = igraph::edge_attr(g, "name")
  tx_to_update = list()
  for ( edges in edges_list) { 
    epath_name = paste0(edges, collapse = ",")
    tx_to_update[[epath_name]] <- list()
    for (e in edges) {
      print(e)
      tx_to_update[[epath_name]] <- append(tx_to_update[[epath_name]], list(trans[grase::edge_in_txpath(g, trans, eid[e])]))
    }
    names(tx_to_update[[epath_name]]) = edges
  }
  tx_to_update
}





#' @export
update_txmat_after_bubble_collapse <- function(g, tx_to_update, subnamelist, txmat) {

  txnames = names(tx_to_update)
  print("before:")
  print(txmat[txnames,subnamelist])
  for (tx_name in txnames) {
    edge_list = tx_to_update[[tx_name]] 
    for (e in edge_list) {
      i = which(subnamelist == e[1]) 
      j = which(subnamelist == e[2])
      txmat[tx_name,as.character(as.numeric(subnamelist[i]):as.numeric(subnamelist[j]))] = TRUE
    }
  }
  print("after:")
  print(txmat[txnames,subnamelist])
  txmat 
}


#' @export
update_txmat_after_bubble_collapse2 <- function(g, tx_to_update, representative_tx, source, sink, txmat) {
  start = which(colnames(txmat)==source)
  end = which(colnames(txmat)==sink)
  nodes = colnames(txmat)
  replacing_pattern = txmat[representative_tx[1], nodes[start:end]] 
  for (tx_list in tx_to_update) {
    for (tx in tx_list) {
      print("before")
      print(txmat[tx, nodes[start:end]])
      txmat[tx, nodes[start:end]] = replacing_pattern
      print("after")
      print(txmat[tx, nodes[start:end]])
    }
  }
  txmat 
}



#' @export
update_txpaths_after_bubble_collapse <- function(g, tx_to_update, epath) {

  updated_edges = c()
  for (x in 1:(length(subnamelist)-1)) {
    e = igraph::E(g)[ .from(subnamelist[x]) & .to(subnamelist[x+1]) & (igraph::E(g)$ex_or_in != 'ex_part') ]
    updated_edges = c(updated_edges, e)
    for (tx in rownames(txmat)) {
      igraph::edge_attr(g, tx, index = e) <- (txmat[tx,subnamelist[x]] && txmat[tx,subnamelist[x+1]])
    }
  }
  txnames = names(tx_to_update)
  print ("before")
  print (igraph::edge_attr(g, index=updated_edges)[txnames])
  for (tx_name in txnames) {
    edge_list = tx_to_update[[tx_name]] 
    for (e in edge_list) {
      i = which(subnamelist == e[1]) 
      j = which(subnamelist == e[2])
      for (x in i:(j-1)) {
        e = igraph::E(g)[ .from(subnamelist[x]) & .to(subnamelist[x+1]) & (igraph::E(g)$ex_or_in != 'ex_part') ]
        igraph::edge_attr(g, tx_name, index = e) <- TRUE
      }
    }
  }
  print ("after")
  print (igraph::edge_attr(g, index=updated_edges)[txnames])
  return (g)
}




#' @export
update_txpaths_after_bubble_collapse2 <- function(g, tx_list, epath) {
  eid = as.integer(igraph::E(g))
  names(eid) = igraph::edge_attr(g, "name")
  for (idx in 1:length(tx_list)) {
    tx_to_update = tx_list[[idx]][[1]]
    for (e in epath) {
      for (tx in tx_to_update) {
        igraph::edge_attr(g, tx, index = eid[e]) <- TRUE
      }
    }
  }
  return (g)
}







#' Compute focal exons from graph and splicing structure
#' @export
focal_exons_gene_powerset <- function(gene, g, sg, outdir, max_path = 30, collapse_bubbles=TRUE) {


#options(keep.source = TRUE)
#Rprof("mem.line.out", line.profiling = TRUE, memory.profiling = TRUE)

  focalexons_df <- data.frame()
  g <- grase::set_edge_names(g)

  txpaths <- grase::txpath_from_edgeattr(g)
  #txmat <- grase::make_matrix_from_txpath_igraph(g, txpaths)
  #txmat_orig = txmat 
  #trans <- rownames(txmat)
  trans <- grase::transcripts_from_igraph(g)
  bubbles_df <- grase::detect_bubbles_igraph(g)     # high_mem
  if ( nrow(bubbles_df) == 0) {
    print("single transcript, no bubbles")
    return (NULL) 
  }

  #bubbles_ordered = grase::bubble_ordering(g, bubbles_df)
  bubbles_ordered = grase::bubble_ordering3(g, bubbles_df)
  bubbles_orig = bubbles_ordered
  bubbles_ordered = bubbles_ordered[,1:ncol(bubbles_df)]


  for (bubble_idx in 1:nrow(bubbles_ordered)) {
  
    ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
    g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex_part']) # have to delete all edges in one command
    g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == 'ex']) # have to delete all edges in one command

    ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
    names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")

    if (bubble_idx > nrow(bubbles_ordered)) {
        print ("collapsed all bubbles. exiting the loop 1")
        break
    }
    source <- bubbles_ordered$source[[bubble_idx]]
    sink <- bubbles_ordered$sink[[bubble_idx]]
    if (source == "R" && sink == "L") next
    print (paste("source: ", source, " sink: ", sink))

    #txpaths <- lapply(SplicingGraphs::txpath(sg), as.character)
    txpaths <- grase::txpath_from_edgeattr(g)
    #tx_exonpaths = lapply(txpaths, grase::from_vpath_to_exon_path, g=g)         # high_mem
    #tx_ex_parts = lapply(tx_exonpaths, grase::from_exon_path_to_exonic_part_path, g=g) #high_mem
    tx_exonpaths = lapply(txpaths, grase::from_vpath_to_exon_path_simple, g=g_exon)         # high_mem
    tx_ex_parts = lapply(tx_exonpaths, grase::from_epath_to_expart_path_simple, g=g_expart) #high_mem

    parsed_paths <- lapply(bubbles_ordered$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_ordered$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])

    vpaths = lapply(parsed_paths, function(vec) c(source, vec, sink))
    epaths = lapply(vpaths, grase::from_vpath_to_exon_path_simple, g=g_exon) 

    ex_part_paths <- lapply(epaths, grase::from_epath_to_expart_path_simple, g=g_expart)
    #ex_part_paths <- lapply(ex_part_paths, function(x) {return (x[igraph::edge_attr(g)$ex_or_in[x] == 'ex_part'])})
    ex_part_paths <- lapply(ex_part_paths, function(x) {return (x[ex_or_in_gexpart[x] == 'ex_part'])})

    #names(ex_part_paths)  = sapply(parsed_partitions, function(vec) paste(vec, collapse=","))
    names(ex_part_paths)  = 1:length(ex_part_paths) 

    print (paste("ex_part_paths len:", length(ex_part_paths)))
    contain_dag <- grase::path_subset_relation(ex_part_paths) 
    print (paste("contain_dag len:", length(igraph::E(contain_dag))))
    if (length(ex_part_paths) > max_path) {
        print("bubble with more than 30 alternative paths: skipping")
        next;
    }
    #valid_splits <- grase::valid_partitions(contain_dag)
    #valid_splits <- grase::remove_symmetric_splits(valid_splits)

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
      retval = grase::find_focal_and_ref_exparts_for_split(g=g, source=source, sink=sink,  split=split, parsed_partitions=parsed_partitions, parsed_paths=parsed_paths, tx_ex_parts=tx_ex_parts, focalexons_df=focalexons_df, gene=gene) 
      focalexons_df = retval$focalexons_df
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
        print("no edges to remove: next")
        next;
      }
      tx_to_update = txs_with_epath[idx_to_remove] 
      edges_to_remove = unique(unlist(epaths[idx_to_remove]))
      g <- igraph::delete_edges(g, edges_to_remove) # have to delete all edges in one command
      # never delete edges by index in a loop. the index is not what you think after iteration. 
      #for (edges in epaths[idx_to_remove]) {
      #  g <- igraph::delete_edges(g, edges)
      #}
      epaths = lapply(vpaths, grase::from_vpath_to_exon_path, g=g) 
      epaths_to_keep = epaths[idx_to_keep]
      g <- grase::update_txpaths_after_bubble_collapse2(g, tx_to_update, epaths_to_keep) 

      # update bubbles_df
      g <- grase::set_txpath_to_vertex_attr(g)
      bubbles_updated <- grase::detect_bubbles_igraph(g)            # high_mem
      if (nrow(bubbles_updated) == 0) {
        print ("collapsed all bubbles. exiting the loop 2")
        break
      }
      bubbles_updated_ordered = grase::bubble_ordering3(g, bubbles_updated)     
      finished_sources = bubbles_ordered$source[(1:bubble_idx)]
      finished_sinks = bubbles_ordered$sink[(1:bubble_idx)]
      finished_bubbles_idx = (bubbles_updated_ordered$source %in% finished_sources) & (bubbles_updated_ordered$sink %in% finished_sinks)
      bubbles_updated_ordered = bubbles_updated_ordered[!finished_bubbles_idx,]
      bubbles_ordered = rbind.data.frame(bubbles_ordered[1:bubble_idx,1:ncol(bubbles_df)],bubbles_updated_ordered[,1:ncol(bubbles_df)])
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

#Rprof(NULL)
  focalexons_filtered
}


