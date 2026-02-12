# Common utility functions for graph analysis and bubble processing


.onLoad <- function(libname, pkgname) {
  op <- options()
  op.my_package <- list(
    my_package.debug = FALSE
  )
  toset <- !(names(op.my_package) %in% names(op))
  if (any(toset)) options(op.my_package[toset])
}

log_debug <- function(msg) {
  if (getOption("my_package.debug", default = FALSE)) {
    message("[DEBUG] ", msg)
  }
}

#' Count items in comma-separated string
#' @export
count_items <- function(x) {
  if (x == "") return(0)
  length(stringr::str_split(x, ",")[[1]])
}

#' Intersect all elements in a list efficiently
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

#' Check if set a is a proper subset of b
is_proper_subset <- function(a, b) {
  length(a) < length(b) && all(a %in% b)
}

#' Create subset relation graph from sets
#' @export
path_subset_relation <- function(sets) {
  n <- length(sets)
  if (n == 0) return(igraph::graph.empty())
  
  # Create adjacency list
  edges <- c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (is_proper_subset(sets[[i]], sets[[j]])) {
        edges <- c(edges, i, j)
      } else if (is_proper_subset(sets[[j]], sets[[i]])) {
        edges <- c(edges, j, i)
      }
    }
  }
  
  g <- igraph::graph(edges, n = n, directed = TRUE)
  igraph::V(g)$name <- names(sets) %||% as.character(1:n)
  g
}

#' Generate valid binary partitions from DAG
#' @export
valid_partitions <- function(dag) {
  nodes <- igraph::V(dag)$name
  n     <- length(nodes)

  # get a topological ordering of dag
  topo <- igraph::topo_sort(dag, mode = "out")$name

  # for each node, precompute its direct predecessors in dag
  preds <- lapply(topo, function(v) {
    igraph::neighbors(dag, v, mode = "in")$name
  })
  names(preds) <- topo

  valid_splits <- vector("list", 0)
  split_idx    <- 1L

  # recursively build for downward-closed sets (lower sets) in dag
  recurse <- function(i, current_set) {
    if (i > n) {
      k <- length(current_set)
      if (k > 0 && k < n) {
        valid_splits[[split_idx]] <<- list(
          group1 = sort(current_set),
          group2 = sort(setdiff(nodes, current_set))
        )
        split_idx <<- split_idx + 1L
      }
      return()
    }

    v <- topo[i]
    # Include v only if all its subset-predecessors are already in the set
    if (all(preds[[v]] %in% current_set)) {
      recurse(i + 1L, c(current_set, v))
    }
    # Also consider skipping v
    recurse(i + 1L, current_set)
  }

  recurse(1L, character(0))
  valid_splits
}

#' Generate valid partition indices from DAG
#' @export
valid_partitions_idx <- function(dag) {
  nodes <- igraph::V(dag)$name
  n     <- length(nodes)

  # get a topological ordering of dag
  topo <- igraph::topo_sort(dag, mode = "out")$name

  # for each node, precompute its direct predecessors in dag
  preds <- lapply(topo, function(v) {
    igraph::neighbors(dag, v, mode = "in")$name
  })
  names(preds) <- topo

  valid_splits <- vector("list", 0)
  split_idx    <- 1L

  # recursively build for downward-closed sets (lower sets) in dag
  recurse <- function(i, current_set) {
    if (i > n) {
      k <- length(current_set)
      if (k > 0 && k < n) {
        current_idx <- which(nodes %in% current_set)
        valid_splits[[split_idx]] <<- current_idx
        split_idx <<- split_idx + 1L
      }
      return()
    }

    v <- topo[i]
    # Include v only if all its subset-predecessors are already in the set
    if (all(preds[[v]] %in% current_set)) {
      recurse(i + 1L, c(current_set, v))
    }
    # Also consider skipping v
    recurse(i + 1L, current_set)
  }

  recurse(1L, character(0))
  valid_splits
}

#' Generate valid partitions with powerset approach (slower)
#' @export
valid_partitions_slow <- function(g, max_powerset) {
  num_nodes <- igraph::vcount(g)
  if (num_nodes > max_powerset) {
    message(paste("too many nodes (", num_nodes, ") for powerset. returning NULL"))
    return(NULL)
  }
  
  nodes <- igraph::V(g)$name
  all_subsets <- gtools::combinations(num_nodes, rep = FALSE, simplify = FALSE, 
                                     v = 1:num_nodes, set = FALSE)
  
  valid_splits <- list()
  for (i in 1:length(all_subsets)) {
    subset_idx <- all_subsets[[i]]
    if (length(subset_idx) == 0 || length(subset_idx) == num_nodes) next
    
    group1 <- nodes[subset_idx]
    group2 <- nodes[-subset_idx]
    valid_splits[[length(valid_splits) + 1]] <- list(group1 = group1, group2 = group2)
  }
  
  valid_splits
}

#' Remove symmetric splits (where group1/group2 are swapped)
#' @export
remove_symmetric_splits <- function(splits){
  unique_splits <- list()
  for (split in splits) {
    # Canonicalize by ensuring group1 has smaller first element
    if (length(split$group1) > 0 && length(split$group2) > 0) {
      if (split$group1[1] > split$group2[1]) {
        split <- list(group1 = split$group2, group2 = split$group1)
      }
    }
    
    # Check if we've seen this split before
    is_duplicate <- any(sapply(unique_splits, function(existing) {
      identical(existing$group1, split$group1) && identical(existing$group2, split$group2)
    }))
    
    if (!is_duplicate) {
      unique_splits[[length(unique_splits) + 1]] <- split
    }
  }
  unique_splits
}

#' Get topological interval for bubble
#' @export
get_bubble_topo_interval <- function(bubble, topo_idx) {
  source_idx <- topo_idx[bubble$source]
  sink_idx <- topo_idx[bubble$sink]
  c(source_idx, sink_idx)
}

#' Get bubble depths from intervals
#' @export
get_bubble_depths <- function(intervals) {
  depths <- integer(length(intervals))
  for (i in seq_along(intervals)) {
    interval <- intervals[[i]]
    depth <- 0
    for (j in seq_along(intervals)) {
      if (i != j) {
        other_interval <- intervals[[j]]
        if (other_interval[1] >= interval[1] && other_interval[2] <= interval[2]) {
          depth <- depth + 1
        }
      }
    }
    depths[i] <- depth
  }
  depths
}

#' Order bubbles by topological depth
#' @export
bubble_ordering <- function(g, bubbles_df) {
  topo <- igraph::topo_sort(g, mode = "out")
  topo_idx <- setNames(seq_along(topo), topo$name)
  
  intervals <- lapply(1:nrow(bubbles_df), function(i) {
    get_bubble_topo_interval(bubbles_df[i,], topo_idx)
  })
  
  depths <- get_bubble_depths(intervals)
  bubbles_df$depth <- depths
  bubbles_df[order(depths, decreasing = TRUE),]
}

#' Alternative bubble ordering method
#' @export
bubble_ordering2 <- function(g, bubbles_df) {
  topo <- igraph::topo_sort(g, mode = "out")
  topo_idx <- setNames(seq_along(topo), topo$name)
  
  bubbles_df$source_idx <- topo_idx[bubbles_df$source]
  bubbles_df$sink_idx <- topo_idx[bubbles_df$sink]
  bubbles_df$span <- bubbles_df$sink_idx - bubbles_df$source_idx
  
  bubbles_df[order(bubbles_df$span, decreasing = TRUE),]
}

#' Third bubble ordering method (most sophisticated)
#' @export
bubble_ordering3 <- function(g, bubbles_df) {
  if (nrow(bubbles_df) == 0) return(bubbles_df)
  
  topo <- igraph::topo_sort(g, mode = "out")
  topo_idx <- setNames(seq_along(topo), topo$name)
  
  # Calculate metrics for each bubble
  bubbles_df$source_idx <- topo_idx[bubbles_df$source]
  bubbles_df$sink_idx <- topo_idx[bubbles_df$sink]
  bubbles_df$span <- bubbles_df$sink_idx - bubbles_df$source_idx
  
  # Count number of paths per bubble
  bubbles_df$num_paths <- sapply(bubbles_df$paths, function(paths_str) {
    length(strsplit(paths_str, "\\},\\{")[[1]])
  })
  
  # Order by span (descending), then by source position (ascending)
  bubbles_df[order(-bubbles_df$span, bubbles_df$source_idx),]
}

#' Order bubbles innermost-first for correct collapsing
#'
#' Orders bubbles by span ascending (smallest/innermost first), then by
#' source position ascending as tiebreaker. This ensures inner bubbles are
#' processed and collapsed before the outer bubbles that contain them.
#' @param g An igraph splicing graph
#' @param bubbles_df A data frame of detected bubbles
#' @return The reordered bubbles_df with added columns: source_idx, sink_idx, span, num_paths
#' @export
bubble_ordering4 <- function(g, bubbles_df) {
  if (nrow(bubbles_df) == 0) return(bubbles_df)

  topo <- igraph::topo_sort(g, mode = "out")
  topo_idx <- setNames(seq_along(topo), topo$name)

  # Calculate metrics for each bubble
  bubbles_df$source_idx <- topo_idx[bubbles_df$source]
  bubbles_df$sink_idx <- topo_idx[bubbles_df$sink]
  bubbles_df$span <- bubbles_df$sink_idx - bubbles_df$source_idx

  # Count number of paths per bubble
  bubbles_df$num_paths <- sapply(bubbles_df$paths, function(paths_str) {
    length(strsplit(paths_str, "\\},\\{"))
  })

  # Order by span ascending (innermost first), then by source position ascending
  bubbles_df[order(bubbles_df$span, bubbles_df$source_idx),]
}

#' Check if edge is in transcript path
#' @export
edge_in_txpath <- function(g, trans, e){
   sapply(igraph::edge_attr(g)[trans], `[`, e)
}

#' Update transcript paths after collapsing a bubble
#'
#' After removing edges from a collapsed bubble, reassign the affected
#' transcripts to the kept paths by setting edge attributes.
#' @param g An igraph splicing graph
#' @param tx_list List of transcript groups to update (from find_tx_with_epath)
#' @param epath Edge paths to keep (list of edge name vectors)
#' @return The updated igraph graph
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

#' Find transcripts that use given edge paths
#' @export
find_tx_with_epath <- function(g, trans, edges_list) {
  eid = as.integer(igraph::E(g))
  names(eid) = igraph::edge_attr(g, "name")
  tx_to_update = list()
  for ( edges in edges_list) { 
    epath_name = paste0(edges, collapse = ",")
    tx_to_update[[epath_name]] <- list()
    for (e in edges) {
      log_debug(e)
      tx_to_update[[epath_name]] <- append(tx_to_update[[epath_name]], list(trans[grase::edge_in_txpath(g, trans, eid[e])]))
    }
    names(tx_to_update[[epath_name]]) = edges
  }
  tx_to_update
}