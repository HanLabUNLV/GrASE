# R/bubbles.R




#' get transcript names from igraph edge attributes
#' @export
transcripts_from_igraph <- function(g) {
  attrs <- igraph::edge_attr_names(g)
  excluded <- c("sgedge_id", "ex_or_in", "from_pos", "to_pos", "dexseq_fragment")
  trans <- setdiff(attrs, excluded)
}



#' Generate transcript-paths from igraph edge attributes
#' @export
txpath_from_edgeattr <- function(g, type="exon") {

  # Identify transcript attributes on edges (exclude known graph attrs)
  trans <- grase::transcripts_from_igraph(g)
  # Determine root and leaf vertices
  root <- if ("R" %in% igraph::V(g)$name) "R" else stop("No root 'R' vertex found")
  leaf <- if ("L" %in% igraph::V(g)$name) "L" else stop("No leaf 'L' vertex found")

  if (type == "exon") {
    g_tmp = g
    # lose all expart edges
    ex_parts = igraph::E(g_tmp)[igraph::edge_attr(g_tmp)$ex_or_in == 'ex_part']
    g_tmp = igraph::delete_edges(g_tmp, ex_parts)
  }
  else if (type == "ex_part") {
    g_tmp = g
    # lose all expart edges
    ex_parts = igraph::E(g_tmp)[igraph::edge_attr(g_tmp)$ex_or_in == 'ex']
    g_tmp = igraph::delete_edges(g_tmp, ex_parts)
  }
  else {
    print("error: type should be either exon or expart")
    return (NULL)
  }
  txpaths <- setNames(vector("list", length(trans)), trans)
  for (tx in trans) {
    tx_bools <- igraph::edge_attr(g_tmp, tx)
    exon_ids = igraph::E(g_tmp)[tx_bools]
    subg  <- igraph::subgraph_from_edges(g_tmp, exon_ids, delete.vertices = FALSE)
    #print ("subg") 
    #print (igraph::E(subg))
    spath <- igraph::shortest_paths(subg, from = root, to = leaf, mode = "out")$vpath[[1]]
    txpaths[[tx]] <- igraph::V(g_tmp)$name[spath]
    #print(paste("txpath for ", tx))
    #print(txpaths[[tx]])
  }
  txpaths = lapply(txpaths, function(x) { x[2:(length(x)-1)] }) # get rid of first and last (R and L)
  txpaths
}


#' @export
set_txpath_to_vertex_attr <- function(g)
{
  txpaths <- grase::txpath_from_edgeattr(g)
  trans <- names(txpaths)
  for (tx in trans) {
    vlist <- txpaths[[tx]]
    g <- igraph::set_vertex_attr(g, tx, value = rep(FALSE,length(igraph::V(g)))) 
    g <- igraph::set_vertex_attr(g, tx, index = vlist, value=TRUE)
    g <- igraph::set_vertex_attr(g, tx, index = 'R', value=TRUE)
    g <- igraph::set_vertex_attr(g, tx, index = 'L', value=TRUE)
  }
  g 
}


 

#' Build transcript-path logical matrix from index lists
#' @export
make_matrix_from_txpath_igraph <- function(g, txpath_vertex_list) {
   
  attrs <- igraph::vertex_attr_names(g)
  excluded <- c("name", "position", "sg_id", "id")
  trans <- setdiff(attrs, excluded)

  if (length(trans) == 0) {
    print( "transcript vertex attributes not set yet")
    print( "running set_txpath_to_vertex_attr first to set the attributes")
    g <- grase::set_txpath_to_vertex_attr(g) 
 
    attrs <- igraph::vertex_attr_names(g)
    excluded <- c("name", "position", "sg_id", "id")
    trans <- setdiff(attrs, excluded)
  }
  newmat <- matrix(unlist(igraph::vertex_attr(g)[trans]), nrow=length(trans), byrow=TRUE)
  colnames(newmat) <- igraph::V(g)$name
  rownames(newmat) <- trans
  newmat[,"R"] = TRUE
  newmat[,"L"] = TRUE
  newmat
}

#' identifies nodes that are exon starts and ends and marks them with 1 and 2 respectively.  
#' copied from SplicingGraphs
#' @export
get_sgnodetypes <- function (txpathmat, check.matrix = FALSE) 
{
    ans <- integer(ncol(txpathmat))
    names(ans) <- colnames(txpathmat)
    for (i in seq_len(nrow(txpathmat))) {
        idx <- which(txpathmat[i, , drop = FALSE])
        if (length(idx) <= 2L) 
            (next)()
        idx <- idx[-c(1L, length(idx))]
        exon_start_idx <- idx[c(TRUE, FALSE)]
        exon_end_idx <- idx[c(FALSE, TRUE)]
        if (check.matrix) {
            if (any(ans[exon_start_idx] == 2L) || any(ans[exon_end_idx] == 
                1L)) 
                stop("invalid matrix of transcript paths: ", 
                  "some columns in 'txpathmat' seem to correspond ", 
                  "at the same time to an exon start and an exon end")
        }
        ans[exon_start_idx] <- 1L
        ans[exon_end_idx] <- 2L
    }
    ans
}

#' copied from SplicingGraphs
#' @export
is_bubble <- function (txpathmat, i, j) 
{
    txbase <- txpathmat[, i] & txpathmat[, j]   # marks whether i and j are both present (i.e. there is a path from i to j) in each transcript
    if (sum(txbase) <= 1L) 
        return(FALSE)  # only one transcript has the path. no bubble
    for (k in (i + 1L):(j - 1L)) if (all(txpathmat[, k] >= txbase))  # marks whether the route going from i to j has any gaps
        return(FALSE)
    TRUE
}






#' copied from SplicingGraphs
#' @export
get_bubble_variants <- function (txpathmat, sgnodetypes, i, j) 
{
    txbase <- txpathmat[, i] & txpathmat[, j] # which transcripts are passing i and j
    bubble_submat <- txpathmat[txbase, (i + 1L):(j - 1L), drop = FALSE] # columns in between i and j for txbase
    bubble_submat <- bubble_submat[, colSums(bubble_submat) != 0L, drop = FALSE] # only keep nodes that are present in txbase
    bubble_submat_rownames <- rownames(bubble_submat)
    bubble_submat_colnames <- colnames(bubble_submat) # now bubble_submat has only transcript and only nodes that are present between i and j
    ans_path <- IRanges::CharacterList(lapply(seq_len(nrow(bubble_submat)), 
        function(i) bubble_submat_colnames[bubble_submat[i, ]]))
    ans_rpath <- IRanges::IntegerList(lapply(seq_len(nrow(bubble_submat)), 
        function(i) which(bubble_submat[i, ])))
    ans_code <- sapply(seq_len(length(ans_path)), function(k) {
        path <- ans_path[[k]]
        path_len <- length(path)
        if (path_len == 0L) 
            return("0")
        types <- c("-", "^")[sgnodetypes[path]]
        if (length(types) != path_len) 
            stop("some splicing sites in variant path ", "are of type 0")
        if (i == 1L) 
            types[1L] <- "["
        if (j == ncol(txpathmat)) 
            types[length(types)] <- "]"
        paste0(ans_rpath[[k]], types, collapse = "")
    })
    oo <- order(ans_code)
    bubble_submat_rownames <- bubble_submat_rownames[oo]
    ans_path <- ans_path[oo]
    ans_code <- ans_code[oo]
    ans_code1 <- ans_code[-length(ans_code)]
    ans_code2 <- ans_code[-1L]
    is_not_dup <- c(TRUE, ans_code1 != ans_code2)
    ans_path <- ans_path[is_not_dup]
    ans_code <- ans_code[is_not_dup]
    ans_partition <- unname(S4Vectors::splitAsList(bubble_submat_rownames, 
        cumsum(is_not_dup)))
    S4Vectors::DataFrame(partition = ans_partition, path = ans_path, code = ans_code)
}




#' bubble detection between i and j
#' @export
detect_bubbles_i_j <- function(i, j, txpathmat, sgnodetypes)
{
    if (!is_bubble(txpathmat, i, j)) 
        return (NULL)
    bubble_variants <- get_bubble_variants(txpathmat, 
        sgnodetypes, i, j)
    bubble_d <- nrow(bubble_variants)
    if (bubble_d <= 1L) 
        return (NULL)
    ans_source <- colnames(txpathmat)[i]
    ans_sink <- colnames(txpathmat)[j]
    ans_d <- bubble_d
    bubble_partitions <- bubble_variants[, "partition"]
    bubble_partitions <- sapply(bubble_partitions, base::paste, 
        collapse = ",")
    bubble_partitions <- paste0("{", bubble_partitions, 
        "}")
    ans_partitions <- IRanges::CharacterList(bubble_partitions)
    bubble_paths <- bubble_variants[, "path"]
    bubble_paths <- sapply(bubble_paths, base::paste, 
        collapse = ",")
    bubble_paths <- paste0("{", bubble_paths, "}")
    ans_paths <- IRanges::CharacterList(bubble_paths)
    return (list(ans_source = ans_source, ans_sink = ans_sink, ans_d = ans_d, ans_partitions = ans_partitions, ans_paths = ans_paths ))

}

#' Main bubble detection from matrix
#' @export
detect_bubbles_from_mat <- function (txpathmat) 
{
    prev_locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", prev_locale))
    sgnodetypes <- get_sgnodetypes(txpathmat)
    ans_source <- ans_sink <- ans_AScode <- character(0)
    ans_d <- integer(0)
    ans_partitions <- ans_paths <- IRanges::CharacterList()
    ncol <- ncol(txpathmat)
    for (i in 1:(ncol - 2L)) {
        for (j in (i + 2L):ncol) {
          retval = detect_bubbles_i_j(i, j, txpathmat, sgnodetypes)
          if (is.null(retval)) 
            next
          
          ans_source <- c(ans_source, retval$ans_source)
          ans_sink <- c(ans_sink, retval$ans_sink)
          ans_d <- c(ans_d, retval$ans_d)
          ans_partitions <- c(ans_partitions, retval$ans_partitions)
          ans_paths <- c(ans_paths, retval$ans_paths)
        }
    }
    S4Vectors::DataFrame(source = ans_source, sink = ans_sink, d = ans_d, 
        partitions = ans_partitions, paths = ans_paths
        )
}


#' @export
get_sgnodetypes_igraph <- function (g_tmp, check.graph=FALSE) 
{
  nodes = igraph::V(g_tmp)
   
  ans <- integer(length(nodes))
  names(ans) <- names(nodes)
  ex_edges = igraph::E(g_tmp)[igraph::edge_attr(g_tmp, "ex_or_in") == "ex"]

  startends = igraph::ends(g_tmp, ex_edges) 
  starts = startends[,1]
  ends = startends[,2]
  if (check.graph) {
    if (length(intersect(starts, ends)) > 0) {
      stop("invalid matrix of transcript paths: ", 
        "some columns in 'txpathmat' seem to correspond ", 
        "at the same time to an exon start and an exon end")
    }
  } 
  ans[starts] <- 1L
  ans[ends] <- 2L
  ans
}



#' @export
bubble_paths_igraph <- function(g, v_start, v_end) {

  ex_parts = igraph::E(g)[igraph::edge_attr(g)$ex_or_in == 'ex_part']
  if (length(ex_parts) > 0) {
    stop("ex_part edges found on the graph",
      "purge the ex_part edges before running bubble detection")
  }
  
  paths = igraph::all_simple_paths(g, from=v_start, to=v_end)
  paths
}



#' @export
get_bubble_variants_igraph <- function(g, v_start_idx, v_end_idx) {

  # ensure numeric vertex IDs
  if (v_end_idx - v_start_idx < 2) {
    return(list(partition = list(), path = list()))
  }

  # which vertex attrs are "transcripts"?
  all_attr_names <- igraph::vertex_attr_names(g)
  excluded  <- c("name", "position", "sg_id", "id")
  trans     <- setdiff(all_attr_names, excluded)
  all_attrs <- igraph::vertex_attr(g)
  attr_trans <- all_attrs[trans]

  # get the vertex names up front
  v_idx_range <- seq(v_start_idx, v_end_idx)

  # 1) Pull out just the startrow and endrow in one vapply
  two_rows <- vapply(
    attr_trans,
    function(vec) vec[c(v_start_idx, v_end_idx)],
    numeric(2L)
  )
  start_row <- two_rows[1L, ]
  end_row   <- two_rows[2L, ]

  # 2) Which transcripts appear in BOTH start & end?
  base_mask <- (start_row != 0) & (end_row != 0)
  if (!any(base_mask)) {
    return(list(partition = list(), path = list()))
  }
  base_idx   <- which(base_mask)
  base_names <- trans[ base_idx ]
  attr_base  <- all_attrs[ base_names ]

  # 3) Check each internal vertex one at a time
  internal_vertices <- setdiff(seq(v_start_idx, v_end_idx), c(v_start_idx, v_end_idx))
  for (i_vert in internal_vertices) {
    blocked <- TRUE
    for (vec in attr_base) {
     if (vec[i_vert] == 0) {
       blocked <- FALSE
       break
     }
    }
    if (blocked) {
     return(list(partition = list(), path = list()))
    }
  }


  # 1) Identify internal vertices (exclude the two endpoints)
  #internal_vertices <- v_idx_range[-c(1, length(v_idx_range))]
  n_int            <- length(internal_vertices)
  n_base           <- length(base_names)

  # If there are no internal vertices, everything downstream is trivial:
  if (n_int == 0L) {
    # In this case, no internal node can block anything, but also
    # there are no columns to group by.  So partition = empty, path = empty.
    partitions   <- list()
    paths_unique <- list()
    return(list(partition = partitions, path = paths_unique))
  }

  # 2) Build each transcripts bitpattern over all internal vertices, one row at a time
  #    row_strs_raw[ i ] will be a string like "01011" of length = n_int,
  #    indicating which internal_vertices[j] contain base_names[i].
  row_strs_raw <- vapply(
    base_names,
    function(attr_name) {
      # Fetch count/presence of this transcript at ALL internal vertices:
      vals <- as.numeric(
        igraph::vertex_attr(g, attr_name, index = internal_vertices)
      )
      # Convert to bits (0/1) and collapse into a single string for this row:
      paste0(as.integer(vals != 0), collapse = "")
    },
    character(1L)
  )
  names(row_strs_raw) <- base_names

  # 3) Determine which columns (positions 1..n_int) are all zero across every row.
  #    We only want to keep columns j for which at least one row_strs_raw[i][j] == "1".
  keep_col_flags <- vapply(
    seq_len(n_int),
    function(j) {
      any(substr(row_strs_raw, j, j) == "1")
    },
    logical(1L)
  )
  if (!any(keep_col_flags)) {
    # If no column has a 1, then all internal vertices are irrelevant: return empty.
    return(list(partition = list(), path = list()))
  }
  keep_positions <- which(keep_col_flags)

  # 4) Build the final, column pruned row_strs: keep only bits in positions `keep_positions`
  #    row_strs_final[i] is now a shorter string (length = length(keep_positions)).
  row_strs_final <- vapply(
    row_strs_raw,
    function(bitstr) {
      bits <- strsplit(bitstr, "")[[1]]
      paste0(bits[keep_positions], collapse = "")
    },
    character(1L)
  )
  names(row_strs_final) <- names(row_strs_raw)  # same base_names

  # 5) Group transcripts by identical bitpatterns to form `partitions`
  patterns   <- unique(row_strs_final)
  partitions <- setNames(
    lapply(patterns, function(pat) {
      names(row_strs_final)[ row_strs_final == pat ]
    }),
    patterns
  )

  # 6) Build `paths_unique`: for each distinct pattern, list the internal vertex names
  #    that correspond to a 1 in that pattern string.
  #    First, figure out the _names_ of the kept internal vertices:
  v_names = igraph::V(g)$name
  kept_vertex_names <- v_names[ internal_vertices ][ keep_positions ]

  # Next, for each transcripts final pattern, record which kept vertices it visits.
  paths <- lapply(
    seq_along(row_strs_final),
    function(i) {
      bitstr <- row_strs_final[i]
      bits   <- strsplit(bitstr, "")[[1]]
      kept_vertex_names[ which(bits == "1") ]
    }
  )
  names(paths) <- row_strs_final

  # 7) Now collapse to _unique_ paths (avoiding duplicates)
  paths_unique <- paths[ !duplicated(paths) ]

  #   partitions   : named list where each name = a bitpattern string
  #                  and each element = vector of transcripts with that pattern
  #
  #   paths_unique : list of unique internalvertex sets (character vectors),
  #                  one entry per distinct bitpattern

  return(list(partition = partitions, path = paths_unique))
}








get_bubble_variants_igraph_himem <- function (g, v_start, v_end) 
{
  attrs <- igraph::vertex_attr_names(g)
  excluded <- c("name", "position", "sg_id", "id")
  trans <- setdiff(attrs, excluded)

  v_names = igraph::V(g)$name
  txbase <- unlist(igraph::vertex_attr(g, index=v_start)[trans]) & unlist(igraph::vertex_attr(g, index=v_end)[trans])
  txbase_names <- trans[txbase] 
  internal_nodes = igraph::V(g)[v_names[(as.integer(v_start)+1):(as.integer(v_end)-1)]]$name    # high_mem
 
  for (n in internal_nodes) { 
    tx_internal = unlist(igraph::vertex_attr(g, index=n)[trans])      # high_mem
    if (all(tx_internal >= txbase))
      return (list())
  } 

  bubble_submat <- igraph::vertex_attr(g, index=internal_nodes)[txbase_names]  
  if (length(bubble_submat) == 0)
    return (list())
  bubble_submat <- matrix(unlist(bubble_submat), nrow=length(bubble_submat), byrow=TRUE)
  rownames(bubble_submat) = txbase_names
  colnames(bubble_submat) = internal_nodes
  bubble_submat <- bubble_submat[, colSums(bubble_submat) != 0L, drop = FALSE] # only keep nodes that are present in txbase
  if (length(bubble_submat) == 0)
    return (list())

  bubble_submat <- bubble_submat+0
  bubble_submat_colnames = colnames(bubble_submat)
  row_strs   <- apply(bubble_submat, 1, paste, collapse = "")
  patterns <- unique(row_strs)

  partitions = lapply(seq_len(length(patterns)), function(i) names(row_strs[row_strs == patterns[i]]))
  names(partitions) = patterns
  paths = lapply(seq_len(nrow(bubble_submat)), 
        function(i) bubble_submat_colnames[bubble_submat[i, ]==1])
  names(paths) = row_strs
  paths_unique <- paths[ !duplicated(paths) ]

  list(partition = partitions, path = paths_unique)
}


#' @export
get_bubble_variants_igraph_slow <- function (g, bubble_paths, v_start, v_end) 
{
  internal_paths = lapply(bubble_paths, function(x) { x[2:(length(x)-1)] }) # get rid of first and last (source and sink)

  partitions = list()
  txbase <- trans[unlist(igraph::vertex_attr(g, index=v_start)[trans]) & unlist(igraph::vertex_attr(g, index=v_end)[trans])]
  for (path in internal_paths) {
    tx_pass = igraph::vertex_attr(g, index=path)[txbase]
    tx_pass = names(which(unlist(tx_pass)))   
    partitions = append(partitions, list(tx_pass))
  } 
  list(partition = partitions, path = internal_paths)
}


#' bubble detection between v_start and v_end
#' @export
detect_bubbles_i_j_igraph <- function(v_start_idx, v_end_idx, g, v_start_name, v_end_name)
{
  
    #bubble_paths = bubble_paths_igraph(g, v_start, v_end)
    #print(paste("check", v_start, v_end))
    bubble_variants <- grase::get_bubble_variants_igraph(g, v_start_idx, v_end_idx)
    bubble_d <- length(bubble_variants$partition)
    if (bubble_d <= 1L) 
        return (NULL)
    #print(paste("bubble", v_start, v_end))
    ans_source <- v_start_name
    ans_sink <- v_end_name
    ans_d <- bubble_d
    bubble_partitions <- bubble_variants$partition
    bubble_partitions <- sapply(bubble_partitions, base::paste, 
        collapse = ",")
    bubble_partitions <- paste0("{", bubble_partitions, 
        "}")
    ans_partitions <- IRanges::CharacterList(bubble_partitions)
    bubble_paths <- bubble_variants$path
    bubble_paths <- sapply(bubble_paths, base::paste, 
        collapse = ",")
    bubble_paths <- paste0("{", bubble_paths, "}")
    ans_paths <- IRanges::CharacterList(bubble_paths)
    return (list(ans_source = ans_source, ans_sink = ans_sink, ans_d = ans_d, ans_partitions = ans_partitions, ans_paths = ans_paths ))

}




#' Main bubble detection from matrix
#' @export
detect_bubbles_igraph <- function (g) 
{
  # set transcript vertex_attr 
  if (length(igraph::vertex_attr_names(g)) <= 4)  {
    g <- grase::set_txpath_to_vertex_attr(g)
  }
 
  g_tmp = g
  # lose all expart edges
  ex_parts = igraph::E(g_tmp)[igraph::edge_attr(g_tmp)$ex_or_in == 'ex_part']
  g_tmp = igraph::delete_edges(g_tmp, ex_parts)
  v_names <- igraph::V(g_tmp)$name

  #sgnodetypes <- grase::get_sgnodetypes_igraph(g_tmp)
  ans_source <- ans_sink <- ans_AScode <- character(0)
  ans_d <- integer(0)
  ans_partitions <- ans_paths <- IRanges::CharacterList()
    
  n_nodes = length(v_names)
  for (i in 1:(n_nodes - 2L)) {
    for (j in (i + 2L):n_nodes) {
      v_start_idx <- i
      v_end_idx   <- j
      retval = grase::detect_bubbles_i_j_igraph(i, j, g_tmp, v_names[i], v_names[j])
      if (is.null(retval)) 
        next
      
      ans_source <- c(ans_source, retval$ans_source)
      ans_sink <- c(ans_sink, retval$ans_sink)
      ans_d <- c(ans_d, retval$ans_d)
      ans_partitions <- c(ans_partitions, retval$ans_partitions)
      ans_paths <- c(ans_paths, retval$ans_paths)
    }
  }
  S4Vectors::DataFrame(source = ans_source, sink = ans_sink, d = ans_d, 
    partitions = ans_partitions, paths = ans_paths
    )
}



