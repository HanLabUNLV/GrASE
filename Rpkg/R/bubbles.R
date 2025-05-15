# R/bubbles.R


#' Generate transcript-paths from igraph edge attributes
#' @export
txpath_from_edgeattr <- function(g, type="exon") {
  # Identify transcript attributes on edges (exclude known graph attrs)
  attrs <- igraph::edge_attr_names(g)
  excluded <- c("sgedge_id", "ex_or_in", "from_pos", "to_pos", "dexseq_fragment")
  trans <- setdiff(attrs, excluded)
  # Determine root and leaf vertices
  root <- if ("R" %in% igraph::V(g)$name) "R" else stop("No root 'R' vertex found")
  leaf <- if ("L" %in% igraph::V(g)$name) "L" else stop("No leaf 'L' vertex found")

  txpaths <- setNames(vector("list", length(trans)), trans)
  for (tx in trans) {
    tx_bools <- igraph::edge_attr(g, tx)
    ex_bools <- igraph::edge_attr(g, 'ex_or_in') != 'ex_part' 
    expart_bools <- igraph::edge_attr(g, 'ex_or_in') != 'ex'
    exon_ids = igraph::E(g)[which(tx_bools & ex_bools)]
    print (paste("exon_ids", paste(exon_ids)))
    expart_ids = igraph::E(g)[which(tx_bools & expart_bools)]
    if (type == "exon") {
      subg  <- igraph::subgraph_from_edges(g, exon_ids, delete.vertices = FALSE)
      print (paste("subg", subg))
    }
    else if (type == "expart") {
      subg  <- igraph::subgraph_from_edges(g, expart_ids, delete.vertices = FALSE)
    }
    else {
      print("error: type should be either exon or expart")
      return (NULL)
    }
    spath    <- igraph::shortest_paths(subg, from = root, to = leaf, mode = "out")$vpath[[1]]
    txpaths[[tx]] <- igraph::V(g)$name[spath]
    print(txpaths)
  }
  txpaths = lapply(txpaths, function(x) { x[2:(length(x)-1)] }) # get rid of first and last (R and L)
  txpaths
}



#' Build transcript-path logical matrix from index lists
#' @export
make_matrix_from_txpath_igraph <- function(txpath_vertex_list) {
  SSids <- as.numeric( sort(unique(unlist(txpath_vertex_list))))
  SSids <- SSids[!is.na(SSids)]
  SSids <- SSids[order(SSids)]
  cols  <- c("R", as.character(SSids), "L")
  mat <- matrix(FALSE,
                nrow = length(txpath_vertex_list),
                ncol = length(cols),
                dimnames = list(names(txpath_vertex_list), cols))
  mat[, "R"] <- TRUE
  mat[, "L"] <- TRUE
  for (i in seq_along(txpath_vertex_list)) {
    vids <- txpath_vertex_list[[i]]
    cmatch <- match(as.character(vids), cols)
    mat[i, cmatch] <- TRUE
  }
  mat
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
    txbase <- txpathmat[, i] & txpathmat[, j]
    bubble_submat <- txpathmat[txbase, (i + 1L):(j - 1L), drop = FALSE]
    bubble_submat <- bubble_submat[, colSums(bubble_submat) != 
        0L, drop = FALSE]
    bubble_submat_rownames <- rownames(bubble_submat)
    bubble_submat_colnames <- colnames(bubble_submat)
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
    bubble_AScode <- paste(bubble_variants[, "code"], 
        collapse = ",")
    ans_AScode <- bubble_AScode
    return (list(ans_source = ans_source, ans_sink = ans_sink, ans_d = ans_d, ans_AScode = ans_AScode, ans_partitions = ans_partitions, ans_paths = ans_paths ))

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
          ans_AScode <- c(ans_AScode, retval$ans_AScode)
        }
    }
    S4Vectors::DataFrame(source = ans_source, sink = ans_sink, d = ans_d, 
        partitions = ans_partitions, paths = ans_paths, AScode = ans_AScode 
        )
}



#' Wrapper: detect bubbles on igraph DAG
#' @export
detect_bubbles_igraph <- function(g) {
  if (!igraph::is_directed(g)) stop("Graph must be directed DAG")
  # generate txpaths
  tx_exon_paths <- grase::txpath_from_edgeattr(g)
  #tx_expart_paths <- grase::txpath_from_edgeattr(g, "expart")
  txmat <- grase::make_matrix_from_txpath_igraph(tx_exon_paths)
  bubble_df <- grase::detect_bubbles_from_mat(txmat)
  bubble_df
}

