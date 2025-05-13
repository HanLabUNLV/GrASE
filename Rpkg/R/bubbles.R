# bubble_detection_igraph.R
# Pure-R replication of SplicingGraphs path & bubble extraction for igraph DAGs

library(igraph)

# 1. Generate transcript-paths from igraph edge attributes
txpath_from_edgeattr <- function(g) {
  # Identify transcript attributes on edges (exclude known graph attrs)
  attrs <- edge_attr_names(g)
  excluded <- c("sgedge_id", "ex_or_in", "from_pos", "to_pos", "dexseq_fragment")
  trans <- setdiff(attrs, excluded)
  # Determine root and leaf vertices
  root <- if ("R" %in% V(g)$name) "R" else stop("No root 'R' vertex found")
  leaf <- if ("L" %in% V(g)$name) "L" else stop("No leaf 'L' vertex found")

  txpaths <- setNames(vector("list", length(trans)), trans)
  for (tx in trans) {
    bools <- edge_attr(g, tx)
    eids  <- E(g)[which(bools)]
    subg  <- subgraph_from_edges(g, eids, delete.vertices = FALSE)
    sp    <- shortest_paths(subg, from = root, to = leaf, mode = "out")$vpath[[1]]
    txpaths[[tx]] <- V(g)$name[sp]
  }
  txpaths
}

# 2. Build transcript-path logical matrix from index lists
make_matrix_from_txpath_igraph <- function(txpath_idx_list) {
  SSids <- sort(unique(unlist(txpath_idx_list)))
  cols  <- c("R", as.character(SSids), "L")
  mat <- matrix(FALSE,
                nrow = length(txpath_idx_list),
                ncol = length(cols),
                dimnames = list(names(txpath_idx_list), cols))
  mat[, "R"] <- TRUE
  mat[, "L"] <- TRUE
  for (i in seq_along(txpath_idx_list)) {
    vids <- txpath_idx_list[[i]]
    cmatch <- match(as.character(vids), cols)
    mat[i, cmatch] <- TRUE
  }
  mat
}

# 3. Bubble-extraction helpers
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

is_bubble <- function (txpathmat, i, j) 
{
    txbase <- txpathmat[, i] & txpathmat[, j]
    if (sum(txbase) <= 1L) 
        return(FALSE)
    for (k in (i + 1L):(j - 1L)) if (all(txpathmat[, k] >= txbase)) 
        return(FALSE)
    TRUE
}




get_bubble_variants <- function (txpathmat, sgnodetypes, i, j) 
{
    txbase <- txpathmat[, i] & txpathmat[, j]
    bubble_submat <- txpathmat[txbase, (i + 1L):(j - 1L), drop = FALSE]
    bubble_submat <- bubble_submat[, colSums(bubble_submat) != 
        0L, drop = FALSE]
    bubble_submat_rownames <- rownames(bubble_submat)
    bubble_submat_colnames <- colnames(bubble_submat)
    ans_path <- CharacterList(lapply(seq_len(nrow(bubble_submat)), 
        function(i) bubble_submat_colnames[bubble_submat[i, ]]))
    ans_rpath <- IntegerList(lapply(seq_len(nrow(bubble_submat)), 
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
    ans_partition <- unname(splitAsList(bubble_submat_rownames, 
        cumsum(is_not_dup)))
    DataFrame(partition = ans_partition, path = ans_path, code = ans_code)
}



# 4. Main bubble detection from matrix
detect_bubbles_from_mat <- function (txpathmat) 
{
    prev_locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", prev_locale))
    sgnodetypes <- get_sgnodetypes(txpathmat)
    ans_source <- ans_sink <- ans_AScode <- character(0)
    ans_d <- integer(0)
    ans_partitions <- ans_paths <- CharacterList()
    ncol <- ncol(txpathmat)
    for (i in 1:(ncol - 2L)) {
        for (j in (i + 2L):ncol) {
            if (!is_bubble(txpathmat, i, j)) 
                next
            bubble_variants <- get_bubble_variants(txpathmat, 
                sgnodetypes, i, j)
            bubble_d <- nrow(bubble_variants)
            if (bubble_d <= 1L) 
                next
            ans_source <- c(ans_source, colnames(txpathmat)[i])
            ans_sink <- c(ans_sink, colnames(txpathmat)[j])
            ans_d <- c(ans_d, bubble_d)
            bubble_partitions <- bubble_variants[, "partition"]
            bubble_partitions <- sapply(bubble_partitions, base::paste, 
                collapse = ",")
            bubble_partitions <- paste0("{", bubble_partitions, 
                "}")
            ans_partitions <- c(ans_partitions, CharacterList(bubble_partitions))
            bubble_paths <- bubble_variants[, "path"]
            bubble_paths <- sapply(bubble_paths, base::paste, 
                collapse = ",")
            bubble_paths <- paste0("{", bubble_paths, "}")
            ans_paths <- c(ans_paths, CharacterList(bubble_paths))
            bubble_AScode <- paste(bubble_variants[, "code"], 
                collapse = ",")
            ans_AScode <- c(ans_AScode, bubble_AScode)
        }
    }
    DataFrame(source = ans_source, sink = ans_sink, d = ans_d, 
        partitions = ans_partitions, paths = ans_paths, AScode = ans_AScode, 
        description = unname(ASCODE2DESC[ans_AScode]))
}



# 5. Wrapper: detect bubbles on igraph DAG
detect_bubbles_igraph <- function(g, tx_exons, node2pos) {
  if (!is_directed(g)) stop("Graph must be directed DAG")
  # generate txpaths
  txpaths <- txpath_from_edgeattr(g)
  name2vid <- setNames(seq_along(V(g)$name), V(g)$name)
  idx_list <- lapply(txpaths, function(p) name2vid[p])
  txmat <- make_matrix_from_txpath_igraph(idx_list)
  detect_bubbles_from_mat(txmat)
}

