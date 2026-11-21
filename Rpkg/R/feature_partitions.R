## Partitions of a bubble's paths read off their features: an alternative to the
## lower-set enumeration in bipartition_paths() (valid_partitions_idx +
## remove_symmetric_splits) that gives the same partitions. bipartition_paths() is
## unchanged; adopting this there replaces its enumeration with feature_partitions_idx()
## and makes bipartition_paths_by_feature() redundant.

#' Partitions of a bubble's paths induced by their features
#'
#' A feature (an exonic part or an intron edge) is in the distinct set of side 1 exactly
#' when every path of side 1 carries it and no path of side 2 does, that is, when side 1
#' is the set of paths that carry it. Every partition with a non-empty distinct set is
#' therefore the paths carrying some feature against the rest, and each such feature
#' gives one. The partitions are read off the features in time linear in their number,
#' where the lower sets of the subset order can number 2^K - 2 for K paths.
#'
#' The sides follow \code{remove_symmetric_splits()}: path ids sorted as strings, and
#' group1 the side holding the first of them, so the rows written for a partition match
#' those of the lower-set enumeration.
#'
#' @param route_features named list, per path, of its features (exonic part edges and
#'   intron edges); names are the path ids
#' @return list of list(group1, group2) of path ids, one per distinct feature pattern
#' @export
#' @examples
#' rf <- list("1" = c("a", "b", "x:in"), "2" = c("a", "c"), "3" = c("a", "b", "y:in"))
#' grase::feature_partitions_idx(rf)   # {1,3}|{2} by b and c, {1}|{2,3} by x:in, {1,2}|{3} by y:in
feature_partitions_idx <- function(route_features) {
  ids <- names(route_features)
  if (is.null(ids)) ids <- as.character(seq_along(route_features))
  n <- length(ids)
  feats <- unique(unlist(route_features, use.names = FALSE))
  if (n < 2 || !length(feats)) return(list())
  memb <- matrix(vapply(route_features, function(r) feats %in% r, logical(length(feats))),
                 nrow = length(feats))
  k <- rowSums(memb)
  memb <- memb[k > 0 & k < n, , drop = FALSE]
  if (!nrow(memb)) return(list())
  first <- which(ids == sort(ids)[1])
  flip <- !memb[, first]
  memb[flip, ] <- !memb[flip, , drop = FALSE]
  memb <- unique(memb)
  lapply(seq_len(nrow(memb)), function(i)
    list(group1 = sort(ids[memb[i, ]]), group2 = sort(ids[!memb[i, ]])))
}

#' Bipartitions of every bubble, with the partitions read off the paths' features
#'
#' \code{bipartition_paths()} with \code{feature_partitions_idx()} in place of the
#' lower-set enumeration, and no bubble collapsing. The partitions kept are the same;
#' rows of one bubble may come out in a different order, which renumbers events.
#'
#' @inheritParams bipartition_paths
#' @export
bipartition_paths_by_feature <- function(gene, g, outdir, max_path = 15, max_span = Inf) {

  kept_rows <- list(); kept_in1 <- list(); kept_in2 <- list()
  g <- grase::set_edge_names(g)
  bubbles_df <- grase::detect_bubbles_igraph(g)
  if (nrow(bubbles_df) == 0) {
    message(paste("single transcript, no bubbles ", gene))
    return(NULL)
  }
  bubbles_ordered <- grase::bubble_ordering4(g, bubbles_df)

  ex_or_in_vec <- igraph::edge_attr(g, "ex_or_in")
  g_exon <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == "ex_part"])
  g_expart <- igraph::delete_edges(g, igraph::E(g)[ex_or_in_vec == "ex"])
  ex_or_in_gexpart <- igraph::edge_attr(g_expart, "ex_or_in")
  names(ex_or_in_gexpart) <- igraph::edge_attr(g_expart, "name")
  dexseq_frag <- igraph::edge_attr(g_expart, "dexseq_fragment")
  names(dexseq_frag) <- igraph::edge_attr(g_expart, "name")
  txpaths <- grase::txpath_from_edgeattr(g)
  tx_ex_parts <- lapply(lapply(txpaths, grase::from_vpath_to_exon_path_simple, g = g_exon),
                        grase::from_epath_to_expart_path_simple, g = g_expart)

  message(sprintf("[%s] %s: processing %d bubbles", format(Sys.time(), "%H:%M:%S"), gene,
                  nrow(bubbles_ordered)))
  for (bubble_idx in seq_len(nrow(bubbles_ordered))) {
    source <- bubbles_ordered$source[[bubble_idx]]
    sink <- bubbles_ordered$sink[[bubble_idx]]
    num_paths <- bubbles_ordered$num_paths[[bubble_idx]]
    span <- bubbles_ordered$span[[bubble_idx]]
    if (source == "R" && sink == "L") next
    if (num_paths > max_path || span > max_span) {
      message(paste("bubble with ", num_paths, " alternative paths and span ", span,
                    " exonic parts: skipping ", gene, source, sink))
      next
    }
    parsed_paths <- lapply(bubbles_ordered$paths[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    parsed_partitions <- lapply(bubbles_ordered$partitions[[bubble_idx]], function(x) strsplit(gsub("[{}]", "", x), ",")[[1]])
    vpaths <- lapply(parsed_paths, function(vec) c(source, vec, sink))
    epaths <- lapply(vpaths, grase::from_vpath_to_exon_path_simple, g = g_exon)
    ex_part_paths <- lapply(epaths, grase::from_epath_to_expart_path_simple, g = g_expart)
    ex_part_paths <- lapply(ex_part_paths, function(x) x[ex_or_in_gexpart[x] == "ex_part"])
    names(ex_part_paths) <- seq_along(ex_part_paths)
    route_features <- mapply(function(p, e) c(p, e[grepl(":in$", e)]), ex_part_paths, epaths,
                             SIMPLIFY = FALSE)
    names(route_features) <- names(ex_part_paths)

    for (split in feature_partitions_idx(route_features)) {
      retval <- grase::find_diff_and_ref_exparts_for_split(
        g = g, source = source, sink = sink, split = split,
        parsed_partitions = parsed_partitions, parsed_paths = parsed_paths,
        tx_ex_parts = tx_ex_parts, bipartitions_df = data.frame(), gene = gene,
        g_exon = g_exon, g_expart = g_expart,
        ex_or_in_gexpart = ex_or_in_gexpart, dexseq_frag = dexseq_frag,
        route_epaths = epaths, route_ex_parts = ex_part_paths)
      if (!retval$has_distinct) next     # never taken: every feature partition has one
      kept_rows[[length(kept_rows) + 1L]] <- retval$new_row
      kept_in1[[length(kept_in1) + 1L]] <- retval$in_setdiff1
      kept_in2[[length(kept_in2) + 1L]] <- retval$in_setdiff2
    }
  }

  bipartitions_df <- if (length(kept_rows)) do.call(rbind, kept_rows) else data.frame()
  in_distinct1 <- as.character(unlist(kept_in1)); in_distinct2 <- as.character(unlist(kept_in2))

  ## from here as in bipartition_paths(): dedup on the measured sets, closest reference first
  if (nrow(bipartitions_df) > 0) {
    bipartitions_df <- bipartitions_df %>% dplyr::mutate(
      ref_part_cnt = sapply(ref_ex_part, count_items) + (source == "R") + (sink == "L"),
      .terminal    = source == "R" | sink == "L",
      .key1 = ifelse(setdiff1 != "", setdiff1, paste0("in:", in_distinct1)),
      .key2 = ifelse(setdiff2 != "", setdiff2, paste0("in:", in_distinct2)))
    bipartitions_filtered <- bipartitions_df %>%
      dplyr::distinct(ref_ex_part, .key1, .key2, .keep_all = TRUE) %>%
      dplyr::group_by(.key1, .key2) %>%
      dplyr::mutate(.prox = mapply(ref_proximity_to_setdiffs,
                                   ref_ex_part, setdiff1, setdiff2)) %>%
      dplyr::filter(.prox == min(.prox, na.rm = TRUE)) %>%
      dplyr::filter(ref_part_cnt == min(ref_part_cnt)) %>%
      dplyr::arrange(.terminal, .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(-ref_part_cnt, -.prox, -.terminal, -.key1, -.key2)
    bipartitions_df <- dplyr::select(bipartitions_df, -.key1, -.key2)
  } else {
    bipartitions_filtered <- bipartitions_df
  }

  utils::write.table(bipartitions_df, file.path(outdir, paste0(gene, ".bipartition.all.txt")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(bipartitions_filtered, file.path(outdir, paste0(gene, ".bipartition.txt")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  bipartitions_filtered
}
