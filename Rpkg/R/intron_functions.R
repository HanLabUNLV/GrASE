
#' Parse gene-to-chromosome map from a DEXSeq GFF file
#'
#' Reads only the aggregate_gene lines (one per gene) for efficiency.
#'
#' @param gff_path Path to a combined DEXSeq GFF file.
#' @return Named character vector: names are gene_ids, values are chromosomes.
#' @export
parse_gff_chr_map <- function(gff_path) {
  lines <- readLines(gff_path)
  lines <- lines[grepl("\taggregate_gene\t", lines, fixed=TRUE)]
  if (length(lines) == 0L) stop("no aggregate_gene lines found in: ", gff_path)
  parts    <- strsplit(lines, "\t")
  chr_col  <- vapply(parts, `[[`, character(1L), 1L)
  attr_col <- vapply(parts, `[[`, character(1L), 9L)
  gene_ids <- sub('.*gene_id "([^"]+)".*', '\\1', attr_col)
  map <- chr_col
  names(map) <- gene_ids
  map
}


#' Consecutive node pairs over every route of a path field
#'
#' A path field holds one or more routes separated by commas, each a "-"
#' separated node sequence ("R-1-5-6, R-3-5-6"). Splitting on "-" alone would
#' manufacture a pair spanning the comma boundary, so routes are split first.
#' Non-numeric nodes (the virtual "R"/"L" terminals) are dropped: they carry no
#' intron.
#'
#' @param field Character scalar, the path1 or path2 value.
#' @return Two-column integer matrix of (from, to) sg_id pairs; zero rows if none.
#' @export
route_node_pairs <- function(field) {
  out <- list()
  for (rt in strsplit(as.character(field), ",")[[1]]) {
    n <- suppressWarnings(as.integer(trimws(strsplit(trimws(rt), "-")[[1]])))
    if (length(n) < 2L) next
    for (k in seq_len(length(n) - 1L))
      if (!is.na(n[k]) && !is.na(n[k + 1L])) out[[length(out) + 1L]] <- c(n[k], n[k + 1L])
  }
  if (!length(out)) return(matrix(integer(0), ncol = 2L))
  do.call(rbind, out)
}


#' Label intronic edges of one bipartition as distinct or shared
#'
#' Uses the precomputed gene graph cache from \code{precompute_gene_graph} to
#' classify each intronic ("in") edge as distinct to path1, distinct to path2, or
#' shared.
#'
#' PATH RULE (default when \code{pairs1}/\code{pairs2} are supplied): an intron
#' is distinct to side 1 iff it lies on at least one of side 1's ROUTES and on
#' none of side 2's. This mirrors how the exonic distinct sets are built in
#' \code{bipartition_analysis.R}, which derives them from the route node
#' sequences rather than from transcript membership.
#'
#' TRANSCRIPT RULE (fallback when no route pairs are given): an intron is
#' distinct to side 1 iff SOME side-1 transcript uses it and NO side-2
#' transcript does, restricted to introns whose endpoints are both bubble
#' vertices. This was the original rule. It is a strict SUPERSET of the path
#' rule: measured over 2,272 sides, the two agree on 95.8% overall and on
#' 100.0% of the sides that are actually substituted (those with an empty exonic
#' distinct set). Every disagreement is a transcript in \code{transcripts1}
#' whose route through the bubble is not among the listed \code{path1} routes.
#'
#' Either way the junction identity is the same: "chr:start:end" from the edge
#' from_pos/to_pos, and a side's count sums over ALL its distinct introns, which
#' is what makes the intron D set an equivalent measure to the exonic one -- the
#' routes partition themselves among those introns, so every transcript of the
#' partition contributes.
#'
#' @param ge  Precomputed gene graph list from \code{precompute_gene_graph}.
#' @param tx1_set Character vector of transcript IDs in path1.
#' @param tx2_set Character vector of transcript IDs in path2.
#' @param chr Chromosome string (e.g. "chr1").
#' @param bubble_verts Integer vector of sg_id values for all vertices in the
#'   bubble (source through sink, both paths). NULL disables the filter. Used by
#'   the transcript rule only; the path rule is already route-bounded.
#' @param pairs1,pairs2 Two-column integer matrices of consecutive sg_id node
#'   pairs over every route of path1 / path2 (see
#'   \code{label_all_bipartition_introns}). Supplying both selects the path rule.
#' @return Named list with elements distinct1, distinct2, shared -- each a
#'   comma-separated string of "chr:start:end" junction identifiers, or NA.
#' @export
label_bipartition_introns <- function(ge, tx1_set, tx2_set, chr, bubble_verts = NULL,
                                      pairs1 = NULL, pairs2 = NULL) {
  na_result <- list(distinct1=NA_character_, distinct2=NA_character_, shared=NA_character_)

  in_idx <- which(ge$ex_or_in == "in")
  if (length(in_idx) == 0L) return(na_result)

  if (!is.null(bubble_verts) && !is.null(ge$vx_sg_from)) {
    keep <- ge$vx_sg_from[in_idx] %in% bubble_verts & ge$vx_sg_to[in_idx] %in% bubble_verts
    in_idx <- in_idx[keep]
  }
  if (length(in_idx) == 0L) return(na_result)

  use_paths <- !is.null(pairs1) && !is.null(pairs2) && !is.null(ge$vx_sg_from)
  if (use_paths) {
    ## an intron is "on" a side iff its endpoint pair is a consecutive pair on
    ## one of that side's routes
    ekey <- paste(ge$vx_sg_from[in_idx], ge$vx_sg_to[in_idx])
    pkey <- function(P) if (is.null(P) || !nrow(P)) character(0) else
                        unique(paste(P[, 1L], P[, 2L]))
    in_tx1 <- ekey %in% pkey(pairs1)
    in_tx2 <- ekey %in% pkey(pairs2)
    if (!any(in_tx1) && !any(in_tx2)) return(na_result)
  } else {
    tx1_mask <- ge$tx_cols %in% tx1_set
    tx2_mask <- ge$tx_cols %in% tx2_set
    if (!any(tx1_mask) || !any(tx2_mask)) return(na_result)

    sub_mat <- ge$tx_mat[in_idx, , drop=FALSE]
    in_tx1 <- rowSums(sub_mat[, tx1_mask, drop=FALSE]) > 0L
    in_tx2 <- rowSums(sub_mat[, tx2_mask, drop=FALSE]) > 0L
  }

  make_junctions <- function(sel_idx) {
    if (length(sel_idx) == 0L) return(NA_character_)
    fp <- ge$from_pos[sel_idx]
    tp <- ge$to_pos[sel_idx]
    ok <- !is.na(fp) & !is.na(tp)
    if (!any(ok)) return(NA_character_)
    jids <- paste0(chr, ":", pmin(fp[ok], tp[ok]), ":", pmax(fp[ok], tp[ok]))
    paste(sort(unique(jids)), collapse=",")
  }

  list(
    distinct1 = make_junctions(in_idx[in_tx1 & !in_tx2]),
    distinct2 = make_junctions(in_idx[!in_tx1 & in_tx2]),
    shared    = make_junctions(in_idx[in_tx1 & in_tx2])
  )
}


#' Label intronic edges for all rows of a bipartition splits data frame
#'
#' For each bipartition, loads the gene's graphml, calls
#' \code{label_bipartition_introns}, and appends three columns:
#' \code{intron_distinct1}, \code{intron_distinct2}, \code{intron_shared}.
#' Each column holds comma-separated "chr:start:end" junction identifiers.
#'
#' @param splits_df Data frame with bipartition splits (must have columns:
#'   gene, transcripts1, transcripts2).
#' @param graphml_dir Directory containing per-gene .dexseq.graphml files,
#'   named as \code{<gene_id>.dexseq.graphml}.
#' @param chr_map Named character vector from \code{parse_gff_chr_map}.
#' @return \code{splits_df} extended with intron_distinct1, intron_distinct2,
#'   intron_shared columns.
#' @export
label_all_bipartition_introns <- function(splits_df, graphml_dir, chr_map) {
  splits_df$intron_distinct1 <- NA_character_
  splits_df$intron_distinct2 <- NA_character_
  splits_df$intron_shared    <- NA_character_

  for (gid in unique(splits_df$gene)) {
    graphml_path <- file.path(graphml_dir, paste0(gid, ".graphml"))
    if (!file.exists(graphml_path)) {
      warning("graphml not found for gene: ", gid)
      next
    }
    g   <- igraph::read_graph(graphml_path, format="graphml")
    ge  <- precompute_gene_graph(g)
    chr <- chr_map[[gid]]
    if (is.null(chr) || is.na(chr)) {
      warning("chromosome not found for gene: ", gid)
      next
    }

    has_paths <- "path1" %in% names(splits_df) && "path2" %in% names(splits_df)
    rows <- which(splits_df$gene == gid)
    for (i in rows) {
      t1_raw <- splits_df$transcripts1[i]
      t2_raw <- splits_df$transcripts2[i]
      if (is.na(t1_raw) || is.na(t2_raw)) next
      tx1_set <- trimws(unlist(strsplit(t1_raw, ",")))
      tx2_set <- trimws(unlist(strsplit(t2_raw, ",")))
      tx1_set <- tx1_set[nchar(tx1_set) > 0L]
      tx2_set <- tx2_set[nchar(tx2_set) > 0L]
      bubble_verts <- NULL; pairs1 <- NULL; pairs2 <- NULL
      if (has_paths) {
        p1_raw <- splits_df$path1[i]
        p2_raw <- splits_df$path2[i]
        if (!is.na(p1_raw) && !is.na(p2_raw)) {
          p1 <- suppressWarnings(as.integer(unlist(strsplit(p1_raw, "-"))))
          p2 <- suppressWarnings(as.integer(unlist(strsplit(p2_raw, "-"))))
          bubble_verts <- unique(c(p1, p2))
          bubble_verts <- bubble_verts[!is.na(bubble_verts)]
          ## consecutive node pairs per ROUTE -- split on "," FIRST so a pair is
          ## never formed across two routes ("...-6, R-3-..." must not yield 6->R)
          pairs1 <- route_node_pairs(p1_raw)
          pairs2 <- route_node_pairs(p2_raw)
        }
      }
      lbl <- label_bipartition_introns(ge, tx1_set, tx2_set, chr, bubble_verts,
                                       pairs1, pairs2)
      splits_df$intron_distinct1[i] <- lbl$distinct1
      splits_df$intron_distinct2[i] <- lbl$distinct2
      splits_df$intron_shared[i]    <- lbl$shared
    }
  }
  splits_df
}
