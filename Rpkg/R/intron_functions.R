
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


#' Consecutive node pairs, kept separate PER ROUTE
#'
#' Like \code{route_node_pairs} but preserves route identity, which the cut
#' rule in \code{label_bipartition_introns} needs: it takes the first distinct
#' intron along each route, so the routes cannot be flattened together.
#'
#' @param field Character scalar, the path1 or path2 value.
#' @return List of two-column integer matrices, one per route.
#' @export
route_node_pairs_by_route <- function(field) {
  routes <- list()
  for (rt in strsplit(as.character(field), ",")[[1]]) {
    n <- suppressWarnings(as.integer(trimws(strsplit(trimws(rt), "-")[[1]])))
    if (length(n) < 2L) next
    out <- list()
    for (k in seq_len(length(n) - 1L))
      if (!is.na(n[k]) && !is.na(n[k + 1L])) out[[length(out) + 1L]] <- c(n[k], n[k + 1L])
    routes[[length(routes) + 1L]] <- if (length(out)) do.call(rbind, out)
                                     else matrix(integer(0), ncol = 2L)
  }
  routes
}


#' Consecutive node pairs over every route of a path field
#'
#' A path field holds one or more routes separated by commas, each a "-"
#' separated node sequence ("R-1-5-6, R-3-5-6"). Splitting on "-" alone would
#' manufacture a pair spanning the comma boundary, so routes are split first.
#' Non-numeric nodes (the virtual "R"/"L" terminals) are dropped: they carry no
#' intron. Route identity is NOT preserved -- use
#' \code{route_node_pairs_by_route} when it matters.
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
#' AGGREGATION. The historical rule (\code{rule="union"}) summed over ALL of a
#' side's distinct introns, on the assumption that the routes partition
#' themselves among those introns so each contributes once. That assumption is
#' false in practice: measured over the 3,442 junction-substituted significant
#' TSS/TTS sides in the DICE activation panel, 83.8% have at least one
#' transcript crossing two or more of the summed introns, so those molecules
#' are counted two or more times and D is inflated. Because the inflation
#' factor is the composition-weighted mean junctions-per-molecule, it shifts
#' with isoform composition and does NOT cancel in the between-condition
#' contrast.
#'
#' \code{rule="cut"} takes the FIRST distinct intron along each route, so every
#' route contributes exactly one junction and the sum is molecule-proportional,
#' matching the exonic intersection semantics. It is unbiased but discards the
#' route's other k-1 junctions, so it is needlessly high-variance.
#'
#' \code{rule="mean"} (THE DEFAULT) keeps every distinct intron on a route,
#' grouped, and \code{sum_sj_counts} averages within the group before summing
#' across routes. Same expectation as the cut with roughly 1/k the variance,
#' because it uses all k measurements of that route's abundance. Junction ids
#' are emitted route-grouped: "|" within a route, "," between routes. A group of
#' one averages to itself, so ungrouped input behaves exactly as before.
#' Routes with no distinct intron are counted in \code{uncovered1} /
#' \code{uncovered2}; a side with uncovered > 0 is only partly represented by
#' its junction measure.
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
                                      pairs1 = NULL, pairs2 = NULL,
                                      routes1 = NULL, routes2 = NULL,
                                      rule = c("mean", "cut", "union")) {
  rule <- match.arg(rule)
  na_result <- list(distinct1=NA_character_, distinct2=NA_character_,
                    shared=NA_character_, uncovered1=NA_integer_,
                    uncovered2=NA_integer_)

  in_idx <- which(ge$ex_or_in == "in")
  if (length(in_idx) == 0L) return(na_result)

  if (!is.null(bubble_verts) && !is.null(ge$vx_sg_from)) {
    keep <- ge$vx_sg_from[in_idx] %in% bubble_verts & ge$vx_sg_to[in_idx] %in% bubble_verts
    in_idx <- in_idx[keep]
  }
  if (length(in_idx) == 0L) return(na_result)

  use_paths <- !is.null(pairs1) && !is.null(pairs2) && !is.null(ge$vx_sg_from)
  ekey <- if (!is.null(ge$vx_sg_from))
            paste(ge$vx_sg_from[in_idx], ge$vx_sg_to[in_idx]) else character(0)
  if (use_paths) {
    ## an intron is "on" a side iff its endpoint pair is a consecutive pair on
    ## one of that side's routes
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

  ## CUT RULE -- one intron per route.
  ## A side's molecules all enter at the bubble source and leave at the sink,
  ## so summing over every distinct intron counts a molecule once per intron it
  ## crosses. Taking the FIRST distinct intron along each route instead gives a
  ## cut: every route contributes exactly one junction, so the sum is
  ## molecule-proportional and comparable to the exonic (intersection) set.
  ## `uncovered` counts routes with no distinct intron at all -- those
  ## molecules are invisible to the junction measure and the side should not be
  ## substituted on its strength.
  ## `all_per_route = FALSE` gives the CUT (first distinct intron per route).
  ## TRUE gives every distinct intron per route, kept grouped, for the MEAN
  ## rule: averaging within a route has the same expectation as taking one of
  ## its junctions but roughly 1/k the variance, because it uses all k
  ## measurements of the same route abundance instead of discarding k-1.
  cut_for_side <- function(routes, keep, all_per_route = FALSE) {
    if (is.null(routes) || !length(routes))
      return(list(sel = integer(0), groups = list(), unc = NA_integer_))
    ok <- ekey[keep]
    idx_keep <- in_idx[keep]
    sel <- integer(0); groups <- list(); unc <- 0L
    for (P in routes) {
      if (!nrow(P)) { unc <- unc + 1L; next }
      hits <- integer(0)
      for (k in seq_len(nrow(P))) {
        key <- paste(P[k, 1L], P[k, 2L])
        j <- which(ok == key)
        if (length(j)) {
          hits <- c(hits, idx_keep[j[1L]])
          if (!all_per_route) break
        }
      }
      if (!length(hits)) { unc <- unc + 1L; next }
      sel <- c(sel, hits)
      groups[[length(groups) + 1L]] <- unique(hits)
    }
    list(sel = unique(sel), groups = groups, unc = unc)
  }

  ## Route-grouped junction ids: "|" separates junctions WITHIN a route, ","
  ## separates routes. sum_sj_counts() averages within a group and sums across
  ## groups, so a group of one reduces exactly to the old behaviour.
  make_junction_groups <- function(groups) {
    if (!length(groups)) return(NA_character_)
    out <- character(0)
    for (g in groups) {
      fp <- ge$from_pos[g]; tp <- ge$to_pos[g]
      ok2 <- !is.na(fp) & !is.na(tp)
      if (!any(ok2)) next
      jids <- paste0(chr, ":", pmin(fp[ok2], tp[ok2]), ":", pmax(fp[ok2], tp[ok2]))
      out <- c(out, paste(sort(unique(jids)), collapse = "|"))
    }
    if (!length(out)) return(NA_character_)
    paste(out, collapse = ",")
  }

  d1 <- in_tx1 & !in_tx2
  d2 <- !in_tx1 & in_tx2
  if (rule %in% c("cut", "mean") && use_paths &&
      !is.null(routes1) && !is.null(routes2)) {
    per_route <- rule == "mean"
    c1 <- cut_for_side(routes1, d1, per_route)
    c2 <- cut_for_side(routes2, d2, per_route)
    emit <- if (per_route) function(x) make_junction_groups(x$groups)
            else           function(x) make_junctions(x$sel)
    return(list(
      distinct1  = emit(c1),
      distinct2  = emit(c2),
      shared     = make_junctions(in_idx[in_tx1 & in_tx2]),
      uncovered1 = c1$unc,
      uncovered2 = c2$unc
    ))
  }

  list(
    distinct1  = make_junctions(in_idx[d1]),
    distinct2  = make_junctions(in_idx[d2]),
    shared     = make_junctions(in_idx[in_tx1 & in_tx2]),
    uncovered1 = NA_integer_,
    uncovered2 = NA_integer_
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
label_all_bipartition_introns <- function(splits_df, graphml_dir, chr_map,
                                          rule = c("mean", "cut", "union")) {
  rule <- match.arg(rule)
  splits_df$intron_distinct1 <- NA_character_
  splits_df$intron_distinct2 <- NA_character_
  splits_df$intron_shared    <- NA_character_
  splits_df$intron_uncovered1 <- NA_integer_
  splits_df$intron_uncovered2 <- NA_integer_

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
      routes1 <- NULL; routes2 <- NULL
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
          routes1 <- route_node_pairs_by_route(p1_raw)
          routes2 <- route_node_pairs_by_route(p2_raw)
        }
      }
      lbl <- label_bipartition_introns(ge, tx1_set, tx2_set, chr, bubble_verts,
                                       pairs1, pairs2, routes1, routes2, rule)
      splits_df$intron_distinct1[i]  <- lbl$distinct1
      splits_df$intron_distinct2[i]  <- lbl$distinct2
      splits_df$intron_shared[i]     <- lbl$shared
      splits_df$intron_uncovered1[i] <- lbl$uncovered1
      splits_df$intron_uncovered2[i] <- lbl$uncovered2
    }
  }
  splits_df
}
