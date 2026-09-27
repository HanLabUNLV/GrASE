## Edge-coverage diagnostic for alternative TSS / TTS calls.
##
## A TSS/TTS bipartition is defined from annotation: one side's transcripts
## begin (or end) at a different terminal exon. That does not establish that
## transcription actually starts or stops there -- signal continuing in from
## outside, for example a retained intron, produces the same proportion change.
##
## The diagnostic asks the direct question. If the boundary is real, coverage
## steps down at the free edge of the alternative terminal exon. If the signal
## is continuation, coverage runs through it.
##
##   step = adj_n / (adj_n + d_n)     0   coverage stops dead at the boundary
##                                    0.5 perfect continuation
##
## This is an ANNOTATION on existing calls. It does not change which sides are
## called, only how a call should be read.


#' Which end of a bipartition side is free
#'
#' A TSS bubble is anchored at the graph source ("R"), a TTS bubble at the sink
#' ("L"). Combined with gene strand this fixes which end of the side's distinct
#' set faces away from the transcript body.
#'
#' @param kind "TSS" or "TTS".
#' @param strand "+" or "-".
#' @return -1 if the free end is at lower coordinates, +1 if higher.
#' @export
free_end_direction <- function(kind, strand) {
  stopifnot(kind %in% c("TSS", "TTS"), strand %in% c("+", "-"))
  if ((kind == "TSS") == (strand == "+")) -1L else 1L
}


#' Distance from a boundary to the nearest exonic part
#'
#' @param parts Named list of exonic parts, part number to c(start, end),
#'   1-based inclusive -- as returned by \code{graph_exonic_parts} reshaped by
#'   part, or read from the flattened GFF.
#' @param boundary Integer coordinate.
#' @param direction -1 or +1.
#' @return Number of bases available before the next part, or a large value if
#'   the boundary faces the edge of the gene.
#' @export
flank_width <- function(parts, boundary, direction) {
  if (direction < 0) {
    nb <- vapply(parts, `[`, numeric(1), 2)
    nb <- nb[nb < boundary]
    if (!length(nb)) return(.Machine$integer.max)
    as.integer(boundary - max(nb) - 1L)
  } else {
    nb <- vapply(parts, `[`, numeric(1), 1)
    nb <- nb[nb > boundary]
    if (!length(nb)) return(.Machine$integer.max)
    as.integer(min(nb) - boundary - 1L)
  }
}


#' Flanking bin just outside a boundary
#'
#' @param boundary Integer coordinate of the free edge.
#' @param direction -1 or +1.
#' @param width Nominal bin width in bases.
#' @param avail Available flank width, from \code{flank_width}.
#' @param min_width Below this the side is not scorable.
#' @return List with start, end, width, tier ("full", "clipped",
#'   "unscorable").
#' @export
flank_bin <- function(boundary, direction, width = 100L, avail = width,
                      min_width = 50L) {
  if (is.na(avail) || avail < min_width)
    return(list(start = NA_integer_, end = NA_integer_, width = NA_integer_,
                tier = "unscorable"))
  w <- min(width, avail)
  if (direction < 0) {
    s <- boundary - w; e <- boundary - 1L
  } else {
    s <- boundary + 1L; e <- boundary + w
  }
  list(start = as.integer(s), end = as.integer(e), width = as.integer(w),
       tier = if (w >= width) "full" else "clipped")
}


#' Edge step statistic
#'
#' Counts are put on a common per-100bp footing so bins of different lengths
#' are comparable. pi is computed by the package on RAW counts, so a length
#' mismatch between the flanking bin and the distinct set would otherwise bias
#' the comparison directly.
#'
#' \code{pi_shadow_n} is returned for inspection but should NOT be used for
#' scoring: it degenerates whenever \code{d_n} and \code{s_n} differ greatly in
#' magnitude, because both pi terms pin near 0 or 1 and their difference
#' collapses. Mid-exon calibration bins show this clearly.
#'
#' @param adj,d,s Counts in the flanking bin, distinct set and reference.
#' @param len_adj,len_d,len_s Their lengths in bases.
#' @return List with step, pi_n, pi_shadow_n and the three normalized counts.
#' @export
edge_step <- function(adj, d, s, len_adj, len_d, len_s) {
  n <- function(x, l) if (is.na(x) || is.na(l) || l <= 0) NA_real_ else x * 100 / l
  adj_n <- n(adj, len_adj); d_n <- n(d, len_d); s_n <- n(s, len_s)
  safe <- function(a, b) if (is.na(a) || is.na(b) || (a + b) == 0) NA_real_ else a / (a + b)
  list(step        = safe(adj_n, d_n),
       pi_n        = safe(d_n, s_n),
       pi_shadow_n = safe(adj_n, s_n),
       adj_n = adj_n, d_n = d_n, s_n = s_n)
}


#' Threshold the edge step against calibration controls
#'
#' Thresholds are never hardcoded. Each dataset calibrates against its own
#' control bins: positives at the shared TSS of single-promoter genes (a real
#' boundary, expected step near 0) and negatives at the midpoint of a long
#' exon (pure continuation, expected step near 0.5).
#'
#' The test is ASYMMETRIC and the returned rates say by how much. On the DICE
#' activation panel, with the read-aware guard in \code{boundary_support}
#' (min_reads = 2), negatives were misread 4.0% and positives 26.2% --
#' specificity ~96%, sensitivity ~74%. Always report both rates alongside any
#' classification: \code{"boundary"} is specific and can be relied on, while
#' \code{"continuation"} is weaker and must not be used to declare an
#' individual call false.
#'
#' Do not read a high positive-misread rate as the test being useless -- check
#' the guard first. With a bare one-base floor the rate was 38.9%, because at
#' depth ~1 a single read in the flank forces a continuation call. The
#' aggregate boundary fraction is insensitive to the guard (27.6% at one base
#' vs 25.6% at five reads), so only per-side confidence is affected.
#'
#' @param pos,neg Numeric vectors of step values for the control bins.
#' @return List with threshold and the two misread rates.
#' @export
calibrate_boundary <- function(pos, neg) {
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  if (!length(pos) || !length(neg))
    return(list(threshold = NA_real_, pos_misread = NA_real_,
                neg_misread = NA_real_, n_pos = length(pos), n_neg = length(neg)))
  thr <- (stats::median(pos) + stats::median(neg)) / 2
  list(threshold   = thr,
       pos_misread = mean(pos >= thr),
       neg_misread = mean(neg <  thr),
       n_pos = length(pos), n_neg = length(neg))
}


#' Label sides as boundary or continuation
#'
#' The guard must be READ-AWARE, not a bare count floor. A single read in a
#' 100bp window contributes about \code{read_length} of base coverage, so at
#' low depth one stray read flips step from 0 to 0.5. With a floor of one base
#' the positive controls were misread 38.9% of the time; requiring the flank to
#' hold at least 2-5 reads before continuation can be called drops that to
#' 20-26% at unchanged specificity (~4%). The estimated boundary fraction is
#' insensitive to the choice (27.6% vs 25.6%), so this only affects per-side
#' confidence, not the aggregate.
#'
#' @param step Numeric vector of step values.
#' @param d_n Normalized distinct-set counts, for the coverage guard.
#' @param cal Result of \code{calibrate_boundary}.
#' @param adj_n Normalized flanking-bin counts. When supplied, a flank holding
#'   fewer than \code{min_reads} reads is treated as empty, i.e. a boundary,
#'   rather than letting one read force a continuation call.
#' @param min_reads Read floor for both D and the flank.
#' @param read_length Bases of coverage one read contributes.
#' @return Character vector.
#' @export
boundary_support <- function(step, d_n, cal, adj_n = NULL, min_reads = 2,
                             read_length = 100) {
  floor_n <- min_reads * read_length
  out <- rep(NA_character_, length(step))
  out[is.na(step)] <- "no_coverage"
  low <- !is.na(d_n) & d_n < floor_n & is.na(out)
  out[low] <- "low_coverage"
  ok <- is.na(out)
  if (is.na(cal$threshold)) {
    out[ok] <- "uncalibrated"
    return(out)
  }
  lab <- ifelse(step[ok] < cal$threshold, "boundary", "continuation")
  if (!is.null(adj_n)) {
    ## a flank below the read floor is empty, not continuing
    thin <- !is.na(adj_n[ok]) & adj_n[ok] < floor_n
    lab[thin] <- "boundary"
  }
  out[ok] <- lab
  out
}

## ---------------------------------------------------------------------------
## Graph-native geometry.
##
## The splice graph already states everything the diagnostic needs, so none of
## this has to be rederived from the flattened GFF:
##
##   vertices        `position` -- a BOUNDARY in increasing-coordinate space
##   ex_part edges   `dexseq_fragment` plus a boolean per transcript; their own
##                   from_pos/to_pos are NA, coordinates come from the endpoints
##   R / L edges     one per distinct TSS / TTS
##
## Conventions, verified against the GFF on a plus-strand gene (FOS,
## ENSG00000170345) and a minus-strand one (LCP2, ENSG00000043462):
##
##   part extent   [min(pos_from, pos_to), max(pos_from, pos_to) - 1]
##   TSS           pos - (strand == "-")
##   TTS           pos - (strand == "+")
##
## Transcript membership is NOT needed. An R edge points at its route's FIRST
## node, so any exonic part 5' of that node is by definition not on that path --
## there is nothing to check. And for the flank it is enough to know an exonic
## part EXISTS beyond the boundary: if one does it belongs to another
## transcript, which is exactly the confound, and whose it is does not change
## the answer.
##
## The terminal offsets are NOT symmetric by accident: `position` marks the
## lower bound of a feature, so a terminus that is the feature's UPPER bound is
## stored one past it. A one-base error here puts the flanking bin inside the
## exon, so the tests pin both strands.


#' Exonic parts from the splice graph
#'
#' Replaces reading the flattened GFF. Coordinates come from the endpoint
#' vertices because `ex_part` edges carry NA in their own from_pos/to_pos.
#'
#' @param g An igraph splice graph read from a per-gene \code{.graphml}.
#' @return Data frame with part (integer), start, end -- 1-based inclusive.
#' @export
graph_exonic_parts <- function(g) {
  ea  <- igraph::edge_attr(g)
  if (is.null(ea$ex_or_in) || is.null(ea$dexseq_fragment))
    return(data.frame(part = integer(0), start = numeric(0), end = numeric(0)))
  el  <- igraph::as_edgelist(g, names = FALSE)
  pos <- suppressWarnings(as.numeric(igraph::vertex_attr(g, "position")))
  k   <- which(ea$ex_or_in == "ex_part")
  if (!length(k))
    return(data.frame(part = integer(0), start = numeric(0), end = numeric(0)))
  a <- pos[el[k, 1L]]; b <- pos[el[k, 2L]]
  out <- data.frame(part  = suppressWarnings(as.integer(ea$dexseq_fragment[k])),
                    start = pmin(a, b),
                    end   = pmax(a, b) - 1L)
  out <- out[!is.na(out$part) & !is.na(out$start) & !is.na(out$end), , drop = FALSE]
  out[order(out$start), , drop = FALSE]
}


#' Gene strand from the graph
#'
#' \code{map_DEXSeq_from_gff()} records strand as a GRAPH attribute at build
#' time (\code{graph_utils.R}, \code{g$strand <- strand}), alongside
#' \code{gene}. Read it rather than inferring anything.
#'
#' The fallback below exists only for a graph built before that attribute was
#' written: transcription runs R to L, so R-edge positions sit below L-edge
#' positions on the plus strand and above them on the minus. Verified against
#' the GFF on 8 genes of both orientations, but the stored attribute is
#' authoritative.
#'
#' CHROMOSOME is not stored on the graph and must still come from the
#' annotation (see \code{parse_gff_chr_map}).
#'
#' @param g An igraph splice graph.
#' @return "+" or "-", or NA if neither the attribute nor the fallback resolves.
#' @export
graph_strand <- function(g) {
  st <- tryCatch(igraph::graph_attr(g, "strand"), error = function(e) NULL)
  if (!is.null(st) && length(st) && !is.na(st[1]) && st[1] %in% c("+", "-"))
    return(as.character(st[1]))
  ea <- igraph::edge_attr(g)
  if (is.null(ea$ex_or_in)) return(NA_character_)
  el  <- igraph::as_edgelist(g, names = FALSE)
  pos <- suppressWarnings(as.numeric(igraph::vertex_attr(g, "position")))
  r <- which(ea$ex_or_in == "R"); l <- which(ea$ex_or_in == "L")
  if (!length(r) || !length(l)) return(NA_character_)
  rp <- pos[el[r, 2L]]; lp <- pos[el[l, 1L]]
  rp <- rp[!is.na(rp)]; lp <- lp[!is.na(lp)]
  if (!length(rp) || !length(lp)) return(NA_character_)
  if (stats::median(rp) < stats::median(lp)) "+" else "-"
}


#' Chromosome from the graph
#'
#' \code{map_DEXSeq_from_gff()} records \code{chrom} as a graph attribute at
#' build time, alongside \code{strand} and \code{gene}. Graphs built before
#' that was added do not carry it, so callers should fall back to
#' \code{parse_gff_chr_map()} when this returns NA and regenerate when
#' convenient.
#'
#' @param g An igraph splice graph.
#' @return Chromosome name, or NA for a graph built without it.
#' @export
graph_chrom <- function(g) {
  v <- tryCatch(igraph::graph_attr(g, "chrom"), error = function(e) NULL)
  if (is.null(v) || !length(v) || is.na(v[1]) || !nzchar(as.character(v[1])))
    return(NA_character_)
  as.character(v[1])
}


#' Distinct TSS or TTS positions from the graph's R / L edges
#'
#' These are the graph's own statement of where transcription starts or stops,
#' so they replace inferring a boundary by walking exonic parts.
#'
#' @param g An igraph splice graph.
#' @param kind "TSS" (R edges) or "TTS" (L edges).
#' @param strand "+" or "-".
#' @return Data frame with sg_id and boundary, the terminal base itself.
#' @export
graph_terminal_positions <- function(g, kind = c("TSS", "TTS"), strand) {
  kind <- match.arg(kind)
  stopifnot(strand %in% c("+", "-"))
  ea  <- igraph::edge_attr(g)
  if (is.null(ea$ex_or_in))
    return(data.frame(sg_id = character(0), boundary = numeric(0)))
  el  <- igraph::as_edgelist(g, names = FALSE)
  pos <- suppressWarnings(as.numeric(igraph::vertex_attr(g, "position")))
  sg  <- as.character(igraph::vertex_attr(g, "sg_id"))
  if (kind == "TSS") {
    k <- which(ea$ex_or_in == "R"); v <- if (length(k)) el[k, 2L] else integer(0)
    off <- if (strand == "-") 1L else 0L
  } else {
    k <- which(ea$ex_or_in == "L"); v <- if (length(k)) el[k, 1L] else integer(0)
    off <- if (strand == "+") 1L else 0L
  }
  if (!length(v)) return(data.frame(sg_id = character(0), boundary = numeric(0)))
  d <- data.frame(sg_id = sg[v], boundary = pos[v] - off)
  d[!is.na(d$boundary), , drop = FALSE]
}


#' Free end of a bipartition side, read from the graph
#'
#' For each of the side's routes, take the terminal node and look up its R (or
#' L) edge boundary; the side's free end is the outermost of those. This is the
#' graph stating the answer directly, rather than
#' reconstructing it by walking exonic-part contiguity in coordinate space.
#'
#' @param g An igraph splice graph.
#' @param path_field The side's \code{path1} / \code{path2} string, routes
#'   comma-separated and nodes "-" separated.
#' @param kind "TSS" or "TTS".
#' @param strand "+" or "-".
#' @return Single boundary coordinate, or NA if no route terminus resolves.
#' @export
graph_side_boundary <- function(g, path_field, kind = c("TSS", "TTS"), strand) {
  kind <- match.arg(kind)
  term <- graph_terminal_positions(g, kind, strand)
  if (!nrow(term)) return(NA_real_)
  routes <- strsplit(as.character(path_field), ",")[[1]]
  ends <- vapply(routes, function(rt) {
    n <- trimws(strsplit(trimws(rt), "-")[[1]])
    n <- n[nzchar(n)]
    if (!length(n)) return(NA_real_)
    ## the node adjacent to the virtual terminal: first for TSS, last for TTS
    nd <- if (kind == "TSS") n[if (n[1] == "R" && length(n) > 1L) 2L else 1L]
          else               n[if (utils::tail(n, 1) == "L" && length(n) > 1L)
                                 length(n) - 1L else length(n)]
    b <- term$boundary[term$sg_id == nd]
    if (!length(b)) NA_real_ else b[1L]
  }, numeric(1))
  ends <- ends[!is.na(ends)]
  if (!length(ends)) return(NA_real_)
  ## outermost: away from the transcript body
  if ((kind == "TSS") == (strand == "+")) min(ends) else max(ends)
}
