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


#' Outer edge of the contiguous run carrying a side's distinct set
#'
#' The distinct set under-represents an alternative terminal exon: it is an
#' intersection across the side's transcripts, so it commonly stops short of
#' the exon's free end. Walking outward through parts that are contiguous with
#' it and share a transcript with it recovers the true edge.
#'
#' This matters a great deal in practice. Measured on the DICE activation
#' panel, taking the distinct set's own edge leaves a median flank of 1 bp and
#' only 9,972 of 20,887 sides testable; extending the run gives a median flank
#' of 903 bp and 18,412 testable.
#'
#' @param parts Named list or data frame of exonic parts: integer part number
#'   to c(start, end), 1-based inclusive.
#' @param tx Named list mapping part number to a character vector of
#'   transcript IDs.
#' @param distinct_parts Integer vector of part numbers in the distinct set.
#' @param direction -1 to walk toward lower coordinates, +1 toward higher.
#' @return Single integer, the outer coordinate of the run.
#' @export
extend_contiguous_run <- function(parts, tx, distinct_parts, direction) {
  stopifnot(direction %in% c(-1L, 1L, -1, 1))
  if (!length(distinct_parts)) return(NA_integer_)
  sdtx <- unique(unlist(tx[as.character(distinct_parts)]))
  cur <- if (direction < 0)
    min(vapply(parts[as.character(distinct_parts)], `[`, numeric(1), 1)) else
    max(vapply(parts[as.character(distinct_parts)], `[`, numeric(1), 2))
  seen <- character(0)
  repeat {
    nxt <- names(parts)[vapply(parts, function(p)
      if (direction < 0) p[2] == cur - 1 else p[1] == cur + 1, logical(1))]
    nxt <- setdiff(nxt, seen)
    nxt <- nxt[vapply(nxt, function(n) any(tx[[n]] %in% sdtx), logical(1))]
    if (!length(nxt)) return(as.integer(cur))
    seen <- c(seen, nxt[1])
    cur <- if (direction < 0) parts[[nxt[1]]][1] else parts[[nxt[1]]][2]
  }
}


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
#' @param parts As for \code{extend_contiguous_run}.
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
