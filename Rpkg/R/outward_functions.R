## Outward units for the TERMINAL-event boundary arm.
##
## THE QUESTION THESE SERVE. A GrASE terminal-event call says a side's
## distinct set differs in usage between conditions. It does not say where
## transcription begins or ends. To ask that, walk outward from the side's
## annotated termini and test, with the same shared reference, whether signal
## stops. These functions supply the geometry for that walk.
##
## WHY A TILING. Exonic parts of a gene are disjoint, so the gene span minus
## its exonic parts gives intronic parts that complete a DISJOINT TILING: every
## base of the gene belongs to exactly one unit. That is what makes "the unit
## immediately outward of this part" well defined at every step, with no gaps
## to fall through and no overlapping units to choose between.
##
## WHY PER ROUTE, NOT PER SIDE. A side has one terminus per ROUTE, not one per
## side. Taking a single position per side -- as a previous implementation did
## -- lands on the outermost terminus only about half the time, so the walk has
## to begin from every terminus and rank them.
##
## Part ids are compared as zero-padded three-digit strings throughout, because
## that is how DEXSeq writes `exonic_part_number` and how the intronic bin ids
## built here embed them.


## Zero-padded three-digit part id. Idempotent: "001" and 1L both give "001".
.part_id <- function(x) {
  x <- suppressWarnings(as.integer(x))
  ifelse(is.na(x), NA_character_, sprintf("%03d", x))
}


#' Does walking outward mean decreasing coordinate?
#'
#' A terminal bubble is anchored at one of the graph's two terminals: the root
#' `R`, where transcription starts, or the leaf `L`, where it ends. A side is
#' labelled `"5prime"` or `"3prime"` accordingly. Combined with gene strand this
#' fixes which direction leads away from the transcript body: outward is the
#' direction of decreasing coordinate for a `"5prime"` side on the plus strand
#' and for a `"3prime"` side on the minus strand.
#'
#' Trivial, but it is the predicate every other step depends on, and getting it
#' backwards inverts every verdict silently rather than erroring.
#'
#' @param kind Character, "5prime" (the side is anchored at the graph root
#'   `R`, so its free end is the transcription start) or "3prime" (anchored at
#'   the leaf `L`). Vectorised.
#' @param strand Character, "+" or "-". Vectorised, recycled against `kind`.
#' @return Logical vector: `TRUE` where outward means decreasing coordinate.
#' @examples
#' outward_decreasing(c("5prime", "5prime", "3prime", "3prime"), c("+", "-", "+", "-"))
#' @export
outward_decreasing <- function(kind, strand) {
  kind <- as.character(kind); strand <- as.character(strand)
  n <- max(length(kind), length(strand))
  if (!length(kind) || !length(strand)) return(logical(0))
  kind <- rep_len(kind, n); strand <- rep_len(strand, n)
  bad <- !kind %in% c("5prime", "3prime") | !strand %in% c("+", "-")
  if (any(bad))
    stop("kind must be \"5prime\" or \"3prime\" and strand + or -; bad entries at ",
         paste(which(bad), collapse = ", "), call. = FALSE)
  (kind == "5prime" & strand == "+") | (kind == "3prime" & strand == "-")
}


#' Terminus exonic parts of a bipartition side, one per route
#'
#' The terminus of a route is its first exonic part after the graph root (a
#' `"5prime"` side) or its last before the leaf (a `"3prime"` side). A path field
#' lists one route per comma-separated element, vertices within a route
#' separated by "-", so the terminal vertex is the second from the beginning
#' for `"5prime"` and the second from the end for `"3prime"`. That vertex is
#' then resolved to an exonic part through the `ex_part` edge sourced
#' (`"5prime"`) or sunk (`"3prime"`) there.
#'
#' Returns EVERY route's terminus, deduplicated, not a single position for the
#' side. Which terminus carries the differential effect is the question the walk
#' exists to answer, so it cannot be collapsed here.
#'
#' @param path_field Character scalar: the side's `path1` or `path2` value,
#'   e.g. `"1-2-5-9,1-3-6-9"`. `NA` or `""` give an empty result.
#' @param kind Character scalar, "5prime" or "3prime": which graph terminal the
#'   side is anchored at, the root `R` or the leaf `L`.
#' @param from_part Named character vector mapping an `ex_part` edge's SOURCE
#'   vertex to its part number. Used for a `"5prime"` side.
#' @param to_part Named character vector mapping an `ex_part` edge's SINK vertex
#'   to its part number. Used for a `"3prime"` side.
#' @return Character vector of zero-padded part ids, possibly empty. Routes
#'   whose terminal vertex carries no `ex_part` edge are dropped.
#' @export
side_termini <- function(path_field, kind, from_part, to_part) {
  kind <- match.arg(as.character(kind), c("5prime", "3prime"))
  if (length(path_field) != 1L)
    stop("path_field must be a single value", call. = FALSE)
  if (is.na(path_field) || !nzchar(trimws(path_field))) return(character(0))
  routes <- trimws(unlist(strsplit(as.character(path_field), ",", fixed = TRUE)))
  routes <- routes[nzchar(routes)]
  if (!length(routes)) return(character(0))
  verts <- vapply(strsplit(routes, "-", fixed = TRUE), function(v) {
    ## a route needs at least a terminal node and the root or leaf beside it
    if (length(v) < 2L) return(NA_character_)
    if (kind == "5prime") v[2L] else v[length(v) - 1L]
  }, character(1))
  verts <- verts[!is.na(verts)]
  if (!length(verts)) return(character(0))
  map <- if (kind == "5prime") from_part else to_part
  parts <- unname(map[verts])
  parts <- parts[!is.na(parts) & nzchar(parts)]
  unique(.part_id(parts))
}


#' The tiling unit immediately outward of an exonic part
#'
#' Given one gene's exonic parts and its intronic bins, return the unit that
#' abuts `part` on its outward side. Three cases, in order:
#'
#' * the neighbouring exonic part, when the tiling is contiguous there (the
#'   neighbour ends exactly one base before this part starts, or starts one
#'   base after it ends);
#' * otherwise the intronic bin lying between the two parts, taken from the end
#'   that faces `part`, so a long intron split into two end bins contributes the
#'   bin actually adjacent to the walk;
#' * otherwise nothing -- `part` faces the gene edge, or the gap between it and
#'   its neighbour produced no bin.
#'
#' The third case is not a failure and must not be read as one. A missing unit
#' means the outward space was never measured, which is different from having
#' measured it and found no signal.
#'
#' @param part Part ids to walk outward from. Vectorised.
#' @param decreasing Logical: is outward the direction of decreasing
#'   coordinate? See [outward_decreasing()]. Recycled against `part`.
#' @param geom Data frame of ONE gene's exonic parts, with columns `part`,
#'   `start`, `end`. Order does not matter; it is sorted here.
#' @param intron_bins Data frame of ONE gene's intronic bins, with columns
#'   `bin_id`, `abut_part_lo`, `abut_part_hi` and `touches` ("both", "lo" or
#'   "hi"), as returned by [intronic_parts()]. May be `NULL` or empty.
#' @return Data frame with one row per `part` and columns `unit_type`
#'   ("exon", "intron" or `NA`) and `unit_id` (a part id, a bin id, or `NA`).
#' @export
outward_unit <- function(part, decreasing, geom, intron_bins = NULL) {
  if (!all(c("part", "start", "end") %in% names(geom)))
    stop("geom needs columns part, start, end", call. = FALSE)
  n <- length(part)
  out <- data.frame(unit_type = rep(NA_character_, n),
                    unit_id   = rep(NA_character_, n),
                    stringsAsFactors = FALSE)
  if (!n) return(out)
  decreasing <- rep_len(as.logical(decreasing), n)

  g <- geom[order(as.numeric(geom$start), as.numeric(geom$end)), , drop = FALSE]
  gp <- .part_id(g$part); gs <- as.numeric(g$start); ge <- as.numeric(g$end)
  m <- length(gp)

  has_bins <- !is.null(intron_bins) && nrow(intron_bins) > 0L
  if (has_bins) {
    if (!all(c("bin_id", "abut_part_lo", "abut_part_hi", "touches") %in% names(intron_bins)))
      stop("intron_bins needs columns bin_id, abut_part_lo, abut_part_hi, touches",
           call. = FALSE)
    blo <- .part_id(intron_bins$abut_part_lo)
    bhi <- .part_id(intron_bins$abut_part_hi)
    btc <- as.character(intron_bins$touches)
    bid <- as.character(intron_bins$bin_id)
  }

  ## the bin between parts lo and hi, taken from the end facing the walk:
  ## walking down into the intron we arrive at its HIGH end, and vice versa.
  pick_bin <- function(lo, hi, dec) {
    if (!has_bins) return(NA_character_)
    want <- if (dec) c("both", "hi") else c("both", "lo")
    k <- which(blo == lo & bhi == hi & btc %in% want)
    if (length(k)) bid[k[1L]] else NA_character_
  }

  ## THE INTERGENIC FLANK at a gene edge, when the manifest carries one.
  ## A flank abuts a SINGLE part, so the unused side is NA: IG_before has
  ## abut_part_hi = the gene's lowest part and lo = NA, IG_after has
  ## abut_part_lo = the highest part and hi = NA. Matching on `touches` FIRST
  ## keeps an NA abutting part from colliding with an ordinary bin. Returns NA
  ## when the manifest has no flanks (intergenic_width = 0, the default) or
  ## where mk_intergenic() suppressed one at a chromosome start, so behaviour
  ## without flanks is byte-identical to before.
  pick_ig <- function(p, dec) {
    if (!has_bins) return(NA_character_)
    k <- if (dec) which(btc == "intergenic_before" & bhi == p)
         else      which(btc == "intergenic_after" & blo == p)
    if (length(k)) bid[k[1L]] else NA_character_
  }

  want <- .part_id(part)
  for (i in seq_len(n)) {
    j <- match(want[i], gp)
    if (is.na(j)) next
    if (decreasing[i]) {
      if (j == 1L) {                          # gene edge: flank, if built
        b <- pick_ig(gp[j], TRUE)
        if (!is.na(b)) { out$unit_type[i] <- "intron"; out$unit_id[i] <- b }
        next
      }
      if (ge[j - 1L] == gs[j] - 1L) {         # contiguous tiling
        out$unit_type[i] <- "exon"; out$unit_id[i] <- gp[j - 1L]; next
      }
      b <- pick_bin(gp[j - 1L], gp[j], TRUE)
    } else {
      if (j == m) {                           # gene edge: flank, if built
        b <- pick_ig(gp[j], FALSE)
        if (!is.na(b)) { out$unit_type[i] <- "intron"; out$unit_id[i] <- b }
        next
      }
      if (gs[j + 1L] == ge[j] + 1L) {
        out$unit_type[i] <- "exon"; out$unit_id[i] <- gp[j + 1L]; next
      }
      b <- pick_bin(gp[j], gp[j + 1L], FALSE)
    }
    if (!is.na(b)) { out$unit_type[i] <- "intron"; out$unit_id[i] <- b }
  }
  out
}


#' Intronic parts of a gene, as quantification units
#'
#' The intronic analogue of a DEXSeq exonic part: a maximal region of the gene
#' that is exonic in no isoform. Taken as the gene's exonic parts subtracted
#' from its span, so exonic and intronic parts together tile the gene without
#' gaps or overlaps.
#'
#' A short intron becomes ONE bin spanning it. A long one becomes TWO bins, at
#' its low and high ends, because the interior of a long intron says nothing
#' about whether transcription continued past a boundary at either end, and a
#' single bin would dilute the signal across it. `max_whole >= 2 * end_width`
#' is required so those two bins cannot overlap and double-count reads.
#'
#' @param start,end,part Exonic part coordinates and part numbers for ONE gene,
#'   in any order. `part` is used only to name the bins.
#' @param max_whole Introns up to this length give one bin spanning the whole
#'   intron; longer ones give two end bins.
#' @param end_width Width of each end bin for a long intron.
#' @return Data frame of bins with columns `bin_id`, `intron_id`, `touches`
#'   ("both" for a whole-intron bin, else "lo" or "hi"), `start`, `end`,
#'   `width`, `intron_start`, `intron_end`, `intron_len`, `abut_part_lo` and
#'   `abut_part_hi`. Zero rows when the gene yields no intronic part, with
#'   attribute `"skip"` giving the reason: `"single"` (fewer than two exonic
#'   parts), `"overlap"` (exonic parts overlap, so the tiling premise fails and
#'   the gene must be excluded rather than silently producing negative-width
#'   bins) or `"no_intron"` (all parts abut).
#' @export
intronic_parts <- function(start, end, part, max_whole = 2000L, end_width = 300L,
                           intergenic_width = 0L) {
  ## intergenic_width: flank width for the gene-edge bins, 0 = none (default).
  ## USE 300L, matching end_width, if these are built. An intergenic gap is
  ## unbounded exactly as a long intron is, so the reason a long intron gets an
  ## end_width bin at each end rather than one whole bin applies unchanged; and
  ## of the existing 300 bp intronic end bins only ~22-25% clear the read floor,
  ## so a narrower flank would lose about half the judgeable cases. Settled
  ## 2026-10-21. Not yet exposed by scripts/intronparts.R.
  max_whole <- as.integer(max_whole); end_width <- as.integer(end_width)
  if (is.na(max_whole) || is.na(end_width) || end_width < 1L)
    stop("max_whole and end_width must be positive integers", call. = FALSE)
  if (max_whole < 2L * end_width)
    stop("max_whole (", max_whole, ") must be >= 2 * end_width (",
         2L * end_width, "): otherwise the two end bins overlap and ",
         "double-count reads", call. = FALSE)
  if (length(start) != length(end) || length(start) != length(part))
    stop("start, end and part must be the same length", call. = FALSE)

  empty <- function(reason) {
    d <- data.frame(bin_id = character(0), intron_id = character(0),
                    touches = character(0), start = numeric(0), end = numeric(0),
                    width = numeric(0), intron_start = numeric(0),
                    intron_end = numeric(0), intron_len = numeric(0),
                    abut_part_lo = character(0), abut_part_hi = character(0),
                    stringsAsFactors = FALSE)
    attr(d, "skip") <- reason
    d
  }
  intergenic_width <- as.integer(intergenic_width)
  if (is.na(intergenic_width) || intergenic_width < 0L)
    stop("intergenic_width must be a non-negative integer", call. = FALSE)

  ## INTERGENIC GAPS. Beyond the gene's outermost exonic part there is no
  ## annotated feature, but transcription does not stop at an annotation --
  ## read-through past the 3' end is exactly what a gene-edge terminus cannot
  ## currently be tested against. The intergenic gap is unbounded, so it
  ## contributes ONE end bin of `end_width`, the same rule a long intron gets.
  ##
  ## "Intergenic" names the gap's position relative to THIS gene, not a
  ## guarantee that no other gene is there: a neighbouring gene can overlap it,
  ## and those cases must be marked unusable by the caller rather than counted,
  ## because dexseq_count resolves ambiguity at gene level and would drop the
  ## reads for both genes.
  ## An intergenic gap abuts a single part, so the unused side of the
  ## manifest's abut_part_lo/hi pair is NA.
  ##
  ## NAMING, deliberately NOT "lo"/"hi". For an intron, `touches` says which
  ## END OF THE GAP the bin occupies, so "hi" means the bin sits at the gap's
  ## high-coordinate end and therefore abuts abut_part_hi. An intergenic gap
  ## below the gene occupies its own gap's HIGH end by that rule, which would
  ## make the natural-language reading ("the low one") the exact opposite of
  ## the encoded one. So these use "intergenic_before"/"intergenic_after",
  ## naming the side of the GENE they sit on in COORDINATE order -- before is
  ## the lower coordinate, after the higher, whatever the strand -- so no
  ## reader has to resolve which convention is in play. They are deliberately
  ## not 5prime/3prime: that would be strand-dependent, and this function
  ## receives no strand. abut_part_lo/hi are unchanged: they keep
  ## their coordinate meaning, with the absent side NA.
  ##
  ## These are emitted even for a gene with no introns at all -- a single-exon
  ## gene has two intergenic gaps and zero intronic parts -- so they are built
  ## BEFORE the no-intron returns rather than appended to them.
  mk_intergenic <- function(s, e, p) {
    if (intergenic_width < 1L) return(NULL)
    lo_s <- min(s) - intergenic_width; lo_e <- min(s) - 1
    hi_s <- max(e) + 1;           hi_e <- max(e) + intergenic_width
    lo_p <- p[which.min(s)];      hi_p <- p[which.max(e)]
    keep <- lo_s >= 1
    data.frame(
      bin_id = c(if (keep) sprintf("IG_before_%s", lo_p), sprintf("IG_after_%s", hi_p)),
      intron_id = c(if (keep) sprintf("IG_before_%s", lo_p), sprintf("IG_after_%s", hi_p)),
      touches = c(if (keep) "intergenic_before", "intergenic_after"),
      start = c(if (keep) lo_s, hi_s), end = c(if (keep) lo_e, hi_e),
      width = c(if (keep) lo_e - lo_s + 1, hi_e - hi_s + 1),
      intron_start = NA_real_, intron_end = NA_real_, intron_len = NA_real_,
      abut_part_lo = c(if (keep) NA_character_, hi_p),
      abut_part_hi = c(if (keep) lo_p, NA_character_),
      stringsAsFactors = FALSE)
  }
  with_intergenic <- function(d, ig, reason) {
    out <- if (is.null(ig)) d else rbind(d, ig)
    out <- out[order(out$start, out$end), , drop = FALSE]
    rownames(out) <- NULL
    attr(out, "skip") <- reason
    out
  }

  if (length(start) < 1L) return(empty("single"))
  o <- order(as.numeric(start), as.numeric(end))
  s <- as.numeric(start)[o]; e <- as.numeric(end)[o]; p <- .part_id(part)[o]
  ## the tiling depends on this; a violation would silently produce
  ## negative-width or overlapping intronic parts. A gene that fails it gets
  ## NO intergenic gaps either -- its coordinates cannot be trusted.
  if (length(s) > 1L && any(s[-1L] <= e[-length(e)])) return(empty("overlap"))
  ig <- mk_intergenic(s, e, p)
  if (length(start) < 2L) return(with_intergenic(empty("single"), ig, "single"))

  gs <- e[-length(e)] + 1                    # intron start = previous part end + 1
  ge <- s[-1L] - 1                           # intron end   = next part start - 1
  has <- ge >= gs                            # abutting parts leave no intron
  if (!any(has)) return(with_intergenic(empty("no_intron"), ig, "no_intron"))

  lo <- p[-length(p)][has]; hi <- p[-1L][has]
  gs <- gs[has]; ge <- ge[has]
  ilen <- ge - gs + 1
  whole <- ilen <= max_whole

  bs <- c(gs[whole], gs[!whole],                  ge[!whole] - end_width + 1)
  be <- c(ge[whole], gs[!whole] + end_width - 1,  ge[!whole])
  lp <- c(lo[whole], lo[!whole], lo[!whole])
  hp <- c(hi[whole], hi[!whole], hi[!whole])
  il <- c(ilen[whole], ilen[!whole], ilen[!whole])
  istart <- c(gs[whole], gs[!whole], gs[!whole])
  iend   <- c(ge[whole], ge[!whole], ge[!whole])
  touches <- c(rep("both", sum(whole)), rep("lo", sum(!whole)), rep("hi", sum(!whole)))
  sfx <- ifelse(touches == "both", "", ifelse(touches == "lo", ".lo", ".hi"))

  d <- data.frame(
    bin_id = sprintf("I%s_%s%s", lp, hp, sfx),
    intron_id = sprintf("I%s_%s", lp, hp), touches = touches,
    start = bs, end = be, width = be - bs + 1,
    intron_start = istart, intron_end = iend, intron_len = il,
    abut_part_lo = lp, abut_part_hi = hp, stringsAsFactors = FALSE)
  with_intergenic(d, ig, NA_character_)
}
