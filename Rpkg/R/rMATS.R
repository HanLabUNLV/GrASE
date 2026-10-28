# Graph-based rMATS mapping utilities (split-free) 
#' Compute rMATS event bubble boundaries and alternative region
#'
#' Given an rMATS event type, strand, and one row of a fromGTF data frame,
#' returns the graph positions of the bubble boundaries (b1, b2) and the
#' coordinate range of the alternative region (alt_low, alt_hi).
#' All positions follow the graph +1 convention (0-based rMATS coords + 1).
#'
#' @param etype  Event type: "SE", "RI", "A3SS", or "A5SS".
#' @param strand Gene strand: "+" or "-".
#' @param row    A single-row data frame or list with rMATS fromGTF columns.
#' @return list(b1, b2, alt_low, alt_hi) as integers, or NULL if unrecognised.
#' @export
#' @examples
#' se_row <- list(upstreamEE = 1000L, downstreamES = 2000L,
#'                exonStart_0base = 1100L, exonEnd = 1500L)
#' grase::rmats_event_boundaries("SE", "+", se_row)
rmats_event_boundaries <- function(etype, strand, row) {
  p1 <- function(x) as.integer(x) + 1L   # 0-based -> graph position

  if (etype == "SE") {
    list(b1      = p1(row$upstreamEE),
         b2      = p1(row$downstreamES),
         alt_low = p1(row$exonStart_0base),
         alt_hi  = p1(row$exonEnd))

  } else if (etype == "RI") {
    list(b1      = p1(row$upstreamEE),
         b2      = p1(row$downstreamES),
         alt_low = p1(row$upstreamEE),
         alt_hi  = p1(row$downstreamES))

  } else if (etype == "A3SS") {
    les <- p1(row$longExonStart_0base)
    lee <- p1(row$longExonEnd)
    ses <- p1(row$shortES)
    see <- p1(row$shortEE)
    fes <- p1(row$flankingES)
    fee <- p1(row$flankingEE)
    if (strand == "+") {
      list(b1 = fee, b2 = lee, alt_low = les, alt_hi = ses)
    } else {
      list(b1 = fes, b2 = les, alt_low = see, alt_hi = lee)
    }

  } else if (etype == "A5SS") {
    les <- p1(row$longExonStart_0base)
    lee <- p1(row$longExonEnd)
    ses <- p1(row$shortES)
    see <- p1(row$shortEE)
    fes <- p1(row$flankingES)
    fee <- p1(row$flankingEE)
    if (strand == "+") {
      list(b1 = les, b2 = fes, alt_low = see, alt_hi = lee)
    } else {
      list(b1 = lee, b2 = fee, alt_low = les, alt_hi = ses)
    }

  } else {
    NULL
  }
}


#' Precompute per-gene edge tables from a GrASE splice graph
#'
#' Builds a compact cache of edge data for fast repeated queries over multiple
#' rMATS events on the same gene.  The cache includes a transcript membership
#' matrix, edge-attribute positions (valid for regular ex/in edges), vertex-based
#' positions (valid for ALL edges including ex_part edges whose from_pos/to_pos
#' attributes are stored as "NA"), and DEXSeq fragment labels.
#'
#' @param g An igraph object read from a GrASE graphml file.
#' @return A named list with elements: tx_cols, tx_mat, ex_or_in, from_pos,
#'   to_pos, vx_from, vx_to, dex_frag.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' ge <- grase::precompute_gene_graph(g)
#' names(ge)
#' }
precompute_gene_graph <- function(g) {
  reserved <- c("sgedge_id", "ex_or_in", "from_pos", "to_pos", "dexseq_fragment")
  tx_cols  <- setdiff(igraph::edge_attr_names(g), reserved)

  # GraphML boolean attrs may be read as logical or character; as.logical() handles both.
  tx_mat <- matrix(0L, nrow = igraph::ecount(g), ncol = length(tx_cols))
  colnames(tx_mat) <- tx_cols
  for (tx in tx_cols) {
    v <- igraph::edge_attr(g, tx)
    tx_mat[, tx] <- as.integer(as.logical(v))
  }

  # Vertex positions (integer; NA for the root "R" vertex)
  v_pos  <- suppressWarnings(as.integer(igraph::V(g)$position))
  v_sgid <- suppressWarnings(as.integer(igraph::V(g)$sg_id))

  # Actual from/to positions for every edge via graph topology.
  # Works for ex_part edges too (whose from_pos/to_pos attributes are "NA").
  ee <- igraph::ends(g, igraph::E(g), names = FALSE)

  list(
    tx_cols    = tx_cols,
    tx_mat     = tx_mat,
    ex_or_in   = igraph::E(g)$ex_or_in,
    from_pos   = suppressWarnings(as.integer(igraph::E(g)$from_pos)),
    to_pos     = suppressWarnings(as.integer(igraph::E(g)$to_pos)),
    vx_from    = v_pos[ee[, 1L]],
    vx_to      = v_pos[ee[, 2L]],
    vx_sg_from = v_sgid[ee[, 1L]],
    vx_sg_to   = v_sgid[ee[, 2L]],
    dex_frag   = igraph::E(g)$dexseq_fragment
  )
}


#' Map one rMATS event to DEXSeq fragments via transcript set-difference
#'
#' Uses a precomputed gene graph cache (from \code{precompute_gene_graph}) to
#' identify which DEXSeq exonic-part fragments are specific to the alternative
#' isoform of an rMATS event.  Local transcripts (those incident to both bubble
#' boundary vertices b1 and b2) are split into path1 (those with exonic coverage
#' of the alternative region) and path2 (remaining local transcripts).  The
#' set-difference gives the ex_part edges TRUE in path1 and FALSE in all path2
#' transcripts.
#'
#' @param ge      Precomputed gene graph list from \code{precompute_gene_graph}.
#' @param b1      Integer graph position of the bubble entry boundary.
#' @param b2      Integer graph position of the bubble exit boundary.
#' @param alt_low Integer lower bound of the alternative region (inclusive).
#' @param alt_hi  Integer upper bound of the alternative region (inclusive).
#' @return Sorted character vector of DEXSeq fragment labels (e.g. "E001"),
#'   or \code{character(0)} if no fragments could be identified.
#' @export
#' @examples
#' \dontrun{
#' g <- igraph::read_graph("ENSG00000000003.dexseq.graphml", format = "graphml")
#' ge <- grase::precompute_gene_graph(g)
#' se_row <- list(upstreamEE = 1000L, downstreamES = 2000L,
#'                exonStart_0base = 1100L, exonEnd = 1500L)
#' bounds <- grase::rmats_event_boundaries("SE", "+", se_row)
#' frags <- grase::graph_setdiff_frags(ge,
#'   b1 = bounds$b1, b2 = bounds$b2,
#'   alt_low = bounds$alt_low, alt_hi = bounds$alt_hi
#' )
#' }
graph_setdiff_frags <- function(ge, b1, b2, alt_low, alt_hi) {
  # Edges incident to boundary vertices, by vertex position so that ex_part
  # edges count: a transcript retaining an intron touches its boundary
  # vertices only through ex_part edges.
  fp <- ge$vx_from
  tp <- ge$vx_to
  at_b1 <- which(!is.na(fp) & !is.na(tp) & (fp == b1 | tp == b1))
  at_b2 <- which(!is.na(fp) & !is.na(tp) & (fp == b2 | tp == b2))
  if (length(at_b1) == 0L || length(at_b2) == 0L) return(character(0))

  # Local transcripts: TRUE on at least one edge at b1 AND at least one at b2
  tx_at_b1   <- colSums(ge$tx_mat[at_b1, , drop = FALSE]) > 0L
  tx_at_b2   <- colSums(ge$tx_mat[at_b2, , drop = FALSE]) > 0L
  local_mask <- tx_at_b1 & tx_at_b2
  if (sum(local_mask) == 0L) return(character(0))

  # ex_part edges within the alternative region.
  # Use vertex-based positions (vx_from/vx_to) since ex_part edges store "NA"
  # in their from_pos/to_pos attributes.
  af        <- ge$vx_from
  at_pos    <- ge$vx_to
  is_expart <- ge$ex_or_in == "ex_part"

  in_alt <- which(is_expart & !is.na(af) & !is.na(at_pos) &
                    pmin(af, at_pos) >= alt_low &
                    pmax(af, at_pos) <= alt_hi)
  if (length(in_alt) == 0L) return(character(0))

  # path1: local txs with TRUE on any alt ex_part edge
  path1_mask <- local_mask &
    (colSums(ge$tx_mat[in_alt, , drop = FALSE]) > 0L)
  path2_mask <- local_mask & !path1_mask
  if (sum(path1_mask) == 0L || sum(path2_mask) == 0L) return(character(0))

  # Set-difference: alt ex_part edges TRUE in any path1 tx AND FALSE in all path2 txs
  any_p1       <- rowSums(ge$tx_mat[in_alt, path1_mask, drop = FALSE]) > 0L
  all_p2_false <- rowSums(ge$tx_mat[in_alt, path2_mask, drop = FALSE]) == 0L
  sdiff_idx    <- in_alt[any_p1 & all_p2_false]
  if (length(sdiff_idx) == 0L) return(character(0))

  frags <- ge$dex_frag[sdiff_idx]
  frags <- frags[!is.na(frags) & nchar(frags) > 0L]
  if (length(frags) == 0L) return(character(0))

  sort(unique(paste0("E", sprintf("%03d", as.integer(frags)))))
}


#' The two forms of an rMATS event as splicing-graph vertex positions
#'
#' Positions follow the graph: an exon starts at its 1-based first base and
#' ends at one past its last base, so rMATS 0-based starts and ends both map
#' by adding 1. Each form is a list of position pairs, each a graph edge that a
#' transcript route visits as consecutive vertices (lower position first).
#' rMATS builds an event from junctions, so the pairs of the inclusion form
#' need not lie on one transcript. The retained form of an intron retention is
#' instead an exon edge spanning the intron (`span`).
#'
#' @param etype One of "SE", "A3SS", "A5SS", "RI".
#' @param strand "+" or "-".
#' @param row One row of an rMATS fromGTF table.
#' @return list(incl, skip, span) or NULL for other event types.
#' @keywords internal
rmats_event_forms <- function(etype, strand, row) {
  p1 <- function(x) as.integer(x) + 1L
  if (etype == "SE") {
    b1 <- p1(row$upstreamEE); b2 <- p1(row$downstreamES)
    lo <- p1(row$exonStart_0base); hi <- p1(row$exonEnd)
    return(list(incl = list(c(b1, lo), c(lo, hi), c(hi, b2)),
                skip = list(c(b1, b2)), span = NULL))
  }
  if (etype == "RI") {
    b1 <- p1(row$upstreamEE); b2 <- p1(row$downstreamES)
    return(list(incl = NULL, skip = list(c(b1, b2)), span = c(b1, b2)))
  }
  if (!etype %in% c("A3SS", "A5SS")) return(NULL)
  les <- p1(row$longExonStart_0base); lee <- p1(row$longExonEnd)
  ses <- p1(row$shortES);             see <- p1(row$shortEE)
  fes <- p1(row$flankingES);          fee <- p1(row$flankingEE)
  # the alternative site is the long/short exon's start when the flanking
  # exon lies to its left in genomic order, else its end
  if ((etype == "A3SS") == (strand == "+"))
    list(incl = list(c(fee, les)), skip = list(c(fee, ses)), span = NULL)
  else
    list(incl = list(c(lee, fes)), skip = list(c(see, fes)), span = NULL)
}


#' Match rMATS events to the bipartitions of one gene
#'
#' Both forms of each event (`rmats_event_forms()`) are compared with the
#' transcript routes of the bipartition (`path1`, `path2`). A bipartition is a
#' candidate when every route on one side contains the skipping form, and on
#' the other side every pair of the inclusion form lies on some route and no
#' route contains the skipping form. The inclusion side may hold other versions of
#' the included sequence because rMATS counts exon-body reads, which those
#' versions share, toward inclusion; the skipping count comes only from the
#' skipping junction. Among candidates
#' the one with the narrowest source-to-sink span is taken; a tie is left
#' unmatched.
#'
#' `R` and `L` have no coordinate; they are placed beyond every vertex at the
#' gene's transcription start and end.
#'
#' @param ge Output of `precompute_gene_graph()`.
#' @param bip Data frame of one gene's bipartitions (source, sink, path1, path2).
#' @param etype One of "SE", "A3SS", "A5SS", "RI".
#' @param ev rMATS fromGTF rows of this gene and type.
#' @return Integer vector, the matched row of `bip` per event, NA if none;
#'   attribute `ncand` holds the number of candidates per event.
#' @keywords internal
match_rmats_bipartitions <- function(ge, bip, etype, ev) {
  out <- rep(NA_integer_, nrow(ev))
  ncand <- integer(nrow(ev))
  if (!nrow(bip) || !nrow(ev)) return(structure(out, ncand = ncand))
  n2p <- stats::setNames(c(ge$vx_from, ge$vx_to), as.character(c(ge$vx_sg_from, ge$vx_sg_to)))
  n2p <- n2p[!duplicated(names(n2p))]
  is_ex <- ge$ex_or_in == "ex" & !is.na(ge$vx_from) & !is.na(ge$vx_to)
  ex_key <- paste(pmin(ge$vx_from, ge$vx_to), pmax(ge$vx_from, ge$vx_to))[is_ex]
  # each side: list of routes, each a vector of positions in increasing order
  routes <- function(s) lapply(strsplit(trimws(strsplit(s, ",")[[1]]), "-"), function(v) {
    x <- unname(n2p[v[!v %in% c("R", "L")]]); x <- x[!is.na(x)]
    if (length(x) > 1L && x[1] > x[length(x)]) rev(x) else x
  })
  r1 <- lapply(bip$path1, routes); r2 <- lapply(bip$path2, routes)
  has_run <- function(x, q) grepl(paste0(",", paste(q, collapse = ","), ","),
                                  paste0(",", paste(x, collapse = ","), ","), fixed = TRUE)
  has_span <- function(x, q) {
    k <- which(utils::head(x, -1) < q[1] & x[-1] > q[2])
    any(paste(x[k], x[k + 1L]) %in% ex_key)
  }
  # per route: does it hold every pair of `form` (or the spanning exon edge)?
  on_route <- function(rs, form) vapply(rs, function(x)
    all(vapply(form, function(q) has_run(x, q), logical(1))), logical(1))
  # does the side hold every pair of the inclusion form, on any of its routes?
  side_incl <- function(rs, fm) if (!is.null(fm$span))
    any(vapply(rs, has_span, logical(1), q = fm$span)) else
    all(vapply(fm$incl, function(q) any(vapply(rs, has_run, logical(1), q = q)), logical(1)))
  sp <- unname(n2p[as.character(bip$source)]); kp <- unname(n2p[as.character(bip$sink)])
  rpos <- if (identical(as.character(ev$strand[1]), "-")) Inf else -Inf
  sp[bip$source == "R"] <- rpos; sp[bip$source == "L"] <- -rpos
  kp[bip$sink   == "R"] <- rpos; kp[bip$sink   == "L"] <- -rpos
  w <- abs(kp - sp)
  for (i in seq_len(nrow(ev))) {
    fm <- rmats_event_forms(etype, ev$strand[i], ev[i, ])
    if (is.null(fm)) next
    ok <- function(a, b) {
      length(a) > 0L && length(b) > 0L && all(on_route(b, fm$skip)) &&
        !any(on_route(a, fm$skip)) && side_incl(a, fm)
    }
    cand <- which(vapply(seq_along(r1), function(k)
      ok(r1[[k]], r2[[k]]) || ok(r2[[k]], r1[[k]]), logical(1)))
    ncand[i] <- length(cand)
    if (!length(cand)) next
    best <- cand[w[cand] == min(w[cand])]
    if (length(best) == 1L) out[i] <- best
  }
  structure(out, ncand = ncand)
}
