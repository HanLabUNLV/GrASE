## Tests for the TSS/TTS edge-coverage diagnostic.

## Five contiguous-ish exonic parts on one gene.
##   E1 100-199   E2 200-299   E3 300-349 (contiguous with E2)
##   E4 500-599   E5 700-799
mk_parts <- function() list("1" = c(100, 199), "2" = c(200, 299),
                            "3" = c(300, 349), "4" = c(500, 599),
                            "5" = c(700, 799))
mk_tx <- function() list("1" = c("tA"), "2" = c("tA", "tB"),
                         "3" = c("tB"), "4" = c("tA", "tB"),
                         "5" = c("tA", "tB"))

test_that("free_end_direction covers all four TSS/TTS x strand cases", {
  expect_equal(free_end_direction("TSS", "+"), -1L)
  expect_equal(free_end_direction("TSS", "-"), 1L)
  expect_equal(free_end_direction("TTS", "+"), 1L)
  expect_equal(free_end_direction("TTS", "-"), -1L)
})

test_that("free_end_direction rejects bad input", {
  expect_error(free_end_direction("XXX", "+"))
  expect_error(free_end_direction("TSS", "?"))
})

test_that("flank_width measures to the next part, exclusive", {
  ## from 349 upward the next part starts at 500 -> 150 bases free
  expect_equal(flank_width(mk_parts(), 349L, 1L), 150L)
})

test_that("flank_width is large when the boundary faces the gene edge", {
  expect_gt(flank_width(mk_parts(), 100L, -1L), 1e6 - 1)
})

test_that("flank_bin returns a full-width bin when there is room", {
  b <- flank_bin(349L, 1L, width = 100L, avail = 150L)
  expect_equal(b$tier, "full")
  expect_equal(b$start, 350L); expect_equal(b$end, 449L); expect_equal(b$width, 100L)
})

test_that("flank_bin clips when the flank is short", {
  b <- flank_bin(349L, 1L, width = 100L, avail = 70L)
  expect_equal(b$tier, "clipped"); expect_equal(b$width, 70L)
})

test_that("flank_bin refuses a flank below the minimum", {
  b <- flank_bin(349L, 1L, width = 100L, avail = 20L, min_width = 50L)
  expect_equal(b$tier, "unscorable"); expect_true(is.na(b$start))
})

test_that("flank_bin runs away from the transcript in both directions", {
  up <- flank_bin(200L, -1L, width = 50L, avail = 50L)
  expect_equal(c(up$start, up$end), c(150L, 199L))
  dn <- flank_bin(200L, 1L, width = 50L, avail = 50L)
  expect_equal(c(dn$start, dn$end), c(201L, 250L))
})

test_that("edge_step is 0 when the flank is empty", {
  r <- edge_step(adj = 0, d = 500, s = 400, 100, 200, 100)
  expect_equal(r$step, 0)
})

test_that("edge_step is 0.5 for perfect continuation at equal density", {
  ## same per-base depth inside and outside
  r <- edge_step(adj = 100, d = 200, s = 50, len_adj = 100, len_d = 200, len_s = 100)
  expect_equal(r$step, 0.5)
})

test_that("edge_step normalizes by length, so raw counts alone do not decide", {
  ## adj has half the raw count of d but the same per-base depth
  r <- edge_step(adj = 100, d = 200, s = 10, len_adj = 100, len_d = 200, len_s = 100)
  expect_equal(r$step, 0.5)
})

test_that("edge_step returns NA rather than dividing by zero", {
  r <- edge_step(adj = 0, d = 0, s = 10, 100, 100, 100)
  expect_true(is.na(r$step))
})

test_that("edge_step pi_shadow_n degenerates when d_n and s_n are far apart", {
  ## this is why pi_shadow_n must not be used for scoring
  r <- edge_step(adj = 1000, d = 1000, s = 1, len_adj = 100, len_d = 100, len_s = 100)
  expect_gt(r$pi_n, 0.99)
  expect_gt(r$pi_shadow_n, 0.99)
  expect_lt(abs(r$pi_shadow_n - r$pi_n), 0.01)   # difference collapses
  expect_equal(r$step, 0.5)                      # but step still reads continuation
})

test_that("calibrate_boundary puts the threshold between the control medians", {
  cal <- calibrate_boundary(pos = c(0, 0, 0.1), neg = c(0.5, 0.5, 0.48))
  expect_gt(cal$threshold, 0); expect_lt(cal$threshold, 0.5)
})

test_that("calibrate_boundary reports both misread rates", {
  cal <- calibrate_boundary(pos = c(0, 0, 0, 0.9), neg = c(0.5, 0.5, 0.5, 0.5))
  expect_equal(cal$pos_misread, 0.25)
  expect_equal(cal$neg_misread, 0)
})

test_that("calibrate_boundary returns NA threshold without controls", {
  cal <- calibrate_boundary(pos = numeric(0), neg = c(0.5))
  expect_true(is.na(cal$threshold))
})

test_that("boundary_support labels against the calibrated threshold", {
  cal <- calibrate_boundary(pos = c(0, 0), neg = c(0.5, 0.5))
  lab <- boundary_support(step = c(0.05, 0.48), d_n = c(1e4, 1e4), cal)
  expect_equal(lab, c("boundary", "continuation"))
})

test_that("boundary_support guards against low coverage on a read scale", {
  cal <- calibrate_boundary(pos = c(0, 0), neg = c(0.5, 0.5))
  ## d_n of 100 is ONE read at read_length 100 -- below a 2-read floor
  lab <- boundary_support(step = 0.05, d_n = 100, cal, min_reads = 2)
  expect_equal(lab, "low_coverage")
})

test_that("a flank below the read floor reads as boundary, not continuation", {
  cal <- calibrate_boundary(pos = c(0, 0), neg = c(0.5, 0.5))
  ## step says continuation, but the flank holds a single read -> empty flank
  lab <- boundary_support(step = 0.5, d_n = 1e4, cal, adj_n = 100,
                          min_reads = 2)
  expect_equal(lab, "boundary")
})

test_that("a flank above the read floor is allowed to read as continuation", {
  cal <- calibrate_boundary(pos = c(0, 0), neg = c(0.5, 0.5))
  lab <- boundary_support(step = 0.5, d_n = 1e4, cal, adj_n = 1e4,
                          min_reads = 2)
  expect_equal(lab, "continuation")
})

test_that("boundary_support marks uncalibrated when controls are missing", {
  cal <- calibrate_boundary(pos = numeric(0), neg = numeric(0))
  lab <- boundary_support(step = 0.3, d_n = 1e4, cal)
  expect_equal(lab, "uncalibrated")
})

test_that("boundary_support reports no_coverage for an NA step", {
  cal <- calibrate_boundary(pos = c(0, 0), neg = c(0.5, 0.5))
  expect_equal(boundary_support(NA_real_, 1e4, cal), "no_coverage")
})

## --- graph-native geometry ---------------------------------------------------
## The splice graph states the geometry directly, so none of it needs rederiving
## from the flattened GFF. The conventions are NOT symmetric across strands and
## a one-base error puts the flanking bin inside the exon, so both are pinned.
##
## Vertex `position` is a BOUNDARY in increasing-coordinate space:
##   part extent  [min(pos_from,pos_to), max(pos_from,pos_to) - 1]
##   TSS          pos - (strand == "-")
##   TTS          pos - (strand == "+")

mk_graph <- function(strand = "+") {
  ## two exonic parts and a terminus, wired the way a real graphml is
  if (strand == "+") {
    ## part1 100-199, part2 200-299; TSS at 100, TTS at 299
    pos <- c(100, 200, 300)
    ed  <- rbind(c(1, 2), c(2, 3))          # ex_part edges
    term <- rbind(c(4, 1), c(3, 5))         # R -> v1, v3 -> L
  } else {
    ## minus strand: positions descend along transcription
    pos <- c(300, 200, 100)
    ed  <- rbind(c(1, 2), c(2, 3))
    term <- rbind(c(4, 1), c(3, 5))
  }
  g <- igraph::make_empty_graph(n = 5, directed = TRUE)
  g <- igraph::add_edges(g, c(t(ed)))
  g <- igraph::add_edges(g, c(t(term)))
  igraph::V(g)$position <- c(pos, NA, NA)
  igraph::V(g)$sg_id    <- c("1", "2", "3", "R", "L")
  igraph::E(g)$ex_or_in <- c("ex_part", "ex_part", "R", "L")
  igraph::E(g)$dexseq_fragment <- c("001", "002", NA, NA)
  g
}

test_that("graph_exonic_parts derives part extents on the plus strand", {
  p <- graph_exonic_parts(mk_graph("+"))
  expect_equal(nrow(p), 2)
  expect_equal(p$start, c(100, 200))
  expect_equal(p$end,   c(199, 299))   # max - 1
})

test_that("graph_exonic_parts is strand-agnostic", {
  p <- graph_exonic_parts(mk_graph("-"))
  expect_equal(sort(p$start), c(100, 200))
  expect_equal(sort(p$end),   c(199, 299))
})

test_that("TSS has no offset on the plus strand", {
  t <- graph_terminal_positions(mk_graph("+"), "TSS", "+")
  expect_equal(t$boundary, 100)
})

test_that("TSS is offset by one on the minus strand", {
  ## R edge points at v1, position 300; the TSS itself is 299
  t <- graph_terminal_positions(mk_graph("-"), "TSS", "-")
  expect_equal(t$boundary, 299)
})

test_that("TTS offsets are the mirror of TSS", {
  expect_equal(graph_terminal_positions(mk_graph("+"), "TTS", "+")$boundary, 299)
  expect_equal(graph_terminal_positions(mk_graph("-"), "TTS", "-")$boundary, 100)
})

test_that("graph_side_boundary takes the OUTERMOST route terminus", {
  g <- mk_graph("+")
  ## two routes starting at different nodes: 1 (pos 100) and 2 (pos 200).
  ## On the plus strand the outermost TSS is the smaller coordinate.
  g <- igraph::add_edges(g, c(4, 2))
  igraph::E(g)$ex_or_in[igraph::ecount(g)] <- "R"
  expect_equal(graph_side_boundary(g, "R-1-2-3", "TSS", "+"), 100)
  expect_equal(graph_side_boundary(g, "R-2-3",   "TSS", "+"), 200)
  expect_equal(graph_side_boundary(g, "R-1-2-3, R-2-3", "TSS", "+"), 100)
})

test_that("graph_side_boundary returns NA when no terminus resolves", {
  expect_true(is.na(graph_side_boundary(mk_graph("+"), "R-99", "TSS", "+")))
})

test_that("graph helpers tolerate a graph with no ex_part or R edges", {
  g <- igraph::make_empty_graph(n = 2, directed = TRUE)
  g <- igraph::add_edges(g, c(1, 2))
  igraph::V(g)$position <- c(1, 2); igraph::V(g)$sg_id <- c("1", "2")
  igraph::E(g)$ex_or_in <- "in"; igraph::E(g)$dexseq_fragment <- NA
  expect_equal(nrow(graph_exonic_parts(g)), 0)
  expect_equal(nrow(graph_terminal_positions(g, "TSS", "+")), 0)
})

test_that("graph_strand reads the stored graph attribute, not an inference", {
  ## map_DEXSeq_from_gff records strand as a graph attribute at build time.
  ## It must win even when the R/L geometry would suggest otherwise.
  g <- mk_graph("+")
  g <- igraph::set_graph_attr(g, "strand", "-")
  expect_equal(graph_strand(g), "-")
})

test_that("graph_strand falls back to R/L geometry when the attribute is absent", {
  expect_equal(graph_strand(mk_graph("+")), "+")
  expect_equal(graph_strand(mk_graph("-")), "-")
})

test_that("graph_chrom reads the stored chromosome attribute", {
  g <- igraph::set_graph_attr(mk_graph("+"), "chrom", "chr7")
  expect_equal(graph_chrom(g), "chr7")
})

test_that("graph_chrom returns NA for a graph built before chrom was stored", {
  ## such graphs need the parse_gff_chr_map fallback, so NA must be detectable
  expect_true(is.na(graph_chrom(mk_graph("+"))))
})

test_that("graph_chrom rejects an empty or missing attribute", {
  g <- igraph::set_graph_attr(mk_graph("+"), "chrom", "")
  expect_true(is.na(graph_chrom(g)))
})
