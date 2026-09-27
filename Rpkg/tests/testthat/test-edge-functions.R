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

test_that("extend_contiguous_run walks through a contiguous same-transcript part", {
  ## distinct set is E2; E3 is contiguous (300 == 299+1) and shares tB
  expect_equal(extend_contiguous_run(mk_parts(), mk_tx(), 2L, 1L), 349L)
})

test_that("extend_contiguous_run stops at a gap", {
  ## walking down from E2: E1 ends at 199, E2 starts at 200 -> contiguous,
  ## and E1/E2 share tA, so it extends to 100
  expect_equal(extend_contiguous_run(mk_parts(), mk_tx(), 2L, -1L), 100L)
})

test_that("extend_contiguous_run does not cross a real intron", {
  ## E4 starts at 500; nothing is contiguous with it on either side
  expect_equal(extend_contiguous_run(mk_parts(), mk_tx(), 4L, 1L), 599L)
  expect_equal(extend_contiguous_run(mk_parts(), mk_tx(), 4L, -1L), 500L)
})

test_that("extend_contiguous_run will not walk into a part sharing no transcript", {
  tx <- mk_tx(); tx[["3"]] <- "tZ"   # E3 no longer shares with E2
  expect_equal(extend_contiguous_run(mk_parts(), tx, 2L, 1L), 299L)
})

test_that("extend_contiguous_run returns NA for an empty distinct set", {
  expect_true(is.na(extend_contiguous_run(mk_parts(), mk_tx(), integer(0), 1L)))
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
