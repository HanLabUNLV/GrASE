## Tests for the outward units of the TSS/TTS boundary arm.

## One gene, five exonic parts. 001/002 abut exactly (no intron between them),
## the rest are separated by gaps of differing length:
##   001 100-199 | 002 200-299 | gap 300-399 | 003 400-499 |
##   gap 500-2499 (long) | 004 2500-2599 | gap 2600-2699 | 005 2700-2799
mk_geom <- function() data.frame(
  part  = c("001", "002", "003", "004", "005"),
  start = c(100, 200, 400, 2500, 2700),
  end   = c(199, 299, 499, 2599, 2799),
  stringsAsFactors = FALSE)

mk_bins <- function() intronic_parts(mk_geom()$start, mk_geom()$end,
                                     mk_geom()$part, max_whole = 1000L,
                                     end_width = 300L)


test_that("outward_decreasing covers all four kind x strand cases", {
  expect_equal(outward_decreasing(c("TSS", "TSS", "TTS", "TTS"),
                                  c("+", "-", "+", "-")),
               c(TRUE, FALSE, FALSE, TRUE))
})

test_that("outward_decreasing recycles and returns empty for empty input", {
  expect_equal(outward_decreasing(c("TSS", "TTS"), "+"), c(TRUE, FALSE))
  expect_equal(outward_decreasing(character(0), "+"), logical(0))
})

test_that("outward_decreasing rejects bad input rather than guessing", {
  expect_error(outward_decreasing("TSS", "*"), "strand")
  expect_error(outward_decreasing("start", "+"), "TSS")
})


test_that("side_termini takes the vertex after the root for a TSS", {
  fp <- c("2" = "007", "3" = "009")
  expect_equal(side_termini("1-2-5-9", "TSS", fp, c()), "007")
})

test_that("side_termini takes the vertex before the leaf for a TTS", {
  tp <- c("5" = "011", "6" = "012")
  expect_equal(side_termini("1-2-5-9", "TTS", c(), tp), "011")
})

test_that("side_termini returns one terminus PER ROUTE, deduplicated", {
  fp <- c("2" = "007", "3" = "008")
  ## two routes with different termini, and a third repeating the first
  expect_equal(sort(side_termini("1-2-5-9,1-3-6-9,1-2-7-9", "TSS", fp, c())),
               c("007", "008"))
})

test_that("side_termini pads part ids to three digits", {
  expect_equal(side_termini("1-2-9", "TSS", c("2" = "7"), c()), "007")
})

test_that("side_termini is empty for a blank, NA or too-short path", {
  fp <- c("2" = "007")
  expect_equal(side_termini(NA_character_, "TSS", fp, c()), character(0))
  expect_equal(side_termini("", "TSS", fp, c()), character(0))
  expect_equal(side_termini("4", "TSS", fp, c()), character(0))
})

test_that("side_termini drops a route whose terminal vertex has no ex_part edge", {
  fp <- c("2" = "007")
  expect_equal(side_termini("1-2-9,1-99-9", "TSS", fp, c()), "007")
})


test_that("outward_unit returns the abutting exonic part when the tiling is contiguous", {
  ## 002 starts one base after 001 ends, so walking down from 002 gives 001
  u <- outward_unit("002", TRUE, mk_geom(), mk_bins())
  expect_equal(u$unit_type, "exon")
  expect_equal(u$unit_id, "001")
  ## and walking up from 001 gives 002
  u <- outward_unit("001", FALSE, mk_geom(), mk_bins())
  expect_equal(c(u$unit_type, u$unit_id), c("exon", "002"))
})

test_that("outward_unit returns the intronic bin when there is a gap", {
  ## 300-399 lies between 002 and 003 and is short, so it is a single bin
  u <- outward_unit("003", TRUE, mk_geom(), mk_bins())
  expect_equal(u$unit_type, "intron")
  expect_equal(u$unit_id, "I002_003")
})

test_that("outward_unit takes the bin at the end FACING the walk", {
  ## 500-2499 is long, so it has two end bins. Walking DOWN from 004 arrives at
  ## the intron's high end; walking UP from 003 arrives at its low end.
  expect_equal(outward_unit("004", TRUE,  mk_geom(), mk_bins())$unit_id, "I003_004.hi")
  expect_equal(outward_unit("003", FALSE, mk_geom(), mk_bins())$unit_id, "I003_004.lo")
})

test_that("outward_unit returns NA at a gene edge, in both directions", {
  u <- outward_unit(c("001", "005"), c(TRUE, FALSE), mk_geom(), mk_bins())
  expect_true(all(is.na(u$unit_type)))
  expect_true(all(is.na(u$unit_id)))
})

test_that("outward_unit returns NA for a part the gene does not have", {
  u <- outward_unit("099", TRUE, mk_geom(), mk_bins())
  expect_true(is.na(u$unit_type))
})

test_that("outward_unit returns NA when a gap exists but no bin covers it", {
  ## the gap between 002 and 003 is real, but with no bins there is no unit --
  ## which is unmeasured space, NOT evidence that signal stopped
  u <- outward_unit("003", TRUE, mk_geom(), intron_bins = NULL)
  expect_true(is.na(u$unit_type))
})

test_that("outward_unit is vectorised over part and decreasing", {
  u <- outward_unit(c("002", "003", "001"), c(TRUE, TRUE, TRUE),
                    mk_geom(), mk_bins())
  expect_equal(u$unit_type, c("exon", "intron", NA))
  expect_equal(nrow(u), 3L)
})

test_that("outward_unit accepts unpadded part ids", {
  expect_equal(outward_unit(2L, TRUE, mk_geom(), mk_bins())$unit_id, "001")
})


test_that("intronic_parts gives one bin for a short intron and none where parts abut", {
  b <- mk_bins()
  ## 001/002 abut, so no bin between them
  expect_false("I001_002" %in% b$intron_id)
  short <- b[b$intron_id == "I002_003", ]
  expect_equal(nrow(short), 1L)
  expect_equal(short$touches, "both")
  expect_equal(c(short$start, short$end), c(300, 399))
  expect_equal(short$width, 100)
})

test_that("intronic_parts splits a long intron into two non-overlapping end bins", {
  b <- mk_bins()
  long <- b[b$intron_id == "I003_004", ]
  expect_equal(nrow(long), 2L)
  expect_equal(sort(long$touches), c("hi", "lo"))
  lo <- long[long$touches == "lo", ]; hi <- long[long$touches == "hi", ]
  expect_equal(c(lo$start, lo$end), c(500, 799))
  expect_equal(c(hi$start, hi$end), c(2200, 2499))
  expect_true(hi$start > lo$end)              # the point of max_whole >= 2*end_width
  expect_equal(c(lo$intron_len, hi$intron_len), c(2000, 2000))
})

test_that("intronic_parts refuses settings where the two end bins would overlap", {
  expect_error(intronic_parts(c(1, 100), c(10, 200), c("001", "002"),
                              max_whole = 100L, end_width = 300L),
               "double-count")
})

test_that("intronic_parts excludes a gene whose exonic parts overlap", {
  b <- intronic_parts(c(100, 150), c(200, 300), c("001", "002"))
  expect_equal(nrow(b), 0L)
  expect_equal(attr(b, "skip"), "overlap")
})

test_that("intronic_parts reports why it emitted nothing", {
  expect_equal(attr(intronic_parts(100, 199, "001"), "skip"), "single")
  expect_equal(attr(intronic_parts(c(100, 200), c(199, 299),
                                   c("001", "002")), "skip"), "no_intron")
  expect_true(is.na(attr(mk_bins(), "skip")))
})

test_that("intronic_parts is insensitive to input order", {
  g <- mk_geom()
  o <- c(3, 1, 5, 2, 4)
  expect_equal(intronic_parts(g$start[o], g$end[o], g$part[o],
                              max_whole = 1000L, end_width = 300L),
               mk_bins())
})

test_that("intronic_parts names bins after the parts they lie between", {
  b <- mk_bins()
  expect_true(all(b$bin_id == paste0(b$intron_id,
    ifelse(b$touches == "both", "", paste0(".", b$touches)))))
  expect_true(all(b$abut_part_lo < b$abut_part_hi))
})
