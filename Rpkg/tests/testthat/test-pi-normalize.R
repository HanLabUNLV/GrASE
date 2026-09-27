## Length-normalized pi.
##
## pi from the test is a RAW COUNT ratio, which is correct for the
## beta-binomial likelihood but is not a molecular proportion when the distinct
## set and reference differ in length.

test_that("pi_perbase leaves pi unchanged when lengths are equal", {
  expect_equal(pi_perbase(0.4, 100, 100), 0.4)
  expect_equal(pi_perbase(0.9, 250, 250), 0.9)
})

test_that("pi_perbase deflates pi when the distinct set is longer", {
  ## LCP2 event 26: D = E021 at 221 bp, S = E020 at 67 bp, reported pi 0.548
  v <- pi_perbase(0.548, 221, 67)
  expect_lt(v, 0.548)
  expect_equal(round(v, 3), 0.269)
})

test_that("pi_perbase inflates pi when the reference is longer", {
  expect_gt(pi_perbase(0.3, 50, 500), 0.3)
})

test_that("pi_perbase preserves the endpoints", {
  expect_equal(pi_perbase(0, 221, 67), 0)
  expect_equal(pi_perbase(1, 221, 67), 1)
})

test_that("pi_perbase is monotone in pi", {
  p <- seq(0, 1, by = 0.05)
  v <- pi_perbase(p, 221, 67)
  expect_true(all(diff(v) >= 0))
})

test_that("pi_perbase returns NA for missing or non-positive lengths", {
  expect_true(is.na(pi_perbase(0.5, NA, 67)))
  expect_true(is.na(pi_perbase(0.5, 0, 67)))
  expect_true(is.na(pi_perbase(NA_real_, 221, 67)))
})

test_that("pi_perbase is vectorized elementwise", {
  v <- pi_perbase(c(0.5, 0.5), c(100, 200), c(100, 100))
  expect_equal(v[1], 0.5)
  expect_lt(v[2], 0.5)
})

test_that("feature_length sums a multi-part set", {
  lens <- c("19" = 84, "20" = 67, "21" = 221)
  expect_equal(feature_length("E019,E021", lens), 305)
})

test_that("feature_length returns NA for an empty distinct set", {
  lens <- c("19" = 84)
  expect_true(is.na(feature_length(NA_character_, lens)))
  expect_true(is.na(feature_length("NA", lens)))
  expect_true(is.na(feature_length("", lens)))
})

test_that("feature_length ignores parts absent from the annotation", {
  lens <- c("19" = 84)
  expect_equal(feature_length("E019,E999", lens), 84)
  expect_true(is.na(feature_length("E999", lens)))
})

test_that("a junction-substituted side gets NA, not a bogus per-base pi", {
  ## D is a junction: a point feature with no length
  expect_true(is.na(pi_perbase(0.72, feature_length(NA_character_, c("1" = 10)), 67)))
})
