# Tests for statistical and data-preparation functions.

# ---------------------------------------------------------------------------
# group_by_event
# ---------------------------------------------------------------------------

make_event_data <- function(n_per_event = 6) {
    # Two events, each with n_per_event rows (> 4 required to pass filter).
    n <- n_per_event * 2
    data.frame(
        gene   = rep("GENE1", n),
        event  = rep(c("e1", "e2"), each = n_per_event),
        groups = rep(c("grp1", "grp2"), n / 2),
        diff   = seq(5, by = 5, length.out = n),
        n      = seq(20, by = 5, length.out = n),
        stringsAsFactors = FALSE
    )
}

test_that("group_by_event splits data into one list element per event", {
    dat     <- make_event_data(6)
    grouped <- group_by_event(dat, "diff", "n")
    expect_equal(length(grouped), 2L)
})

test_that("group_by_event adds columns y and n to each group", {
    dat     <- make_event_data(6)
    grouped <- group_by_event(dat, "diff", "n")
    for (grp in grouped) {
        expect_true("y" %in% names(grp))
        expect_true("n" %in% names(grp))
    }
})

test_that("group_by_event filters out events with 4 or fewer samples", {
    # Exactly 4 rows per event — should be excluded (filter: n_samples > 4).
    dat <- data.frame(
        gene   = rep("GENE1", 8),
        event  = rep(c("e1", "e2"), each = 4),
        groups = rep(c("grp1", "grp2"), 4),
        diff   = 1:8,
        n      = 10:17,
        stringsAsFactors = FALSE
    )
    grouped <- group_by_event(dat, "diff", "n")
    expect_equal(length(grouped), 0L)
})

test_that("group_by_event keeps events with exactly 5 samples", {
    dat <- data.frame(
        gene   = rep("GENE1", 5),
        event  = rep("e1", 5),
        groups = c("grp1", "grp1", "grp1", "grp2", "grp2"),
        diff   = 1:5,
        n      = 10:14,
        stringsAsFactors = FALSE
    )
    grouped <- group_by_event(dat, "diff", "n")
    expect_equal(length(grouped), 1L)
})

test_that("group_by_event ignores rows where n is NA", {
    dat        <- make_event_data(6)
    dat$n[1:3] <- NA
    grouped    <- group_by_event(dat, "diff", "n")
    for (grp in grouped) {
        expect_false(any(is.na(grp$n)))
    }
})

# ---------------------------------------------------------------------------
# pvalueAdjustment (no independent filtering)
# ---------------------------------------------------------------------------

test_that("pvalueAdjustment adds a padj column", {
    res <- data.frame(
        gene   = c("g1", "g1", "g2", "g2"),
        event  = c("e1", "e2", "e3", "e4"),
        pvalue = c(0.01, 0.05, 0.10, 0.50)
    )
    out <- pvalueAdjustment(res,
                            independentFiltering = FALSE,
                            alpha                = 0.05,
                            pAdjustMethod        = "BH")
    expect_true("padj" %in% names(out))
    expect_equal(nrow(out), 4L)
})

test_that("pvalueAdjustment BH adjustment matches p.adjust", {
    pvals <- c(0.001, 0.01, 0.05, 0.20, 0.50)
    res   <- data.frame(pvalue = pvals)
    out   <- pvalueAdjustment(res,
                              independentFiltering = FALSE,
                              alpha                = 0.05,
                              pAdjustMethod        = "BH")
    expected_padj <- p.adjust(pvals, method = "BH")
    expect_equal(out$padj, expected_padj)
})

test_that("pvalueAdjustment preserves all input columns", {
    res <- data.frame(
        gene   = "g1",
        event  = "e1",
        pvalue = 0.03,
        effect = 1.5
    )
    out <- pvalueAdjustment(res,
                            independentFiltering = FALSE,
                            alpha                = 0.05,
                            pAdjustMethod        = "BH")
    expect_true("gene"   %in% names(out))
    expect_true("event"  %in% names(out))
    expect_true("effect" %in% names(out))
})

test_that("pvalueAdjustment with Bonferroni matches p.adjust", {
    pvals <- c(0.01, 0.05, 0.10)
    res   <- data.frame(pvalue = pvals)
    out   <- pvalueAdjustment(res,
                              independentFiltering = FALSE,
                              alpha                = 0.05,
                              pAdjustMethod        = "bonferroni")
    expect_equal(out$padj, p.adjust(pvals, method = "bonferroni"))
})

# ---------------------------------------------------------------------------
# moderate_phi_log_scale
# ---------------------------------------------------------------------------

make_phi_df <- function() {
    data.frame(
        gene    = paste0("g", 1:5),
        event   = paste0("e", 1:5),
        phi     = c(2.0, 3.0, 4.0, 5.0, 6.0),
        var_phi = c(0.1, 0.2, 0.3, 0.4, 0.5),
        stringsAsFactors = FALSE
    )
}

test_that("moderate_phi_log_scale adds required output columns", {
    result <- moderate_phi_log_scale(make_phi_df())
    expect_true(all(c("z", "var_z", "w", "z_mod", "phi_mod") %in% names(result)))
})

test_that("moderate_phi_log_scale computes z as log(phi)", {
    phi_df <- make_phi_df()
    result <- moderate_phi_log_scale(phi_df)
    expect_equal(result$z, log(phi_df$phi), tolerance = 1e-10)
})

test_that("moderate_phi_log_scale computes var_z via delta method", {
    phi_df <- make_phi_df()
    result <- moderate_phi_log_scale(phi_df)
    expected_var_z <- phi_df$var_phi / (phi_df$phi^2)
    expect_equal(result$var_z, expected_var_z, tolerance = 1e-10)
})

test_that("moderate_phi_log_scale returns strictly positive phi_mod", {
    result <- moderate_phi_log_scale(make_phi_df())
    expect_true(all(result$phi_mod > 0))
})

test_that("moderate_phi_log_scale shrinkage weights are in [0, 1]", {
    result <- moderate_phi_log_scale(make_phi_df())
    expect_true(all(result$w >= 0 & result$w <= 1))
})

test_that("moderate_phi_log_scale strips phi values outside (0, 1e10)", {
    phi_df <- rbind(
        make_phi_df(),
        data.frame(gene = "bad1", event = "bad1",
                   phi = -1,    var_phi = 0.1),   # negative phi → removed
        data.frame(gene = "bad2", event = "bad2",
                   phi = 2e10,  var_phi = 0.1)    # above limit → removed
    )
    result <- moderate_phi_log_scale(phi_df)
    # The two bad rows are removed before moderation.
    expect_equal(nrow(result), 5L)
})

test_that("moderate_phi_log_scale moderated phi lies between min and max observed phi", {
    phi_df <- make_phi_df()
    result <- moderate_phi_log_scale(phi_df)
    # Shrinkage always pulls estimates towards the mean, so phi_mod should
    # be within the observed range (or very close due to the tail shrinkage).
    expect_true(all(result$phi_mod >= min(phi_df$phi) * 0.5))
    expect_true(all(result$phi_mod <= max(phi_df$phi) * 2.0))
})

# ---------------------------------------------------------------------------
# test_model_wilcoxon
# ---------------------------------------------------------------------------

make_wilcoxon_data <- function(effect = TRUE) {
    # 5 samples per group, clear difference when effect = TRUE.
    if (effect) {
        y_g1 <- c(10L, 12L, 11L, 13L, 10L)   # low proportions ~0.10
        y_g2 <- c(45L, 43L, 44L, 46L, 45L)   # high proportions ~0.45
    } else {
        y_g1 <- c(10L, 11L, 12L, 10L, 11L)   # similar proportions
        y_g2 <- c(10L, 11L, 12L, 10L, 11L)
    }
    data.frame(
        gene   = rep("GENE1", 10),
        event  = rep("e1",    10),
        groups = c(rep("grp1", 5), rep("grp2", 5)),
        y      = c(y_g1, y_g2),
        n      = rep(100L, 10),
        stringsAsFactors = FALSE
    )
}

test_that("test_model_wilcoxon returns a data.frame with expected columns", {
    dd     <- make_wilcoxon_data(effect = TRUE)
    result <- test_model_wilcoxon(dd)
    expect_s3_class(result, "data.frame")
    expect_true("gene"        %in% names(result))
    expect_true("event"       %in% names(result))
    expect_true("p.value"     %in% names(result))
    expect_true("effect_size" %in% names(result))
    expect_equal(result$model, "wilcoxon")
})

test_that("test_model_wilcoxon gives a small p-value for a large effect", {
    dd     <- make_wilcoxon_data(effect = TRUE)
    result <- test_model_wilcoxon(dd)
    expect_lt(result$p.value, 0.05)
})

test_that("test_model_wilcoxon effect_size has correct sign and magnitude", {
    dd     <- make_wilcoxon_data(effect = TRUE)
    result <- test_model_wilcoxon(dd)
    # grp2 proportions (~0.45) > grp1 proportions (~0.10)
    # effect_size = mean(grp2) - mean(grp1) should be positive
    expect_gt(result$effect_size, 0)
    expect_gt(abs(result$effect_size), 0.2)
})

test_that("test_model_wilcoxon returns NULL when only one group is present", {
    dd        <- make_wilcoxon_data(effect = TRUE)
    dd$groups <- "grp1"   # collapse to one group
    result    <- test_model_wilcoxon(dd)
    expect_null(result)
})

test_that("test_model_wilcoxon returns NULL when all counts are zero", {
    dd    <- make_wilcoxon_data(effect = TRUE)
    dd$n  <- 0L
    result <- test_model_wilcoxon(dd)
    expect_null(result)
})

test_that("test_model_wilcoxon gene and event fields match input", {
    dd     <- make_wilcoxon_data(effect = TRUE)
    result <- test_model_wilcoxon(dd)
    expect_equal(result$gene,  "GENE1")
    expect_equal(result$event, "e1")
})

## --- per-base effect-size gating -------------------------------------------
## pi is a raw count ratio, so one absolute |delta_pi| threshold is not
## scale-fair across sides whose distinct set and reference differ in length.
## use_perbase gates on the normalized value where it exists, raw elsewhere.

mk_res <- function(dpi, pb, sd = "E001") {
  data.frame(padj = 1e-6, delta_pi = dpi, delta_pi_perbase = pb,
             lfc_diff_net = 5, comparison = "diff1_vs_ref",
             setdiff1 = sd, setdiff2 = NA_character_,
             d_support = 1e4, stringsAsFactors = FALSE)
}

test_that("use_perbase=FALSE keeps gating on the raw delta_pi", {
  r <- mk_res(dpi = 0.12, pb = 0.03)
  expect_true(is_significant(r, 0.01, 0, min_dpi = 0.1, use_perbase = FALSE))
})

test_that("use_perbase=TRUE gates on the normalized value instead", {
  ## raw clears 0.1 but per-base does not -> now fails
  r <- mk_res(dpi = 0.12, pb = 0.03)
  expect_false(is_significant(r, 0.01, 0, min_dpi = 0.1, use_perbase = TRUE))
})

test_that("use_perbase can also promote a side the raw value failed", {
  ## normalization makes effect sizes larger on average
  r <- mk_res(dpi = 0.08, pb = 0.15)
  expect_false(is_significant(r, 0.01, 0, min_dpi = 0.1, use_perbase = FALSE))
  expect_true(is_significant(r, 0.01, 0, min_dpi = 0.1, use_perbase = TRUE))
})

test_that("the default needs raw and per-base delta_pi >= 0.05 and |effect_size| >= 0.5", {
  r <- mk_res(dpi = 0.06, pb = 0.07); r$effect_size <- 0.6
  expect_true(is_significant(r, 0.01, 0))
  r$effect_size <- 0.4                       # log odds ratio too small
  expect_false(is_significant(r, 0.01, 0))
  r <- mk_res(dpi = 0.03, pb = 0.06); r$effect_size <- 2
  expect_false(is_significant(r, 0.01, 0))   # raw too small
  r <- mk_res(dpi = 0.08, pb = 0.04); r$effect_size <- 2
  expect_false(is_significant(r, 0.01, 0))   # per-base too small
  r <- mk_res(dpi = 0.06, pb = NA_real_, sd = NA_character_); r$effect_size <- 2
  expect_true(is_significant(r, 0.01, 0))    # junction side: raw against min_dpi_sj
  r$effect_size <- NA_real_                  # a missing effect_size passes, as delta_pi does
  expect_true(is_significant(r, 0.01, 0))
  r$effect_size <- 0.1
  expect_true(is_significant(r, 0.01, 0, min_effect = 0))   # 0 drops the condition
})

test_that("a junction side with no per-base value falls back to raw, not NA", {
  ## setdiff empty -> junction-substituted -> no length -> per-base is NA.
  ## It must be gated on raw against min_dpi_sj, never failed on NA.
  r <- mk_res(dpi = 0.07, pb = NA_real_, sd = NA_character_)
  expect_true(is_significant(r, 0.01, 0, min_dpi = 0.1, min_dpi_sj = 0.05,
                             use_perbase = TRUE))
})

test_that("use_perbase='both' fails a side either value fails", {
  expect_false(is_significant(mk_res(0.03, 0.15), 0.01, 0, min_dpi = 0.1, use_perbase = "both"))
  expect_false(is_significant(mk_res(0.15, 0.03), 0.01, 0, min_dpi = 0.1, use_perbase = "both"))
  expect_true(is_significant(mk_res(0.12, 0.15), 0.01, 0, min_dpi = 0.1, use_perbase = "both"))
  expect_error(is_significant(mk_res(0.12, 0.15), 0.01, 0, use_perbase = "or"))
})

test_that("use_perbase is inert when the column is absent", {
  r <- mk_res(dpi = 0.12, pb = 0.03)
  r$delta_pi_perbase <- NULL
  expect_true(is_significant(r, 0.01, 0, min_dpi = 0.1, use_perbase = TRUE))
})

test_that("add_significant threads use_perbase through", {
  r <- mk_res(dpi = 0.12, pb = 0.03)
  expect_true(add_significant(r, 0.01, 0, min_dpi = 0.1, use_perbase = FALSE)$significant)
  expect_false(add_significant(r, 0.01, 0, min_dpi = 0.1,
                               use_perbase = TRUE)$significant)
})

## --- a group with no reads on one side ---------------------------------------
## Its log odds are infinite, which used to fail the whole fit and drop every
## contrast. It is now left out: the other contrasts match a fit without it.

mk_dd <- function(groups, y, n) {
  data.frame(gene = "G", event = 1L, groups = factor(groups), y = y, n = n,
             z_mod = log(20), stringsAsFactors = FALSE)
}
mk_L <- function(trt, ref, lev) {
  C <- matrix(0, 1, length(lev), dimnames = list(NULL, paste0("groups", lev)))
  C[1, paste0("groups", trt)] <- 1; C[1, paste0("groups", ref)] <- -1
  C
}

test_that("an empty group drops only the contrasts that need it", {
  set.seed(1)
  g <- rep(c("A", "B", "C"), each = 6)
  n <- rep(100L, 18)
  y <- c(rbinom(6, 100, 0.3), rbinom(6, 100, 0.6), rep(0L, 6))
  L <- list(A_vs_B = mk_L("A", "B", c("A", "B", "C")), A_vs_C = mk_L("A", "C", c("A", "B", "C")))
  r <- test_model_glmmTMB_EB(mk_dd(g, y, n), L)
  expect_equal(r$contrast, "A_vs_B")
  k <- g != "C"
  r2 <- test_model_glmmTMB_EB(mk_dd(g[k], y[k], n[k]), L["A_vs_B"])
  expect_equal(r$p.value, r2$p.value, tolerance = 1e-6)
  expect_equal(r$effect_size, r2$effect_size, tolerance = 1e-6)
})

test_that("a group with every read on the distinct set is left out too", {
  g <- rep(c("A", "B", "C"), each = 4)
  y <- c(30L, 35L, 28L, 33L, 60L, 55L, 62L, 58L, 100L, 100L, 100L, 100L)
  L <- list(A_vs_B = mk_L("A", "B", c("A", "B", "C")), C_vs_B = mk_L("C", "B", c("A", "B", "C")))
  r <- test_model_glmmTMB_EB(mk_dd(g, y, rep(100L, 12)), L)
  expect_equal(r$contrast, "A_vs_B")
})
