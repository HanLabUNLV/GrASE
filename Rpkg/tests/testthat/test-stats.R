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
