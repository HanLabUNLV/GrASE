# --- Data Preparation Functions ---

#' Group exon count data by gene and event
#' @param dat A data frame of exon count data.
#' @param col_y Character string. Name of the column containing the count of reads supporting the
#'   alternative splicing event (numerator).
#' @param col_n Character string. Name of the column containing the total read count (denominator).
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("bipartition.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' grouped_data <- grase::group_by_event(splitcnts, col_y = "diff", col_n = "n")
#' length(grouped_data)
#' }
group_by_event <- function(dat, col_y, col_n) {
  dat$y <- dat[[col_y]]
  dat$n <- dat[[col_n]]
  dat <- dat[!is.na(dat$n),]
  if ("comparison" %in% names(dat)) {
    # Combined-comparison path: count samples per (gene, event, comparison)
    # and arrange so both comparisons for the same gene/event are adjacent.
    dat <- dat %>% dplyr::add_count(gene, event, comparison, name="n_samples")
    dat <- dat %>% dplyr::filter(n_samples > 4)
    grouped_data <- dat %>%
      dplyr::arrange(gene, event, comparison) %>%
      dplyr::group_by(gene, event, comparison) %>%
      dplyr::group_split()
  } else {
    dat <- dat %>% dplyr::add_count(gene, event, name="n_samples")
    dat <- dat %>% dplyr::filter(n_samples > 4)
    grouped_data <- dat %>%
      dplyr::group_by(gene, event) %>%
      dplyr::group_split()
  }
  return(grouped_data)
}

#' Estimate beta-binomial dispersion (phi) per event using glmmTMB
#' @param dd A data frame of exon count data.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("bipartition.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' grouped_data <- grase::group_by_event(splitcnts, "diff", "n")
#' phi_result <- grase::phi_estimate_glmmTMB(grouped_data[[1]])
#' phi_result
#' }
phi_estimate_glmmTMB <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2) return(NULL)

    # Estimate the single dispersion under the FULL design (mean modeled by
    # ~groups) rather than the intercept-only null (~1). Estimating dispersion
    # under ~1 lets a real group effect inflate the apparent overdispersion,
    # which biases the fixed-dispersion test toward being conservative under
    # the alternative (edgeR/DESeq2 estimate dispersion under the full design
    # for the same reason). The dispersion formula stays ~1, so this is still a
    # single phi -- just estimated controlling for the group means. Falls back
    # to ~1 when only one group is present.
    if ("groups" %in% names(dd)) dd$groups <- droplevels(factor(dd$groups))
    full_design <- "groups" %in% names(dd) && length(unique(dd$groups)) >= 2
    mean_form <- if (full_design) cbind(y, n - y) ~ groups else cbind(y, n - y) ~ 1

    # Binomial null (same mean design): used only to test whether the
    # beta-binomial dispersion is identifiable (significant overdispersion
    # beyond the group means; see 'identifiable' flag below).
    m0 <- tryCatch(
      glmmTMB(mean_form, data = dd, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    m1 <- tryCatch(
      glmmTMB(
        mean_form,
        data = dd,
        family = glmmTMB::betabinomial(link = "logit")
      ),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1))) {
      phi_hat <- as.numeric(sigma(m1))  # dispersion
      vc <- vcov(m1, full = TRUE)
      ## log(phi) and its variance from dispersion model
      if (nrow(vc) < 2 || ncol(vc) < 2) {
        var_phi <- NA_real_
      } else {
        var_log_phi <- vc["disp~(Intercept)", "disp~(Intercept)"]
        ## Delta method: Var(phi) = phi^2 * Var(log(phi))
        var_phi <- (phi_hat^2) * var_log_phi
      }

      # LRT: is the beta-binomial dispersion identifiable (i.e. is there
      # significant overdispersion relative to binomial)?  Events that fail
      # this are indistinguishable from binomial; their phi runs to the
      # +Inf (binomial-limit) boundary and carries no information about the
      # typical dispersion.  We flag them so they can be excluded when
      # estimating the global shrinkage target / prior, while still keeping
      # phi (they are tested with the moderated target, never dropped or
      # sent to a binomial GLM).
      identifiable <- FALSE
      if (!is.null(m0) && !is.na(logLik(m0))) {
        test <- tryCatch(anova(m1, m0, test = "LRT"), error = function(e) NULL)
        if (!is.null(test) && !is.na(test$`Pr(>Chisq)`[2]) &&
            test$`Pr(>Chisq)`[2] < 0.05) {
          identifiable <- TRUE
        }
      }

      return(data.frame(gene = gene, event = event, phi = phi_hat,
                        var_phi = var_phi, identifiable = identifiable))
    }
    return(NULL)
}


#' Estimate MAP log-phi per event using glmmTMB with an empirical Bayes prior
#' @param dd A data frame of exon count data.
#' @param prior_disp A data frame specifying the glmmTMB prior on the dispersion parameter.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("bipartition.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' grouped_data <- grase::group_by_event(splitcnts, "diff", "n")
#' prior_disp <- data.frame(prior = "normal(-1.5, 1.2)", class = "fixef_disp",
#'                          coef = "", stringsAsFactors = FALSE)
#' phi_map <- grase::phi_map_glmmTMB(grouped_data[[1]], prior_disp)
#' phi_map
#' }
phi_map_glmmTMB <- function(dd, prior_disp) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2) return(NULL)

    # MAP dispersion estimated under the FULL design (mean ~groups), matching
    # phi_estimate_glmmTMB, so the dispersion that EBmap fixes in the test is
    # not inflated by an unmodeled group effect. Dispersion formula stays ~1
    # (single phi); the prior on coef="1" still targets the dispersion intercept.
    if ("groups" %in% names(dd)) dd$groups <- droplevels(factor(dd$groups))
    full_design <- "groups" %in% names(dd) && length(unique(dd$groups)) >= 2
    mean_form <- if (full_design) cbind(y, n - y) ~ groups else cbind(y, n - y) ~ 1

    m <- tryCatch(
      glmmTMB(
        mean_form,
        data = dd,
        family = glmmTMB::betabinomial(link = "logit"),
        priors = prior_disp
      ),
      error = function(e) NULL
    )

    if (!is.null(m) && !is.na(logLik(m))) {
      return(data.frame(gene = gene, event = event,
                        z_mod = log(as.numeric(sigma(m)))))
    }
    return(NULL)
}


#' Moderate phi estimates on the log scale using empirical Bayes shrinkage
#' @param phi_table A data frame with columns \code{gene}, \code{event},
#'   \code{phi}, and \code{var_phi}, as returned by
#'   \code{phi_estimate_glmmTMB}.
#' @param trimming_limit Numeric. Upper bound for \code{phi} values; rows with
#'   \code{phi >= trimming_limit} are removed before shrinkage. Default is
#'   \code{1e+10}.
#' @export
#' @examples
#' \dontrun{
#' phi_df <- read.table("phi.glmmtmb.txt", header = TRUE, row.names = NULL)
#' phi_table <- grase::moderate_phi_log_scale(phi_df)
#' head(phi_table[, c("gene", "event", "phi", "phi_mod", "w")])
#' }
moderate_phi_log_scale <- function(phi_table, trimming_limit = 1e+10) {
  # 1. Basic Cleaning
  # Remove non-positive phis before log transformation
  phi_table <- phi_table[phi_table$phi > 0 & phi_table$phi < trimming_limit, ]
  
  # 2. Transform to Log-Space
  # Let z = log(phi)
  phi_table$z <- log(phi_table$phi)
  
  # 3. Transform Variance using the Delta Method
  # Var(log(phi)) is approximately Var(phi) / (phi^2)
  phi_table$var_z <- phi_table$var_phi / (phi_table$phi^2)
  
  # 4. Global Estimates (Using Medians in Log-Space)
  # Only events with an identifiable beta-binomial dispersion (significantly
  # overdispersed vs binomial) are allowed to define the global target and the
  # biological variance.  Boundary/non-identifiable events (phi at the binomial
  # limit) would otherwise inflate both z_bar and tau2_z (bimodal spread) and
  # prevent shrinkage.  All events still receive a z_mod below.
  if ("identifiable" %in% names(phi_table)) {
    use <- which(phi_table$identifiable %in% TRUE)
  } else {
    use <- seq_len(nrow(phi_table))
  }
  # Safety fallback: if too few identifiable events, use all of them.
  if (length(use) < 10) use <- seq_len(nrow(phi_table))

  z_bar      <- median(phi_table$z[use], na.rm = TRUE)
  s2_z       <- var(phi_table$z[use], na.rm = TRUE)
  # Use the median of sampling variances to represent the 'typical' noise
  # this will likely be ~1.0 instead of 51,000,000
  typical_var_z <- median(phi_table$var_z[use], na.rm = TRUE)

  # 5. Biological Variance (tau^2) in Log-Space
  # We subtract the average sampling error from the total observed variance
  tau2_z <- max(s2_z - typical_var_z, 0)
  
  if (tau2_z == 0) {
    warning("Biological variance in log-space is zero. All genes will shrink to the mean.")
  }
  
  # 6. Shrinkage Weights
  # w_g = biological_var / (biological_var + sampling_var_g)
  phi_table$w <- tau2_z / (tau2_z + phi_table$var_z)
  phi_table$w[!is.finite(phi_table$w)] <- 0
  
  # 7. Moderated Log-Phi
  phi_table$z_mod <- (phi_table$w * phi_table$z) + ((1 - phi_table$w) * z_bar)
  
  # 8. Back-transform to Original Scale
  phi_table$phi_mod <- exp(phi_table$z_mod)
  
  return(phi_table)
}


#' Trend-based phi moderation (analogous to DESeq2's dispersion trend).
#'
#' Fits a loess curve of log(phi) ~ log(baseMean) on events that have phi
#' estimates, then uses the trend as the shrinkage target instead of the
#' global mean.  Events with no phi estimate (non-convergent BB fit) are
#' fully shrunk to the trend prediction (weight = 0).
#'
#' @param phi_df      data.frame with columns gene, event, phi, var_phi
#' @param baseMean_df data.frame with columns gene, event, baseMean
#'                    (ALL events, not just those with phi estimates)
#' @param span        loess span parameter (default 0.5)
#' @param trimming_limit upper bound on phi before log-transform (default 1e10)
#' @return data.frame with all events in baseMean_df and columns
#'         z, var_z, z_trend, w, z_mod, phi_mod
#' Moderate phi estimates using a trend-based loess prior over baseMean
#' @export
#' @examples
#' \dontrun{
#' phi_df <- read.table("phi.glmmtmb.txt", header = TRUE, row.names = NULL)
#' splitcnts <- read.table("bipartition.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' baseMean_df <- dplyr::summarise(
#'   dplyr::group_by(splitcnts, gene, event),
#'   baseMean = mean(n, na.rm = TRUE), .groups = "drop"
#' )
#' phi_table <- grase::moderate_phi_trend(phi_df, baseMean_df)
#' head(phi_table[, c("gene", "event", "baseMean", "phi_mod")])
#' }
moderate_phi_trend <- function(phi_df, baseMean_df, span = 0.5,
                               trimming_limit = 1e+10) {
  # -- Step 1: prepare phi estimates --
  # Include "comparison" in join key when both inputs carry it, to avoid
  # a cross-product that would duplicate rows and inflate sample sizes.
  join_cols <- c("gene", "event")
  if ("comparison" %in% names(phi_df) && "comparison" %in% names(baseMean_df))
    join_cols <- c(join_cols, "comparison")

  phi_est <- phi_df[phi_df$phi > 0 & phi_df$phi < trimming_limit, ]
  phi_est <- phi_est %>%
    left_join(baseMean_df, by = join_cols) %>%
    filter(!is.na(baseMean) & baseMean > 0)
  phi_est$z      <- log(phi_est$phi)
  phi_est$var_z  <- phi_est$var_phi / (phi_est$phi^2)
  phi_est$log_bm <- log(phi_est$baseMean)

  # -- Step 2: fit loess trend on reliable estimates --
  # Exclude the noisiest 10 % of estimates when fitting the trend, and restrict
  # to events with an identifiable beta-binomial dispersion when that flag is
  # available (boundary/non-identifiable events at the binomial limit carry no
  # dispersion information and would distort the trend and tau^2).
  if ("identifiable" %in% names(phi_est)) {
    ident <- phi_est$identifiable %in% TRUE
    # Safety fallback: if too few identifiable events, keep all of them.
    if (sum(ident, na.rm = TRUE) < 10) ident <- rep(TRUE, nrow(phi_est))
  } else {
    ident <- rep(TRUE, nrow(phi_est))
  }
  var_z_thresh <- quantile(phi_est$var_z, 0.9, na.rm = TRUE)
  reliable     <- ident & is.finite(phi_est$z) & !is.na(phi_est$var_z) &
                  phi_est$var_z < var_z_thresh
  z_bar <- median(phi_est$z[reliable], na.rm = TRUE)

  trend_fit <- NULL
  if (sum(reliable, na.rm = TRUE) >= 10) {
    trend_fit <- tryCatch(
      loess(z ~ log_bm, data = phi_est[reliable, ], span = span),
      error = function(e) NULL
    )
  }

  # -- Step 3: predict trend for ALL events --
  all_events <- baseMean_df %>%
    mutate(log_bm = log(pmax(baseMean, 1)))

  if (!is.null(trend_fit)) {
    all_events$z_trend <- predict(trend_fit, newdata = all_events)
    # For out-of-range extrapolation, fall back to global median
    all_events$z_trend[is.na(all_events$z_trend)] <- z_bar
  } else {
    all_events$z_trend <- z_bar
  }

  # -- Step 4: estimate biological variance tau^2 --
  # Estimate the typical sampling noise and total variance from identifiable
  # events only, so the bimodal spread introduced by boundary estimates does
  # not inflate tau^2 (which would suppress shrinkage).
  typical_var_z <- median(phi_est$var_z[ident], na.rm = TRUE)
  s2_z          <- var(phi_est$z[ident], na.rm = TRUE)
  tau2_z        <- max(s2_z - typical_var_z, 0)

  # -- Step 5: join phi estimates onto all events, compute shrinkage --
  phi_est_cols <- c(join_cols, "z", "var_z")
  result <- all_events %>%
    left_join(phi_est %>% select(all_of(phi_est_cols)), by = join_cols)

  # Events without a phi estimate get w = 0 fully shrunk to the trend
  result$w <- ifelse(
    is.na(result$var_z), 0,
    tau2_z / (tau2_z + result$var_z)
  )
  result$w[!is.finite(result$w)] <- 0

  result$z_mod   <- result$w * result$z + (1 - result$w) * result$z_trend
  result$phi_mod <- exp(result$z_mod)

  return(result)
}


# --- DM Estimation & Moderation (NEW) ---


#' Estimate Dirichlet-Multinomial precision per event via direct 1-D grid log-likelihood optimization
#'
#' Estimates DM precision via direct 1-D log-likelihood optimization.  Under the intercept-only model the MLE of the proportion vector
#' is \eqn{\hat\pi_j = \sum_i y_{ij} / \sum_{ij} y_{ij}}, so precision
#' estimation reduces to a one-dimensional optimization over \eqn{\log\alpha}.
#' Variance is obtained from the analytical observed Fisher information at the
#' MLE.  uses only base R (\code{lgamma},
#' \code{trigamma}, \code{optimize}).
#'
#' @param dd A data frame for a single gene-event group with columns
#'   \code{gene}, \code{event}, \code{sample}, \code{groups}, \code{type},
#'   and \code{count}.
#' @return A one-row data frame with columns \code{gene}, \code{event},
#'   \code{log_prec} (log-scale MLE of the DM precision \eqn{\alpha}), and
#'   \code{var_log_prec} (estimated variance of \code{log_prec}), or
#'   \code{NULL} on failure.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("multinomial.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' grouped_counts <- dplyr::group_split(dplyr::group_by(splitcnts, gene, event))
#' prec_result <- grase::prec_estimate_plugin_dm(grouped_counts[[1]])
#' prec_result
#' }
prec_estimate_plugin_dm <- function(dd) {
    gene <- unique(dd$gene); event <- unique(dd$event)

    # --- 1. Build count matrix Y ---
    wide_df <- dd %>%
        dplyr::select(sample, groups, type, count) %>%
        pivot_wider(names_from = type, values_from = count, values_fill = 0)
    Y <- as.matrix(wide_df[, setdiff(names(wide_df), c("sample", "groups"))])
    wide_df <- wide_df[rowSums(Y) > 0, ]
    Y <- Y[rowSums(Y) > 0, , drop = FALSE]
    if (nrow(Y) < 2 || ncol(Y) < 2) return(NULL)
    if (sum(colSums(Y) > 0) < 2) return(NULL)

    # --- 2. Closed-form MLE of proportions under null model ---
    pi_hat <- colSums(Y) / sum(Y)   # length-K vector
    n_i    <- rowSums(Y)             # per-sample totals
    N      <- nrow(Y)

    # --- 3. DM log-likelihood as a function of log(alpha) ---
    # ll = N*lgamma(alpha) - sum_i lgamma(n_i + alpha)
    #    + sum_j [ sum_i lgamma(y_ij + alpha*pi_j) - N*lgamma(alpha*pi_j) ]
    dm_ll <- function(log_alpha) {
        alpha    <- exp(log_alpha)
        alpha_pi <- alpha * pi_hat                           # K-vector
        row_terms <- N * lgamma(alpha) - sum(lgamma(n_i + alpha))
        mat_terms <- sum(lgamma(sweep(Y, 2, alpha_pi, "+")) -
                         matrix(lgamma(alpha_pi), nrow = N, ncol = ncol(Y), byrow = TRUE))
        row_terms + mat_terms
    }

    # --- 4. 1-D maximisation; interval covers alpha in [4.5e-5, 1.6e5] ---
    opt <- tryCatch(
        optimize(dm_ll, interval = c(-10, 12), maximum = TRUE),
        error = function(e) NULL
    )
    if (is.null(opt)) return(NULL)

    log_alpha_hat <- opt$maximum
    alpha_hat     <- exp(log_alpha_hat)
    alpha_pi      <- alpha_hat * pi_hat

    # --- 5. Analytical second derivative at the MLE ---
    # d2ll/dalpha2 = N*trigamma(alpha) - sum_i trigamma(n_i + alpha)
    #              + sum_j pi_j^2 * [ sum_i trigamma(y_ij + alpha*pi_j)
    #                                 - N * trigamma(alpha*pi_j) ]
    d2_dalpha2 <- N * trigamma(alpha_hat) - sum(trigamma(n_i + alpha_hat)) +
                  sum(pi_hat^2 * (colSums(trigamma(sweep(Y, 2, alpha_pi, "+"))) -
                                  N * trigamma(alpha_pi)))

    # Log-scale: d2ll/d(log alpha)^2 ~= alpha^2 * d2ll/dalpha2  (gradient ~ 0 at MLE)
    d2_dlogalpha2 <- alpha_hat^2 * d2_dalpha2

    # var_log_prec = 1 / (- d2ll/d(log alpha)^2)
    if (!is.finite(d2_dlogalpha2) || d2_dlogalpha2 >= 0) return(NULL)
    var_log_prec <- 1 / (-d2_dlogalpha2)
    if (!is.finite(var_log_prec) || var_log_prec <= 0) return(NULL)

    # --- 6. Identifiability: DM vs multinomial (alpha -> Inf) LRT ---
    # The multinomial limit is the null with no over-dispersion (rho -> 0).
    # Its log-likelihood is the closed-form multinomial LL at pi_hat. Events
    # whose DM fit is not significantly better than multinomial (or whose alpha
    # ran to the optimizer's upper boundary) are non-identifiable: their
    # precision carries no over-dispersion information and would bias the global
    # shrinkage target if included. We flag them so the moderation step can
    # exclude them from prior estimation while still returning a precision.
    multinom_ll <- sum(Y * log(matrix(pi_hat, nrow = N, ncol = ncol(Y),
                                      byrow = TRUE)), na.rm = TRUE)
    LR      <- 2 * (opt$objective - multinom_ll)
    # Boundary of the optimize() interval (see step 4); a fit that reaches it is
    # unidentified regardless of the LRT.
    at_boundary <- log_alpha_hat >= 12 - 1e-3
    identifiable <- is.finite(LR) && LR > 0 &&
                    pchisq(LR, df = 1, lower.tail = FALSE) < 0.05 &&
                    !at_boundary

    data.frame(gene = gene, event = event,
               log_prec = log_alpha_hat, var_log_prec = var_log_prec,
               identifiable = identifiable)
}


#' Moderate Dirichlet-Multinomial precision estimates using empirical Bayes
#' shrinkage
#' @param prec_table A data frame with columns \code{gene}, \code{event},
#'   \code{log_prec}, and \code{var_log_prec}, as returned by
#'   \code{prec_estimate_plugin_dm}.
#' @export
#' @examples
#' \dontrun{
#' prec_table <- read.table("prec_dm.txt", header = TRUE, row.names = NULL)
#' prec_table <- grase::moderate_prec_log_scale(prec_table)
#' head(prec_table[, c("gene", "event", "log_prec", "log_prec_mod", "rho_mod")])
#' }
moderate_prec_log_scale <- function(prec_table) {
  # 1. Filter for robust estimation of the prior (ignore failed fits)
  # log_prec is logit(rho). Values like -2e6 are numerical garbage.
  # var_log_prec > 1000 indicates essentially no information.
  valid_subset <- prec_table$var_log_prec < 1000 & abs(prec_table$log_prec) < 50
  # Additionally restrict to events with an identifiable DM precision (when the
  # flag is available): events at the multinomial (no-over-dispersion) boundary
  # pass the var/magnitude filter but form a spurious second mode that biases
  # z_bar and tau2. They still receive a moderated precision below.
  if ("identifiable" %in% names(prec_table)) {
    ident_valid <- valid_subset & (prec_table$identifiable %in% TRUE)
    if (sum(ident_valid, na.rm = TRUE) >= 10) valid_subset <- ident_valid
  }

  # If too few valid points, fallback to simple median or original
  if (sum(valid_subset, na.rm=TRUE) < 10) {
      z_bar <- median(prec_table$log_prec, na.rm=TRUE)
      tau2_z <- 0 
  } else {
      # 2. Estimate Prior Parameters from valid subset
      # Use Median for center to be robust against outliers
      z_bar <- median(prec_table$log_prec[valid_subset], na.rm = TRUE)
      
      # Estimate biological variance (tau2)
      # Total Variance of estimates
      s2_z <- var(prec_table$log_prec[valid_subset], na.rm = TRUE)
      # Typical sampling variance (use median to be robust against outliers)
      typical_var_z <- median(prec_table$var_log_prec[valid_subset], na.rm = TRUE)
      
      tau2_z <- max(s2_z - typical_var_z, 0)
  }

  # 3. Calculate Weights for ALL data points
  # w = tau2 / (tau2 + sampling_var)
  # If sampling_var is huge (garbage fit), w -> 0, and we shrink to z_bar
  prec_table$w_prec <- tau2_z / (tau2_z + prec_table$var_log_prec)
  prec_table$w_prec[is.na(prec_table$w_prec)] <- 0
  
  # 4. Moderate
  prec_table$log_prec_mod <- (prec_table$w_prec * prec_table$log_prec) + ((1 - prec_table$w_prec) * z_bar)
  prec_table$prec_mod <- exp(prec_table$log_prec_mod)
  
  # 5. Convert to rho
  # rho = 1 / (1 + A)
  rho_mod <- 1 / (1 + prec_table$prec_mod)
  prec_table$rho_mod <- pmax(pmin(rho_mod, 1 - 1e-6), 1e-6)

  return(prec_table)
}


#' Trend-based Dirichlet-Multinomial precision moderation
#'
#' Fits a loess curve of \code{log(prec) ~ log(baseMean)} on events that have
#' precision estimates, then uses the trend as the shrinkage target instead of
#' the global mean.  Events without a precision estimate (non-convergent fit or
#' not in the estimation subsample) are fully shrunk to the trend prediction
#' (weight = 0).
#'
#' @param prec_df     data.frame with columns \code{gene}, \code{event},
#'   \code{log_prec}, \code{var_log_prec}, as returned by
#'   \code{prec_estimate_plugin_dm}.
#' @param baseMean_df data.frame with columns \code{gene}, \code{event},
#'   \code{baseMean} covering ALL events (not just those with prec estimates).
#' @param span loess span parameter (default 0.5).
#' @return data.frame with one row per event in \code{baseMean_df} and columns
#'   \code{gene}, \code{event}, \code{baseMean}, \code{z_trend}, \code{w_prec},
#'   \code{log_prec_mod}, \code{prec_mod}, \code{rho_mod}.
#' @export
#' @examples
#' \dontrun{
#' prec_df <- read.table("prec_dm.txt", header = TRUE, row.names = NULL)
#' splitcnts <- read.table("multinomial.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' baseMean_df <- splitcnts %>%
#'   dplyr::group_by(gene, event, sample) %>%
#'   dplyr::summarise(total = sum(count), .groups = "drop") %>%
#'   dplyr::group_by(gene, event) %>%
#'   dplyr::summarise(baseMean = mean(total), .groups = "drop")
#' prec_table <- grase::moderate_prec_trend(prec_df, baseMean_df)
#' head(prec_table[, c("gene", "event", "baseMean", "log_prec_mod", "rho_mod")])
#' }
moderate_prec_trend <- function(prec_df, baseMean_df, span = 0.5) {
  # -- Step 1: prepare valid precision estimates --
  # Include "comparison" in join key when both inputs carry it, to avoid
  # a cross-product that would duplicate rows and inflate sample sizes.
  join_cols <- c("gene", "event")
  if ("comparison" %in% names(prec_df) && "comparison" %in% names(baseMean_df))
    join_cols <- c(join_cols, "comparison")

  valid     <- abs(prec_df$log_prec) < 50 & prec_df$var_log_prec < 1000
  # Restrict the trend/prior-fitting set to identifiable events when flagged:
  # boundary (multinomial-limit) estimates would otherwise distort the trend
  # and inflate tau2. All events still receive a moderated precision below.
  if ("identifiable" %in% names(prec_df) &&
      sum(valid & (prec_df$identifiable %in% TRUE), na.rm = TRUE) >= 10) {
    valid <- valid & (prec_df$identifiable %in% TRUE)
  }
  prec_est  <- prec_df[valid, ] %>%
    left_join(baseMean_df, by = join_cols) %>%
    filter(!is.na(baseMean) & baseMean > 0)
  prec_est$log_bm <- log(prec_est$baseMean)

  # -- Step 2: fit loess trend on the most reliable estimates --
  # Exclude noisy top 10% by sampling variance
  var_thresh <- quantile(prec_est$var_log_prec, 0.9, na.rm = TRUE)
  reliable   <- is.finite(prec_est$log_prec) &
                !is.na(prec_est$var_log_prec) &
                prec_est$var_log_prec < var_thresh
  z_bar <- median(prec_est$log_prec[reliable], na.rm = TRUE)

  trend_fit <- NULL
  if (sum(reliable, na.rm = TRUE) >= 10) {
    trend_fit <- loess(log_prec ~ log_bm, data = prec_est[reliable, ], span = span)
  }

  # -- Step 3: predict trend for ALL events --
  all_events <- baseMean_df %>%
    mutate(log_bm = log(pmax(baseMean, 1)))

  if (!is.null(trend_fit)) {
    all_events$z_trend <- predict(trend_fit, newdata = all_events)
    # Fall back to global median for out-of-range extrapolation
    all_events$z_trend[is.na(all_events$z_trend)] <- z_bar
  } else {
    all_events$z_trend <- z_bar
  }

  # -- Step 4: estimate biological variance tau^2 --
  typical_var <- median(prec_est$var_log_prec, na.rm = TRUE)
  s2          <- var(prec_est$log_prec, na.rm = TRUE)
  tau2        <- max(s2 - typical_var, 0)

  # -- Step 5: join estimates onto all events and compute shrinkage weights --
  prec_est_cols <- c(join_cols, "log_prec", "var_log_prec")
  result <- all_events %>%
    left_join(prec_est %>% select(all_of(prec_est_cols)), by = join_cols)

  # Events without an estimate (w=0) are fully shrunk to the trend
  result$w_prec <- ifelse(
    is.na(result$var_log_prec), 0,
    tau2 / (tau2 + result$var_log_prec)
  )
  result$w_prec[!is.finite(result$w_prec)] <- 0

  result$log_prec_mod <- result$w_prec * result$log_prec +
                         (1 - result$w_prec) * result$z_trend
  result$prec_mod     <- exp(result$log_prec_mod)
  rho_mod             <- 1 / (1 + result$prec_mod)
  result$rho_mod      <- pmax(pmin(rho_mod, 1 - 1e-6), 1e-6)

  return(result)
}



# --- Moderated Testing Functions ---

# 1. glmmTMB Beta-Binomial without prior (for comparison)
#' Test differential exon usage with a beta-binomial glmmTMB model, no prior on dispersion.
#' Intended for comparison with \code{test_model_glmmTMB_EB}.
#' @param dd A data frame of exon count data.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("bipartition.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' splitcnts$groups <- factor(splitcnts$groups)
#' grouped_data <- grase::group_by_event(splitcnts, "diff", "n")
#' result <- grase::test_model_glmmTMB_without_prior(grouped_data[[1]])
#' result
#' }
test_model_glmmTMB_without_prior <- function(dd, L, model_label = "betabinom_MLE") {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
    if (mean(dd$y) == 0) return(NULL)

    dd$groups <- droplevels(dd$groups)
    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 0 + groups, data = dd,
              family = glmmTMB::betabinomial(link = "logit")),
      error = function(e) NULL
    )
    if (is.null(m1) || is.na(logLik(m1))) return(NULL)

    beta     <- fixef(m1)$cond
    V        <- vcov(m1)$cond
    coef_nms <- names(beta)
    phi_val  <- sigma(m1)

    ## L is a named list of contrast MATRICES: 1 row for a pairwise trt-vs-ref
    ## contrast, K-1 rows for a K-group omnibus. A 1-row matrix reproduces the
    ## former (est/se)^2 with df = 1 exactly, so pairwise output is unchanged.
    ## effect_size is the contrast estimate for a pair and NA for an omnibus,
    ## which has no single direction -- the pairwise effect sizes of its
    ## constituent pairs carry that information instead.
    results <- lapply(names(L), function(ctr_name) {
      C       <- L[[ctr_name]]
      if (is.null(dim(C))) C <- matrix(C, nrow = 1L, dimnames = list(NULL, names(C)))
      needed  <- colnames(C)[apply(C != 0, 2, any)]
      if (!all(needed %in% coef_nms)) return(NULL)
      C_use   <- C[, coef_nms, drop = FALSE]
      w       <- wald_contrast(beta, V, C_use)
      if (is.null(w) || !is.finite(w$stat)) return(NULL)
      est     <- if (nrow(C_use) == 1L) as.numeric(C_use %*% beta[coef_nms]) else NA_real_
      data.frame(gene = gene, event = event, contrast = ctr_name,
                 LRT = w$stat, p.value = w$p, df = w$df,
                 model = model_label, phi = phi_val, effect_size = est,
                 stringsAsFactors = FALSE)
    })
    bind_rows(Filter(Negate(is.null), results))
}


#' Wald test for a contrast vector or a contrast matrix
#'
#' Generalises the 1-df contrast test to a K-1 df omnibus test over K groups.
#' A contrast VECTOR reproduces the 1-df test exactly -- W = (est/se)^2 with
#' df = 1 -- so pairwise results are unchanged to the last digit. A contrast
#' MATRIX with r independent rows gives the omnibus statistic
#'   W = (C b)' (C V C')^{-1} (C b),  df = rank(C).
#'
#' Uses an eigen-based generalised inverse because C V C' is rank-deficient
#' whenever a group is absent or aliased in this gene/event, which is common:
#' df is then the numerical rank, not nrow(C).
#'
#' @param beta named coefficient vector from the fitted model.
#' @param V coefficient covariance matrix, dimnames matching \code{beta}.
#' @param C contrast vector (named) or matrix (named columns).
#' @return list(stat, df, p), or NULL when no direction is estimable.
#' @export
wald_contrast <- function(beta, V, C) {
  if (is.null(dim(C))) C <- matrix(C, nrow = 1L, dimnames = list(NULL, names(C)))
  keep <- colnames(C)
  if (!all(keep %in% names(beta))) return(NULL)
  est <- C %*% beta[keep]
  M   <- C %*% V[keep, keep, drop = FALSE] %*% t(C)
  ei  <- eigen(M, symmetric = TRUE)
  tol <- max(ei$values) * 1e-10
  pos <- ei$values > tol
  if (!any(pos)) return(NULL)
  U    <- ei$vectors[, pos, drop = FALSE]
  Minv <- U %*% diag(1 / ei$values[pos], sum(pos)) %*% t(U)
  W    <- as.numeric(t(est) %*% Minv %*% est)
  df   <- sum(pos)
  list(stat = W, df = df, p = stats::pchisq(W, df = df, lower.tail = FALSE))
}

# 2. Moderated glmmTMB Beta-Binomial EB (EBapprox or EBmap)
#' Test differential exon usage with a beta-binomial glmmTMB model using fixed moderated dispersion.
#' @param dd A data frame of exon count data with column \code{z_mod} (log-scale moderated phi).
#' @param model_label Character. Label written to the \code{model} column of the result.
#'   Use \code{"betabinom_EBapprox"} for closed-form Gaussian-approximation moderation and
#'   \code{"betabinom_EBmap"} for MAP-optimized moderation. Default is \code{"betabinom_EBapprox"}.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("bipartition.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' phi_df <- read.table("phi.glmmtmb.txt", header = TRUE, row.names = NULL)
#' phi_table <- grase::moderate_phi_log_scale(phi_df)
#' splitcnts_eb <- dplyr::left_join(splitcnts, phi_table, by = c("gene", "event"))
#' splitcnts_eb$groups <- factor(splitcnts_eb$groups)
#' grouped_eb <- grase::group_by_event(splitcnts_eb, "diff", "n")
#' result <- grase::test_model_glmmTMB_EB(grouped_eb[[1]])
#' result
#' }
test_model_glmmTMB_EB <- function(dd, L, model_label = "betabinom_EBapprox") {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
    if (mean(dd$y) == 0) return(NULL)

    target_log_phi <- dd$z_mod[1]
    if (is.na(target_log_phi)) return(NULL)

    dd$groups <- droplevels(dd$groups)
    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 0 + groups,
        data = dd,
        family = glmmTMB::betabinomial(link = "logit"),
        start = list(betadisp = target_log_phi),
        map = list(betadisp = factor(NA))),
      error = function(e) NULL
    )
    if (is.null(m1) || is.na(logLik(m1))) return(NULL)

    beta     <- fixef(m1)$cond
    V        <- vcov(m1)$cond
    coef_nms <- names(beta)
    phi_val  <- sigma(m1)

    ## L is a named list of contrast MATRICES: 1 row for a pairwise trt-vs-ref
    ## contrast, K-1 rows for a K-group omnibus. A 1-row matrix reproduces the
    ## former (est/se)^2 with df = 1 exactly, so pairwise output is unchanged.
    ## effect_size is the contrast estimate for a pair and NA for an omnibus,
    ## which has no single direction -- the pairwise effect sizes of its
    ## constituent pairs carry that information instead.
    results <- lapply(names(L), function(ctr_name) {
      C       <- L[[ctr_name]]
      if (is.null(dim(C))) C <- matrix(C, nrow = 1L, dimnames = list(NULL, names(C)))
      needed  <- colnames(C)[apply(C != 0, 2, any)]
      if (!all(needed %in% coef_nms)) return(NULL)
      C_use   <- C[, coef_nms, drop = FALSE]
      w       <- wald_contrast(beta, V, C_use)
      if (is.null(w) || !is.finite(w$stat)) return(NULL)
      est     <- if (nrow(C_use) == 1L) as.numeric(C_use %*% beta[coef_nms]) else NA_real_
      data.frame(gene = gene, event = event, contrast = ctr_name,
                 LRT = w$stat, p.value = w$p, df = w$df,
                 model = model_label, phi = phi_val, effect_size = est,
                 stringsAsFactors = FALSE)
    })
    bind_rows(Filter(Negate(is.null), results))
}





# 3b. Wilcoxon Rank-Sum Test (Nonparametric)
#' Test for differential exon usage using Wilcoxon rank-sum test
#' @param dd A data frame of exon count data.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("bipartition.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' splitcnts$groups <- factor(splitcnts$groups)
#' grouped_data <- grase::group_by_event(splitcnts, "diff", "n")
#' result <- grase::test_model_wilcoxon(grouped_data[[1]])
#' result
#' }
test_model_wilcoxon <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    
    # Filter valid rows
    dd <- dd[dd$n > 0, ]

    # Ensure there are enough samples and both groups are present
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
    if (mean(dd$y) == 0) return(NULL)
    
    # Calculate proportions
    dd$prop <- dd$y / dd$n

    # Check variation in groups
    counts <- table(dd$groups)
    if (length(counts) < 2 || any(counts < 1)) return(NULL)

    res <- tryCatch(
        wilcox.test(prop ~ groups, data = dd),
        error = function(e) NULL
    )
    
    if (!is.null(res)) {
        # Calculate Delta pi (path proportion; difference in means)
        # Using factor levels to determine order: Level 2 - Level 1
        grps <- levels(factor(dd$groups))
        m1_val <- mean(dd$prop[dd$groups == grps[1]])
        m2_val <- mean(dd$prop[dd$groups == grps[2]])
        eff_size <- m2_val - m1_val

        return(data.frame(
            gene = gene, 
            event = event, 
            LRT = NA, 
            p.value = res$p.value, 
            model = "wilcoxon", 
            phi = NA,
            effect_size = eff_size
        ))
    }
    return(NULL)
}


# 4b. Direct DM LRT with fixed EB precision
#' Test differential exon usage with a direct Dirichlet-multinomial LRT using empirical Bayes fixed precision.
#'
#' Direct Dirichlet-multinomial LRT using empirical Bayes fixed precision.
#' When precision is fixed (EB-moderated), the MLE of the proportion vector under each model
#' has a closed form (empirical marginal counts), so both null and alternative log-likelihoods
#' can be evaluated directly without iterative fitting.
#'
#' @param dd A data frame of exon count data with columns \code{gene}, \code{event},
#'   \code{sample}, \code{groups}, \code{type}, \code{count}, \code{log_prec_mod}.
#' @return A one-row data frame with columns \code{gene}, \code{event}, \code{LRT},
#'   \code{p.value}, \code{model}, \code{effect_size}, or \code{NULL} on failure.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("multinomial.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' prec_table <- grase::moderate_prec_log_scale(
#'   read.table("prec_dm.txt", header = TRUE, row.names = NULL)
#' )
#' splitcnts_eb <- dplyr::left_join(splitcnts, prec_table, by = c("gene", "event"))
#' grouped_eb <- dplyr::group_split(dplyr::group_by(splitcnts_eb, gene, event))
#' result <- grase::test_model_multinomial_plugin_dm_EB(grouped_eb[[1]])
#' result
#' }
test_model_multinomial_plugin_dm_EB <- function(dd) {
    gene <- unique(dd$gene); event <- unique(dd$event)
    log_prec_mod <- dd$log_prec_mod[1]
    if (is.na(log_prec_mod)) return(NULL)
    alpha <- exp(log_prec_mod)

    wide_df <- dd %>% dplyr::select(sample, groups, type, count) %>%
               pivot_wider(names_from = type, values_from = count, values_fill = 0)
    Y <- as.matrix(wide_df[, setdiff(names(wide_df), c("sample", "groups"))])
    keep <- rowSums(Y) > 0
    wide_df <- wide_df[keep, ]
    Y <- Y[keep, , drop = FALSE]

    if (nrow(Y) < 2 || length(unique(wide_df$groups)) < 2) return(NULL)
    if (sum(colSums(Y) > 0) < 2) return(NULL)

    K    <- ncol(Y)
    n_i  <- rowSums(Y)
    grps <- wide_df$groups
    g_levels <- unique(grps)

    # DM log-likelihood for a count matrix Y with row totals n_i,
    # fixed precision alpha, and proportion vector pi
    dm_ll_fixed <- function(Y, n_i, alpha, pi) {
        pi     <- pmax(pi, 1e-300)
        pi     <- pi / sum(pi)
        ap     <- alpha * pi
        sum(lgamma(alpha) - lgamma(n_i + alpha)) +
          sum(lgamma(sweep(Y, 2, ap, "+")) -
              matrix(lgamma(ap), nrow = nrow(Y), ncol = K, byrow = TRUE))
    }

    # Null: common proportions
    pi_null <- colSums(Y) / sum(Y)
    ll_null <- dm_ll_fixed(Y, n_i, alpha, pi_null)
    if (!is.finite(ll_null)) return(NULL)

    # Alt: per-group proportions (closed-form MLE)
    ll_alt <- {
        s <- 0
        for (g in g_levels) {
            idx  <- grps == g
            Y_g  <- Y[idx, , drop = FALSE]
            ni_g <- n_i[idx]
            pi_g <- colSums(Y_g) / sum(Y_g)
            s    <- s + dm_ll_fixed(Y_g, ni_g, alpha, pi_g)
        }
        s
      } 
    if (!is.finite(ll_alt)) return(NULL)

    LR   <- max(2 * (ll_alt - ll_null), 0)
    df   <- (K - 1) * (length(g_levels) - 1)   # = K-1 for 2 groups
    pval <- pchisq(LR, df = df, lower.tail = FALSE)

    # Effect size: log ratio of group proportions for category 1 vs reference (last)
    pi1  <- colSums(Y[grps == g_levels[1], , drop = FALSE]) / sum(Y[grps == g_levels[1], , drop = FALSE])
    pi2  <- colSums(Y[grps == g_levels[2], , drop = FALSE]) / sum(Y[grps == g_levels[2], , drop = FALSE])
    pi1  <- pmax(pi1, 1e-10); pi2 <- pmax(pi2, 1e-10)
    eff_size <- log(pi2[1] / pi2[K]) - log(pi1[1] / pi1[K])

    data.frame(gene = gene, event = event, LRT = LR, p.value = pval,
               model = "dirmult_EBplugin", effect_size = eff_size)
}


# --- Denominator-Effect Filter ---

#' Compute per-event log2 fold changes for diff and ref counts
#'
#' For each (gene, event), computes mean diff and mean ref counts per condition,
#' then the log2 fold change (treatment / reference) for each. Used to distinguish
#' true DTU (diff drives the ratio shift) from denominator-effect false positives
#' (ref drives the ratio shift because a DTE transcript passes through the shared
#' denominator).
#'
#' Library-size bias cancels when comparing abs(lfc_diff) vs abs(lfc_ref): both
#' are computed from the same samples, so any per-sample depth offset is identical
#' in both LFCs and disappears in the comparison.
#'
#' @param splitcnts A data frame with columns gene, event, sample, groups, diff, ref.
#'   Must be bipartition or n_choose_2 format (not multinomial). The groups column
#'   must be a factor with levels c(cond2, cond1) so that LFC = log2(cond1 / cond2).
#' @param pseudocount Integer pseudocount added to means before log2 to avoid log(0).
#'   Default 1L.
#' @return A data frame with columns gene, event, lfc_diff, lfc_ref, lfc_diff_net.
#'   lfc_diff_net = abs(lfc_diff) - abs(lfc_ref): positive values indicate the diff
#'   exon changes more than the ref exon (genuine DTU signal); negative values suggest
#'   the ref exon is the likely driver (denominator effect FP).
#'   Filter criterion: keep events where lfc_diff_net > delta.
#' @export
compute_lfc_summary <- function(splitcnts, pseudocount = 1L,
                                cond_ref = levels(splitcnts$groups)[1],
                                cond_trt = levels(splitcnts$groups)[2]) {

  splitcnts %>%
    dplyr::group_by(gene, event, groups) %>%
    dplyr::summarise(
      mean_diff = mean(diff, na.rm = TRUE),
      mean_ref  = mean(ref,  na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from  = groups,
      values_from = c(mean_diff, mean_ref)
    ) %>%
    dplyr::mutate(
      lfc_diff     = log2((.data[[paste0("mean_diff_", cond_trt)]] + pseudocount) /
                          (.data[[paste0("mean_diff_", cond_ref)]] + pseudocount)),
      lfc_ref      = log2((.data[[paste0("mean_ref_",  cond_trt)]] + pseudocount) /
                          (.data[[paste0("mean_ref_",  cond_ref)]] + pseudocount)),
      lfc_diff_net = abs(lfc_diff) - abs(lfc_ref),
      pi_trt       = .data[[paste0("mean_diff_", cond_trt)]] /
                       (.data[[paste0("mean_diff_", cond_trt)]] + .data[[paste0("mean_ref_", cond_trt)]]),
      pi_ref       = .data[[paste0("mean_diff_", cond_ref)]] /
                       (.data[[paste0("mean_diff_", cond_ref)]] + .data[[paste0("mean_ref_", cond_ref)]]),
      delta_pi     = pi_trt - pi_ref
    ) %>%
    dplyr::select(gene, event, lfc_diff, lfc_ref, lfc_diff_net, pi_trt, pi_ref, delta_pi)
}

#' Annotate test results with denominator-effect flag columns
#'
#' Joins the output of \code{compute_lfc_summary} onto a test results data frame,
#' adding columns lfc_diff, lfc_ref.
#'
#' @param results A test results data frame with columns gene and event.
#' @param lfc_summary A data frame as returned by \code{compute_lfc_summary}.
#' @return results with lfc_diff, lfc_ref columns added.
#' @export
posthoc_lfc_summary <- function(results, lfc_summary) {
  dplyr::left_join(results, lfc_summary, by = c("gene", "event"))
}

#' Add a significant column to test results
#'
#' @param res data frame with padj, lfc_diff_net, delta_pi columns
#' @param padj_thr FDR threshold
#' @param delta minimum lfc_diff_net (directionality filter; events with lfc_diff_net <= delta excluded)
#' @param min_dpi minimum |delta_pi| (path-proportion effect size filter)
#' @export
#' @param min_dpi_sj Numeric. |delta_pi| threshold for sides whose distinct set
#'   came from SPLIT READS rather than exonic parts (merged exon+SJ runs). Split
#'   -read distinct counts sit on a much lower pi scale than the exonic shared
#'   reference -- measured mean pi_ref 0.120 vs 0.332 on the DICE internal arm --
#'   so one absolute threshold is not scale-fair: |delta_pi| >= 0.1 demands a
#'   2.07x odds shift from a junction-sourced side against 1.53x from an exonic
#'   one, and passes them at 3.7% vs 14.7%. 0.05 at pi 0.120 is an odds ratio of
#'   1.47, which matches the exonic side's 1.53 at 0.1. Source is detected from
#'   setdiff1/setdiff2 being NA; when those columns are absent (pre-annotation
#'   call sites) every row uses min_dpi.
#' The GrASE call rule, as a predicate
#'
#' The single authoritative definition of "significant". A call requires an
#' adjusted p-value below \code{padj_thr}, a net log fold change above
#' \code{delta}, and an absolute path-proportion shift of at least the dpi
#' threshold. Returns a logical vector so callers sweeping a threshold grid
#' (PR curves, ROC, metric tables) use the same rule as the tool itself
#' instead of re-implementing it.
#'
#' @param res data frame with \code{padj}, and optionally \code{lfc_diff_net},
#'   \code{delta_pi}, \code{comparison}, \code{setdiff1}, \code{setdiff2}.
#'   Absent columns drop their condition rather than failing.
#' @param padj_thr adjusted p-value threshold.
#' @param delta minimum \code{lfc_diff_net}. Note this is SIGNED and that is
#'   deliberate: \code{lfc_diff_net = abs(lfc_diff) - abs(lfc_ref)}, so the
#'   magnitudes are already inside the quantity and it is positive when the
#'   distinct set moved more than the shared reference. Do not wrap in abs().
#' @param use_perbase Gate on \code{delta_pi_perbase} wherever it is defined,
#'   falling back to the raw \code{delta_pi} elsewhere. pi is a raw count
#'   ratio, so one absolute threshold is not scale-fair across sides whose
#'   distinct set and reference differ in length. The fallback is required, not
#'   optional: a junction-substituted side has a point feature as its distinct
#'   set, so it has no length and no per-base pi.
#' @param min_dpi minimum \code{abs(delta_pi)} for exon-sourced sides.
#' @param min_dpi_sj minimum \code{abs(delta_pi)} for sides whose distinct set
#'   came from split reads. DEFAULTS TO \code{min_dpi} (a uniform gate). Set it
#'   lower for a scale-fair gate: split-read distinct counts sit on a lower pi
#'   scale than the exonic reference, so one absolute threshold demands a larger
#'   odds shift from a junction-sourced side than an exonic one.
#' @param min_reads minimum read support for the tested distinct set, required
#'   in AT LEAST ONE of the two groups being contrasted. Read from a
#'   \code{d_support} column (the max over the contrast's two group means);
#'   no-op when that column is absent.
#'
#'   EITHER, not both: requiring both groups would exclude on/off switches -- a
#'   path absent in one condition by design -- which is the strongest splicing
#'   signal there is and which the ground truth scores as maximal movement. On
#'   the simulation a both-groups bar at 10 reads discarded 169 more true
#'   positives than an either-group bar, for one extra false positive.
#'
#'   Applied to every side, not only split-read-sourced ones: restricting it to
#'   junction sides saved only 16 true positives out of 692 at identical
#'   precision, which does not justify a source-dependent rule.
#' @param padj optional vector to use in place of \code{res$padj} -- for a
#'   two-stage nested-BH sweep, pass \code{pmax(within_gene_BH, padj_gene)}.
#' @return logical vector, one element per row of \code{res}.
#'
#' @details NA handling: a missing \code{lfc_diff_net} or \code{delta_pi}
#'   PASSES its condition, so the call rests on padj alone. A missing padj never
#'   passes. This is the tool's long-standing convention and is preserved here;
#'   on the current simulation no row has either NA, so it is latent.
#' @export
is_significant <- function(res, padj_thr, delta = 0, min_dpi = 0.1,
                           min_dpi_sj = min_dpi, min_reads = 0, padj = NULL,
                           omnibus_pairs = NULL, use_perbase = FALSE) {
  ## use_perbase is last in the signature on purpose: existing positional calls
  ## must keep working.
  pv <- if (is.null(padj)) res$padj else padj
  n  <- nrow(res)

  lfc_ok <- if ("lfc_diff_net" %in% names(res)) {
    (res$lfc_diff_net > delta) | is.na(res$lfc_diff_net)
  } else TRUE

  dpi_ok <- if ("delta_pi" %in% names(res)) {
    thr <- rep(min_dpi, n)
    if (!identical(min_dpi_sj, min_dpi) &&
        all(c("comparison", "setdiff1", "setdiff2") %in% names(res))) {
      sd <- ifelse(grepl("diff1", res$comparison), res$setdiff1, res$setdiff2)
      # nzchar(NA) is TRUE, so test is.na() explicitly rather than relying on it
      thr[is.na(sd) | sd %in% c("NA", "")] <- min_dpi_sj
    }
    ## With use_perbase, gate on the LENGTH-NORMALIZED delta_pi wherever it is
    ## defined and fall back to the raw one elsewhere. pi is a raw count ratio,
    ## so one absolute threshold is not scale-fair across sides whose distinct
    ## set and reference differ in length; per-base makes the threshold mean the
    ## same thing everywhere it can be computed.
    ##
    ## The fallback is not optional: a junction-substituted side has a point
    ## feature as its distinct set, so it has no length and no per-base pi
    ## (30.7% of rows on the DICE activation arm). Requiring per-base
    ## everywhere would fail those sides on NA rather than on effect size.
    ## Note min_dpi_sj already carries a partial scale correction for them.
    dpi <- abs(res$delta_pi)
    if (use_perbase && "delta_pi_perbase" %in% names(res)) {
      pb <- abs(res$delta_pi_perbase)
      dpi <- ifelse(is.na(pb), dpi, pb)
    }
    (dpi >= thr) | is.na(res$delta_pi)
  } else TRUE

  ## read support for the tested distinct set, in EITHER contrasted group.
  ## Applied to every side regardless of source: restricting it to junction
  ## sides saved only 16 of 692 true positives at identical precision.
  sup_ok <- if (min_reads > 0 && "d_support" %in% names(res)) {
    !is.na(res$d_support) & res$d_support >= min_reads
  } else TRUE

  ## K-group omnibus rows have no single direction, so delta_pi and
  ## lfc_diff_net are NA and the two effect-size gates cannot be evaluated on
  ## the row itself. Gate them instead on whether ANY constituent pairwise
  ## contrast clears both, for the same gene/event/side. Those gates stay
  ## pairwise deliberately: generalising delta_pi to a range over K groups
  ## makes it an extreme-order statistic, and a fixed floor then loosens as K
  ## grows (measured on DICE: the share clearing 0.1 rises from 7.5% at K=2 to
  ## 32.5% at K=13). omnibus_pairs maps each omnibus contrast name to the
  ## pairwise contrast names it spans; those rows must be present in `res`.
  ## Identify omnibus rows from the test itself: df > 1 means a K-1 df Wald
  ## contrast, which has no single direction. Do NOT rely on omnibus_pairs to
  ## identify them -- delta_pi and lfc_diff_net are NA on such rows, and both
  ## gates treat NA as a pass (so that frames lacking those columns are not
  ## filtered), so an unidentified omnibus row would clear both vacuously.
  is_omni <- if ("df" %in% names(res)) !is.na(res$df) & res$df > 1 else rep(FALSE, n)
  if (any(is_omni) && is.null(omnibus_pairs)) {
    ## Cannot verify the effect size for these rows: fail them rather than let
    ## the NA-passes-through rule call them significant on padj alone.
    if (length(lfc_ok) == 1L) lfc_ok <- rep(lfc_ok, n)
    if (length(dpi_ok) == 1L) dpi_ok <- rep(dpi_ok, n)
    lfc_ok[is_omni] <- FALSE; dpi_ok[is_omni] <- FALSE
  }
  if (!is.null(omnibus_pairs) && all(c("contrast", "gene", "event") %in% names(res))) {
    is_omni <- is_omni | as.character(res$contrast) %in% names(omnibus_pairs)
    if (any(is_omni)) {
      pass_pair <- lfc_ok & dpi_ok
      if (length(pass_pair) == 1L) pass_pair <- rep(pass_pair, n)
      cmp <- if ("comparison" %in% names(res)) as.character(res$comparison) else rep("", n)
      key <- paste(res$gene, res$event, cmp, res$contrast, sep = "\r")
      lut <- stats::setNames(pass_pair, key)
      for (i in which(is_omni)) {
        prs <- omnibus_pairs[[as.character(res$contrast[i])]]
        hit <- lut[paste(res$gene[i], res$event[i], cmp[i], prs, sep = "\r")]
        ok  <- any(hit %in% TRUE)
        if (length(lfc_ok) > 1L) lfc_ok[i] <- ok else lfc_ok <- replace(rep(lfc_ok, n), i, ok)
        if (length(dpi_ok) > 1L) dpi_ok[i] <- ok else dpi_ok <- replace(rep(dpi_ok, n), i, ok)
      }
    }
  }

  !is.na(pv) & pv < padj_thr & lfc_ok & dpi_ok & sup_ok
}

#' Add the significant column
#'
#' Thin wrapper on \code{\link{is_significant}}, which holds the rule.
#' @inheritParams is_significant
#' @return \code{res} with a logical \code{significant} column.
#' @export
add_significant <- function(res, padj_thr, delta, min_dpi = 0.1, min_dpi_sj = 0.05,
                            min_reads = 10, omnibus_pairs = NULL,
                            use_perbase = FALSE) {
  res$significant <- is_significant(res, padj_thr, delta, min_dpi, min_dpi_sj,
                                    min_reads, omnibus_pairs = omnibus_pairs,
                                    use_perbase = use_perbase)
  res
}

#' Impute missing z_mod values with per-comparison or global median fallback
#'
#' @param df data frame with comparison and z_mod columns
#' @param p_table phi table with z_mod values used for fallback
#' @export
impute_z_mod <- function(df, p_table) {
  if (nrow(df) == 0) return(df)
  comp_name  <- df$comparison[1]
  fallback_z <- median(p_table$z_mod[p_table$comparison == comp_name], na.rm = TRUE)
  if (is.na(fallback_z)) fallback_z <- median(p_table$z_mod, na.rm = TRUE)
  if (is.na(fallback_z)) fallback_z <- 0
  df$z_mod[is.na(df$z_mod)] <- fallback_z
  df
}

#' Run a test function over events in chunked parallel batches with checkpointing
#'
#' @param sc split counts data frame
#' @param test_fn function applied per event
#' @param err_log path to error log file
#' @param chunk_size number of events per parallel batch
#' @param checkpoint_prefix prefix for checkpoint RDS files; NULL disables checkpointing
#' @param mc_cores number of parallel cores
#' @export
run_one_comparison <- function(sc, test_fn, err_log, ..., chunk_size = 10000L,
                               checkpoint_prefix = NULL, mc_cores = 1L) {
  if (nrow(sc) == 0) return(NULL)
  gd <- group_by_event(sc, 'diff', 'n')
  n_events <- length(gd)
  chunk_starts <- seq(1L, n_events, by = chunk_size)
  n_chunks <- length(chunk_starts)
  all_results <- vector("list", n_chunks)
  for (ci in seq_along(chunk_starts)) {
    ckpt_file <- if (!is.null(checkpoint_prefix))
      sprintf("%s_chunk%dof%d.rds", checkpoint_prefix, ci, n_chunks) else NULL
    if (!is.null(ckpt_file) && file.exists(ckpt_file)) {
      all_results[[ci]] <- readRDS(ckpt_file)
      message(sprintf("[%s] LRT chunk %d/%d loaded from checkpoint",
                      format(Sys.time(), "%H:%M:%S"), ci, n_chunks))
      next
    }
    idx <- chunk_starts[ci]:min(chunk_starts[ci] + chunk_size - 1L, n_events)
    chunk_res <- parallel::mclapply(gd[idx], function(dd) {
      result <- tryCatch(test_fn(dd, ...), error = function(e) {
        msg <- sprintf("[%s] ERROR gene=%s event=%s (PID %d): %s\n",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       dd$gene[1], dd$event[1], Sys.getpid(), conditionMessage(e))
        cat(msg, file = err_log, append = TRUE)
        NULL
      })
      if (!is.null(result) && "comparison" %in% names(dd))
        result$comparison <- dd$comparison[1]
      result
    }, mc.cores = mc_cores)
    all_results[[ci]] <- dplyr::bind_rows(Filter(Negate(is.null), chunk_res))
    if (!is.null(ckpt_file)) saveRDS(all_results[[ci]], ckpt_file)
    message(sprintf("[%s] LRT testing chunk %d/%d done (%d events)",
                    format(Sys.time(), "%H:%M:%S"), ci, n_chunks, length(idx)))
  }
  res <- dplyr::bind_rows(all_results)
  if (nrow(res) > 0 && !"comparison" %in% names(res))
    res$comparison <- sc$comparison[1]
  res
}

#' Run phi_map_glmmTMB in parallel for one comparison with checkpointing
#'
#' @param sc split counts data frame for one comparison
#' @param comp_name comparison label (e.g. "diff1_vs_ref")
#' @param prior_fn function taking a grouped event data frame, returning phi_map_glmmTMB result
#' @param err_log path to error log file
#' @param checkpoint_prefix prefix for checkpoint RDS files; NULL disables checkpointing
#' @param mc_cores number of parallel cores
#' @export
run_map_for_comparison <- function(sc, comp_name, prior_fn, err_log,
                                   checkpoint_prefix = NULL, mc_cores = 1L) {
  gd <- group_by_event(sc, 'diff', 'n')
  n_events     <- length(gd)
  chunk_size   <- 10000L
  chunk_starts <- seq(1L, n_events, by = chunk_size)
  n_chunks     <- length(chunk_starts)
  chunks <- vector("list", n_chunks)
  for (ci in seq_along(chunk_starts)) {
    ckpt_file <- if (!is.null(checkpoint_prefix))
      sprintf("%s_chunk%dof%d.rds", checkpoint_prefix, ci, n_chunks) else NULL
    if (!is.null(ckpt_file) && file.exists(ckpt_file)) {
      chunks[[ci]] <- readRDS(ckpt_file)
      message(sprintf("[%s] MAP phi chunk %d/%d loaded from checkpoint (%s)",
                      format(Sys.time(), "%H:%M:%S"), ci, n_chunks, comp_name))
      next
    }
    idx <- chunk_starts[ci]:min(chunk_starts[ci] + chunk_size - 1L, n_events)
    chunks[[ci]] <- parallel::mclapply(gd[idx], function(dd) {
      result <- tryCatch(prior_fn(dd), error = function(e) {
        cat(sprintf("[%s] WARN gene=%s event=%s: %s\n",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    dd$gene[1], dd$event[1], conditionMessage(e)),
            file = err_log, append = TRUE)
        NULL
      })
      if (!is.null(result)) result$comparison <- comp_name
      result
    }, mc.cores = mc_cores)
    if (!is.null(ckpt_file)) saveRDS(chunks[[ci]], ckpt_file)
    message(sprintf("[%s] MAP phi moderation chunk %d/%d done (%d events, %s)",
                    format(Sys.time(), "%H:%M:%S"), ci, n_chunks, length(idx), comp_name))
  }
  dplyr::bind_rows(Filter(Negate(is.null), unlist(chunks, recursive = FALSE)))
}

#' Estimate phi (overdispersion) for all events in one comparison using glmmTMB
#'
#' @param sc split counts data frame for one comparison
#' @param comp_name comparison label used for log file names
#' @param outdir output directory for progress and error logs
#' @param mc_cores number of parallel cores
#' @export
estimate_phi_for_comparison <- function(sc, comp_name, outdir, mc_cores = 1L) {
  if (nrow(sc) == 0) return(NULL)
  grouped_data      <- group_by_event(sc, 'diff', 'n')
  phi_progress_file <- file.path(outdir, paste0("phi_progress_", comp_name, ".log"))
  error_log_file    <- paste0("phi.glmmtmb.errors.", comp_name, ".log")

  n_events     <- length(grouped_data)
  chunk_size   <- 10000L
  chunk_starts <- seq(1L, n_events, by = chunk_size)
  phi_chunks   <- vector("list", length(chunk_starts))
  for (ci in seq_along(chunk_starts)) {
    idx <- chunk_starts[ci]:min(chunk_starts[ci] + chunk_size - 1L, n_events)
    phi_chunks[[ci]] <- parallel::mclapply(grouped_data[idx], function(dd) {
      result <- withCallingHandlers(
        tryCatch({
          phi_estimate_glmmTMB(dd)
        }, error = function(e) {
          msg <- sprintf("[%s] ERROR gene=%s event=%s (PID %d): %s\n",
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                         dd$gene[1], dd$event[1], Sys.getpid(), conditionMessage(e))
          cat(msg, file = error_log_file, append = TRUE)
          NULL
        }),
        warning = function(w) {
          msg <- sprintf("[%s] WARN  gene=%s event=%s (PID %d): %s\n",
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                         dd$gene[1], dd$event[1], Sys.getpid(), conditionMessage(w))
          cat(msg, file = error_log_file, append = TRUE)
          invokeRestart("muffleWarning")
        }
      )
      cat(sprintf("%s\t%s\n", dd$gene[1], dd$event[1]),
          file = phi_progress_file, append = TRUE)
      result
    }, mc.cores = mc_cores)
    message(sprintf("[%s] phi chunk %d/%d done (%d events)",
                    format(Sys.time(), "%H:%M:%S"), ci,
                    length(chunk_starts), length(idx)))
  }
  phi_list     <- unlist(phi_chunks, recursive = FALSE)
  phi_df_comp  <- dplyr::bind_rows(Filter(Negate(is.null), phi_list))
  if (nrow(phi_df_comp) > 0) phi_df_comp$comparison <- comp_name
  phi_df_comp
}
