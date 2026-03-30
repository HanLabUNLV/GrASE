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
  dat <- dat %>% dplyr::add_count(gene, event, name="n_samples")
  dat <- dat %>% dplyr::filter(n_samples > 4)
  grouped_data <- dat %>%
    dplyr::group_by(gene, event) %>%
    dplyr::group_split()
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

    m1 <- tryCatch(
      glmmTMB(
        cbind(y, n - y) ~ 1,
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
      return(data.frame(gene = gene, event = event, phi = phi_hat, var_phi = var_phi))
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
  z_bar      <- median(phi_table$z, na.rm = TRUE)
  s2_z       <- var(phi_table$z, na.rm = TRUE)
  # Use the median of sampling variances to represent the 'typical' noise
  # this will likely be ~1.0 instead of 51,000,000
  typical_var_z <- median(phi_table$var_z, na.rm = TRUE)
  
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
  phi_est <- phi_df[phi_df$phi > 0 & phi_df$phi < trimming_limit, ]
  phi_est <- phi_est %>%
    left_join(baseMean_df, by = c("gene", "event")) %>%
    filter(!is.na(baseMean) & baseMean > 0)
  phi_est$z      <- log(phi_est$phi)
  phi_est$var_z  <- phi_est$var_phi / (phi_est$phi^2)
  phi_est$log_bm <- log(phi_est$baseMean)

  # -- Step 2: fit loess trend on reliable estimates --
  # Exclude the noisiest 10 % of estimates when fitting the trend
  var_z_thresh <- quantile(phi_est$var_z, 0.9, na.rm = TRUE)
  reliable     <- is.finite(phi_est$z) & !is.na(phi_est$var_z) &
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
  typical_var_z <- median(phi_est$var_z, na.rm = TRUE)
  s2_z          <- var(phi_est$z, na.rm = TRUE)
  tau2_z        <- max(s2_z - typical_var_z, 0)

  # -- Step 5: join phi estimates onto all events, compute shrinkage --
  result <- all_events %>%
    left_join(phi_est %>% select(gene, event, z, var_z), by = c("gene", "event"))

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

    data.frame(gene = gene, event = event,
               log_prec = log_alpha_hat, var_log_prec = var_log_prec)
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
  valid     <- abs(prec_df$log_prec) < 50 & prec_df$var_log_prec < 1000
  prec_est  <- prec_df[valid, ] %>%
    left_join(baseMean_df, by = c("gene", "event")) %>%
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
  result <- all_events %>%
    left_join(prec_est %>% select(gene, event, log_prec, var_log_prec),
              by = c("gene", "event"))

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

# 1. Moderated glmmTMB Beta-Binomial with Prior
#' Test differential exon usage with a beta-binomial glmmTMB model using an empirical Bayes prior on dispersion.
#' @param dd A data frame of exon count data.
#' @param prior_disp A data frame specifying the glmmTMB prior on the dispersion parameter, with
#'   columns \code{prior}, \code{class}, and \code{coef} (as produced by \code{glmmTMB::prior}).
#' @param z_bar Numeric. The global empirical Bayes log-phi estimate used as a fixed-dispersion
#'   fallback when the prior-based model does not converge. Pass \code{NULL} to disable the
#'   fallback. Default is \code{NULL}.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("bipartition.internal.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' splitcnts$groups <- factor(splitcnts$groups)
#' grouped_data <- grase::group_by_event(splitcnts, "diff", "n")
#' prior_disp <- data.frame(prior = "normal(-1.5, 1.2)", class = "fixef_disp",
#'                          coef = "", stringsAsFactors = FALSE)
#' result <- grase::test_model_glmmTMB_with_prior(grouped_data[[1]], prior_disp, z_bar = -1.5)
#' result
#' }
test_model_glmmTMB_with_prior <- function(dd, prior_disp, z_bar = NULL) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
    if (mean(dd$y) == 0) return(NULL)

    # Models include 'priors' for empirical Bayes moderation
    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ groups, data = dd, family = glmmTMB::betabinomial(link="logit"), priors = prior_disp, REML=FALSE),
      error = function(e) NULL
    )
    m0 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd, family = glmmTMB::betabinomial(link="logit"), priors = prior_disp, REML=FALSE),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0))) {

        # If phi escaped to an extreme value despite the prior, the optimizer likely
        # diverged rather than reflecting truly binomial data. Fall through to the
        # fixed-dispersion fallback below instead of using a binomial GLM.
        if (sigma(m1) <= 1e5) {
            LR <- 2 * (logLik(m1) - logLik(m0))
            pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
            eff_size <- fixef(m1)$cond[2]
            return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pval, model="betabinomial_glmmTMB_MAP_with_prior", phi=sigma(m1), effect_size=eff_size))
        }
    }

    # Fallback: fix dispersion at the global EB median log-phi and refit.
    # This is conservative and consistent with the glmmTMB_fixedEB approach.
    # A binomial GLM fallback would be anticonservative for overdispersed data.
    if (!is.null(z_bar) && is.finite(z_bar)) {
        m1 <- tryCatch(
          glmmTMB(cbind(y, n - y) ~ groups, data = dd,
                  family = glmmTMB::betabinomial(link = "logit"),
                  start = list(betadisp = z_bar),
                  map   = list(betadisp = factor(NA))),
          error = function(e) NULL
        )
        m0 <- tryCatch(
          glmmTMB(cbind(y, n - y) ~ 1, data = dd,
                  family = glmmTMB::betabinomial(link = "logit"),
                  start = list(betadisp = z_bar),
                  map   = list(betadisp = factor(NA))),
          error = function(e) NULL
        )
        if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0))) {
            LR <- 2 * (logLik(m1) - logLik(m0))
            pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
            eff_size <- fixef(m1)$cond[2]
            return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pval, model="betabinomial_glmmTMB_fixed_median", phi=sigma(m1), effect_size=eff_size))
        }
    }
    return(NULL)
}



# 1b. glmmTMB Beta-Binomial without prior (for comparison)
#' Test differential exon usage with a beta-binomial glmmTMB model, no prior on dispersion.
#' Intended for comparison with \code{test_model_glmmTMB_with_prior}.
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
test_model_glmmTMB_without_prior <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
    if (mean(dd$y) == 0) return(NULL)

    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ groups, data = dd,
              family = glmmTMB::betabinomial(link = "logit")),
      error = function(e) NULL
    )
    m0 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd,
              family = glmmTMB::betabinomial(link = "logit")),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0))) {
        LR <- 2 * (logLik(m1) - logLik(m0))
        pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
        eff_size <- fixef(m1)$cond[2]
        return(data.frame(gene = gene, event = event, LRT = as.numeric(LR),
                          p.value = pval, model = "betabinomial_glmmTMB_no_prior",
                          phi = sigma(m1), effect_size = eff_size))
    }
    return(NULL)
}


# 2. Moderated glmmTMB Beta-Binomial EB
#' Test differential exon usage with a beta-binomial glmmTMB model using empirical Bayes fixed dispersion.
#' @param dd A data frame of exon count data.
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
test_model_glmmTMB_EB <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
    if (mean(dd$y) == 0) return(NULL)

    target_log_phi <- dd$z_mod[1]
    if (is.na(target_log_phi)) return(NULL)

    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ groups, 
        data = dd, 
        family = glmmTMB::betabinomial(link="logit"), 
        start = list(betadisp = target_log_phi),
        map = list(betadisp = factor(NA))), # This freezes the dispersion
      error = function(e) NULL
    )
    m0 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, 
        data = dd, 
        family = glmmTMB::betabinomial(link="logit"), 
        start = list(betadisp = target_log_phi),
        map = list(betadisp = factor(NA))),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0))) {
        LR <- 2 * (logLik(m1) - logLik(m0))
        pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
        eff_size <- fixef(m1)$cond[2]
        return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pval, model="betabinomial_glmmTMB_fixed_EB", phi=sigma(m1), effect_size=eff_size))
    }
    return(NULL)
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
        # Calculate Delta PSI (Difference in Means)
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
               model = "dirmult_moderated_direct", effect_size = eff_size)
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
#' @return A data frame with columns gene, event, lfc_diff, lfc_ref, denom_flag.
#'   denom_flag is TRUE when abs(lfc_ref) > abs(lfc_diff), indicating the ref is
#'   the likely driver (denominator-effect FP).
#' @export
compute_lfc_summary <- function(splitcnts, pseudocount = 1L) {
  lvls     <- levels(splitcnts$groups)
  cond_ref <- lvls[1]   # reference condition (denominator of LFC)
  cond_trt <- lvls[2]   # treatment condition  (numerator  of LFC)

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
      lfc_diff   = log2((.data[[paste0("mean_diff_", cond_trt)]] + pseudocount) /
                        (.data[[paste0("mean_diff_", cond_ref)]] + pseudocount)),
      lfc_ref    = log2((.data[[paste0("mean_ref_",  cond_trt)]] + pseudocount) /
                        (.data[[paste0("mean_ref_",  cond_ref)]] + pseudocount)),
      denom_flag = abs(lfc_ref) > abs(lfc_diff)
    ) %>%
    dplyr::select(gene, event, lfc_diff, lfc_ref, denom_flag)
}

#' Annotate test results with denominator-effect flag columns
#'
#' Joins the output of \code{compute_lfc_summary} onto a test results data frame,
#' adding columns lfc_diff, lfc_ref, and denom_flag.
#'
#' @param results A test results data frame with columns gene and event.
#' @param lfc_summary A data frame as returned by \code{compute_lfc_summary}.
#' @return results with lfc_diff, lfc_ref, denom_flag columns added.
#' @export
flag_denominator_effect <- function(results, lfc_summary) {
  dplyr::left_join(results, lfc_summary, by = c("gene", "event"))
}
