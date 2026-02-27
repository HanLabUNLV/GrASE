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

  # Events without a phi estimate get w = 0 → fully shrunk to the trend
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

#' Estimate Dirichlet-Multinomial precision per event using VGAM
#' @param dd A data frame of exon count data for a single gene-event group,
#'   with columns \code{gene}, \code{event}, \code{sample}, \code{groups},
#'   \code{type}, and \code{count}.
#' @export
#' @examples
#' \dontrun{
#' splitcnts <- read.table("multinomial.exoncnt.combined.txt",
#'                         header = TRUE, row.names = NULL)
#' grouped_counts <- dplyr::group_split(dplyr::group_by(splitcnts, gene, event))
#' prec_result <- grase::prec_estimate_vgam(grouped_counts[[1]])
#' prec_result
#' }
prec_estimate_vgam <- function(dd) {
    gene <- unique(dd$gene); event <- unique(dd$event)
    wide_df <- dd %>% dplyr::select(sample, groups, type, count) %>% pivot_wider(names_from = type, values_from = count, values_fill = 0)
    Y <- as.matrix(wide_df[, setdiff(names(wide_df), c("sample", "groups"))])
    wide_df <- wide_df[rowSums(Y) > 0, ]
    Y <- Y[rowSums(Y) > 0, , drop = FALSE]
    if (nrow(Y) < 2 || ncol(Y) < 2) return(NULL)
    if (sum(colSums(Y) > 0) < 2) return(NULL)

    m_init <- vglm(Y ~ 1, dirmultinomial, data = wide_df)

    log_prec_hat <- as.numeric(coef(m_init)[ncol(Y)])
    vc <- vcov(m_init)
    if (nrow(vc) < ncol(Y)) return(NULL)

    return(data.frame(gene = gene, event = event, log_prec = log_prec_hat, var_log_prec = vc[ncol(Y), ncol(Y)]))
}


#' Moderate Dirichlet-Multinomial precision estimates using empirical Bayes
#' shrinkage
#' @param prec_table A data frame with columns \code{gene}, \code{event},
#'   \code{log_prec}, and \code{var_log_prec}, as returned by
#'   \code{prec_estimate_vgam}.
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




# 3. EB-moderated VGAM beta-binomial
#' Test differential exon usage with a VGAM beta-binomial model initialized from empirical Bayes dispersion.
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
#' result <- grase::test_model_vgam_EB_init(grouped_eb[[1]])
#' result
#' }
test_model_vgam_EB_init <- function(dd) {

  ## per-gene moderated LRT in VGAM
  gene  <- unique(dd$gene)
  event <- unique(dd$event)
  dd <- dd[dd$n > 0, ]
  if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
  if (mean(dd$y) == 0) return(NULL)

  phi_mod <- dd$phi_mod[1]
  if (is.na(phi_mod)) return(NULL)

  phi_hat = dd$phi[1]
  phi_mod = dd$phi_mod[1]
  rho_mod = 1 / (1 + phi_mod)
  # Nudge slightly away from boundaries to prevent logit errors
  rho_mod = min(max(rho_mod, 1e-6), 1 - 1e-6)
  w = dd$w[1]

  m1 <- tryCatch(
    vglm(cbind(y, n - y) ~ groups,
         family = VGAM::betabinomial(irho = rho_mod, zero = 2),
         data = dd),
    error = function(e) NULL
  )

  m0 <- tryCatch(
    vglm(cbind(y, n - y) ~ 1,
         family = VGAM::betabinomial(irho = rho_mod, zero = 2),
         data = dd),
    error = function(e) NULL
  )

  if (is.null(m1) || is.null(m0)) return(NULL)

  intercept_2 <- coef(m1)["(Intercept):2"]
  rho <- VGAM::logitlink(intercept_2, inverse = TRUE)
  phi <- (1 / rho) - 1
  # this can result in extremely large phi estimates. 
  # but found it is not an inssue as it robustly estimates the phi and results in more conservative LRT compared to falling back to binomial glm. 
  # no need to limit phi here unlike in glmmTMB where large phi can cause convergence issues and unreliable LRT.
  cfs <- coef(m1)
  eff_size <- cfs[grep("groups", names(cfs))][1] 

  LR <- 2 * (logLik(m1) - logLik(m0))

  if (as.numeric(LR) == -Inf) {
    eff_size <- NA      # remove abnormally high effect sizes when splice event is perfectly correlated with the group (e.g., Group A has 0 reads for the event in all samples, while Group B has counts). complete separation
  }

  data.frame(
    gene    = gene,
    event   = event,
    LRT     = as.numeric(LR),
    p.value = pchisq(as.numeric(LR), df = 1, lower.tail = FALSE),
    phi_hat = phi_hat,
    phi_mod = phi_mod,
    w       = w,
    rho = rho,
    phi = phi,
    model   = "betabinomial_MLE_EB_anchor",
    effect_size = eff_size,
    stringsAsFactors = FALSE
  )

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


# 4. Moderated VGAM Dirichlet-Multinomial EB
#' Test differential exon usage with a VGAM Dirichlet-multinomial model using empirical Bayes fixed precision.
#' @param dd A data frame of exon count data.
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
#' result <- grase::test_model_multinomial_vgam_EB(grouped_eb[[1]])
#' result
#' }
test_model_multinomial_vgam_EB <- function(dd) {
    gene <- unique(dd$gene); event <- unique(dd$event)
    log_prec_mod <- dd$log_prec_mod[1]
    rho_mod <- dd$rho_mod[1]

    wide_df <- dd %>% dplyr::select(sample, groups, type, count) %>% pivot_wider(names_from = type, values_from = count, values_fill = 0)
    Y <- as.matrix(wide_df[, setdiff(names(wide_df), c("sample", "groups"))])
    wide_df <- wide_df[rowSums(Y) > 0, ]
    Y <- Y[rowSums(Y) > 0, , drop = FALSE]
    
    if (nrow(Y) < 2 || length(unique(wide_df$groups)) < 2) return(NULL)
    if (sum(colSums(Y) > 0) < 2) return(NULL)

    # --- Fix rho using Offset and Constraints ---
    K <- ncol(Y) # Number of predictors = Number of categories (K-1 probs + 1 rho)

    # Offset: Set the K-th predictor (log precision) directly
    off <- matrix(0, nrow = nrow(Y), ncol = K)
    off[, K] <- log_prec_mod
        
    # Constraints: Suppress estimation of the K-th predictor
    # Create a K x (K-1) matrix (Identity with last column removed)
    cm <- diag(K)[, -K, drop = FALSE]
    
    # Apply constraints to Intercept and groups
    clist_full <- list("(Intercept)" = cm, groups = cm)
    clist_null <- list("(Intercept)" = cm)

    m1 <- tryCatch(
        vglm(Y ~ groups, dirmultinomial, data = wide_df, 
             offset = off, constraints = clist_full), 
        error = function(e) NULL
    )
    m0 <- tryCatch(
        vglm(Y ~ 1, dirmultinomial, data = wide_df, 
             offset = off, constraints = clist_null), 
        error = function(e) NULL
    )

    if (!is.null(m1) && !is.null(m0)) {
        LR <- 2 * (logLik(m1) - logLik(m0))
        df_diff <- df.residual(m0) - df.residual(m1)
        
        cfs <- coef(m1)
        eff_size <- cfs[grep("groups", names(cfs))][1]
        
        return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pchisq(as.numeric(LR), df_diff, lower.tail=F), model="dirmult_moderated", effect_size=eff_size))
    }
    return(NULL)
}



# Updated function for Wald tests in a Multinomial/Dirichlet-Multinomial setting
#' Test differential exon usage with Wald tests in a VGAM Dirichlet-multinomial model using empirical Bayes precision.
#' @param dd A data frame of exon count data.
#' @param shape0 Numeric. Shape parameter of the Gamma prior on precision (currently unused in the
#'   fixed-precision implementation; retained for API compatibility).
#' @param rate0 Numeric. Rate parameter of the Gamma prior on precision (currently unused in the
#'   fixed-precision implementation; retained for API compatibility).
#' @param ref_option Character string or \code{NULL}. Name of the reference isoform/exon-part
#'   category. If \code{NULL} or not found in the data, the first available option is used.
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
#' result <- grase::test_model_multinomial_vgam_wald_EB(grouped_eb[[1]],
#'                                                       shape0 = 1, rate0 = 1)
#' result
#' }
test_model_multinomial_vgam_wald_EB <- function(dd, shape0, rate0, ref_option = NULL) {
    gene <- unique(dd$gene); event <- unique(dd$event)
    prec_mod <- dd$prec_mod[1]
    log_prec_mod <- dd$log_prec_mod[1]
  
    # 1. Pivot data to wide format for multinomial modeling
    wide_df <- dd %>% 
      dplyr::select(sample, groups, type, count) %>% 
      pivot_wider(names_from = type, values_from = count, values_fill = 0)
    
    # Identify available options (isoforms/exon parts)
    all_options <- setdiff(names(wide_df), c("sample", "groups"))
    
    # 2. Set the Reference Option
    # If no ref is provided, use the first one available
    if (is.null(ref_option) || !(ref_option %in% all_options)) {
        ref_option <- all_options[1]
    }
    
    # Reorder columns so the reference is the LAST column (VGAM baseline default)
    other_options <- setdiff(all_options, ref_option)
    Y <- as.matrix(wide_df[, c(other_options, ref_option)])
    
    # Filter out empty samples
    valid_rows <- rowSums(Y) > 0
    wide_df <- wide_df[valid_rows, ]
    Y <- Y[valid_rows, , drop = FALSE]
    
    if (nrow(Y) < 2 || length(unique(wide_df$groups)) < 2 || ncol(Y) < 2) return(NULL)
    if (sum(colSums(Y) > 0) < 2) return(NULL)

    # 3. Fix rho using Offset and Constraints
    K <- ncol(Y)
    
    off <- matrix(0, nrow = nrow(Y), ncol = K)
    off[, K] <- log_prec_mod
 
    cm <- diag(K)[, -K, drop = FALSE]
    clist_full <- list("(Intercept)" = cm, groups = cm)

    # 4. Fit the Full Model with moderated precision
    m1 <- tryCatch(
      vglm(Y ~ groups, dirmultinomial, data = wide_df, 
           offset = off, constraints = clist_full), 
      error = function(e) NULL
    )

    if (!is.null(m1)) {
        # 5. Extract Wald Statistics for the 'groups' effect
        # The summary object contains the p-values for each coefficient
        s1 <- summary(m1)
        coef_table <- as.data.frame(s1@coef3)
        
        # Filter for the predictor effect (groups) and exclude intercepts
        # VGAM names these like 'groupsCD8T:1', 'groupsCD8T:2', etc.
        wald_results <- coef_table %>%
          mutate(coef_name = rownames(coef_table)) %>%
          filter(grepl("groups", coef_name)) %>%
          mutate(
            gene = gene,
            event = event,
            reference = ref_option,
            # Map the index (:1, :2) back to the actual option name
            comparison_option = other_options[as.numeric(str_extract(coef_name, "\\d+$"))]
          ) %>%
          dplyr::select(gene, event, comparison_option, reference, 
                        effect_size = Estimate, std_err = `Std. Error`, 
                        z_value = `z value`, p.value = `Pr(>|z|)`)
        
        return(wald_results)
    }
    return(NULL)
}
