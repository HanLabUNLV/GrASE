# --- Data Preparation Functions ---

#' @export
group_by_event <- function(dat, col_y, col_n) {
  dat$y <- dat[[col_y]]
  dat$n <- dat[[col_n]]
  dat <- dat[!is.na(dat$n),]
  dat <- dat %>% add_count(gene, event, name="n_samples")
  dat <- dat %>% filter(n_samples > 4)
  grouped_data <- dat %>%
    group_by(gene, event) %>%
    group_split()
  return(grouped_data)
}

#' @export
phi_estimate_glmmTMB <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2) return(NULL)

    m0 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    m1 <- tryCatch(
      glmmTMB(
        cbind(y, n - y) ~ 1,
        data = dd,
        family = glmmTMB::betabinomial(link = "logit")
      ),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0)) ) {
      test <- tryCatch(anova(m1, m0, test = "LRT"), error = function(e) NULL)
      if (!is.null(test) && test$`Pr(>Chisq)`[2] < 0.05) {
        phi_hat <- as.numeric(sigma(m1))  # dispersion
        vc <- vcov(m1, full = TRUE)
        ## log(phi and its variance from dispersion model
        if ( nrow(vc) < 2 || ncol(vc) < 2) {
          var_phi = NA_real_
        } else {
          var_log_phi <- vc["disp~(Intercept)", "disp~(Intercept)"]
          se_log_phi  <- sqrt(var_log_phi)

          ## Delta method: Var(phi) = phi^2 * Var(log(phi))
          var_phi <- (phi_hat^2) * var_log_phi
        }
        return(data.frame(gene = gene, event = event, phi = phi_hat, var_phi=var_phi))
      }
    }
    return(NULL)
}


## per-gene/event phi_hat and var_phi from VGAM
# DO NOT USE THIS FUNCTION
# VGAM's delta method over 1/rho^2 will always exaggerate variance when phi is large.
# converting VGAM variance to log(phi) scale before using delta method doesn't fix the issue of unstable estimate of variance
#' @export
bb_phi_estimate_VGAM <- function(dd) {

  gene  <- unique(dd$gene)
  event <- unique(dd$event)
  dd <- dd[dd$n > 0, ]
  if (nrow(dd) < 2) return(NULL)

  m_init <- tryCatch(
    vglm(cbind(y, n - y) ~ 1, family = VGAM::betabinomial, data = dd),
    error = function(e) NULL
  )
  if (is.null(m_init)) return(NULL)

  cf <- Coef(m_init)
  rho_hat <- as.numeric(cf[2])
  phi_hat <- (1 / rho_hat) - 1    #phi on original scale

  vc <- tryCatch(vcov(m_init), error = function(e) NULL)
  if (is.null(vc) || nrow(vc) < 2 || ncol(vc) < 2) return(NULL)

  var_rho <- as.numeric(vc[2, 2])

  # --- convert variance to log(phi) scale ---
  dlogphi_drho <- -1 / (rho_hat * (1 - rho_hat))    # derivative d log(phi) / d rho
  var_log_phi  <- (dlogphi_drho^2) * var_rho       # delta method

  # --- convert to phi scale ---
  var_phi <- (phi_hat^2) * var_log_phi

  data.frame(
    gene    = gene,
    event   = event,
    phi_hat = phi_hat,
    var_phi = var_phi
  )
}


#' @export
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


# --- DM Estimation & Moderation (NEW) ---

#' @export
prec_estimate_vgam <- function(dd) {
    gene <- unique(dd$gene); event <- unique(dd$event)
    wide_df <- dd %>% dplyr::select(sample, groups, type, count) %>% pivot_wider(names_from = type, values_from = count, values_fill = 0)
    Y <- as.matrix(wide_df[, setdiff(names(wide_df), c("sample", "groups"))])
    if (nrow(Y) < 2 || ncol(Y) < 2) return(NULL)

    m_init <- tryCatch(vglm(Y ~ 1, dirmultinomial, data = wide_df), error = function(e) NULL)
    if (is.null(m_init)) return(NULL)

    log_prec_hat <- as.numeric(coef(m_init)[ncol(Y)])
    vc <- tryCatch(vcov(m_init), error = function(e) NULL)
    if (is.null(vc) || nrow(vc) < ncol(Y)) return(NULL)

    return(data.frame(gene = gene, event = event, log_prec = log_prec_hat, var_log_prec = vc[ncol(Y), ncol(Y)]))
}


#' @export
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
#' @export
test_model_glmmTMB_with_prior <- function(dd, prior_disp) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
   
    # Models include 'priors' for empirical Bayes moderation
    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ groups, data = dd, family = glmmTMB::betabinomial(link="logit"), priors = prior_disp, REML=TRUE), 
      error = function(e) NULL
    )
    m0 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd, family = glmmTMB::betabinomial(link="logit"), priors = prior_disp, REML=TRUE),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0))) {
        
        # If phi is extremely large (>1e4), the model effectively converged to Binomial.
        # We force fallback to standard GLM for better stability and interpretability.
        # this is necessary: 
        # "Hauck-Donner effect" (or similar separation issues) in the glmmTMB optimizer.
#< ENSG00000285219.2	242	51.7512574983776	6.29964641910099e-13	betabinomial_glmmTMB_MAP_with_prior22280123.2191247	-48.0759117491982 without fallback
#> ENSG00000285219.2	242	0.483132649576106	0.487006771127791	binomial_glm	NA	-16.4594777497934107889c107889 # with fallback
#< ENSG00000285219.2	248	51.7512574983776	6.29964641910099e-13	betabinomial_glmmTMB_MAP_with_prior22280123.2191247	-48.0759117491982 without fallback
#> ENSG00000285219.2	248	0.483132649576106	0.487006771127791	binomial_glm	NA	-16.4594777497934 with fallback

        if (sigma(m1) > 1e5) {
            m1 <- NULL # Trends to else block
        } else {
            LR <- 2 * (logLik(m1) - logLik(m0))
            #LR <- max(0, as.numeric(LR))
            pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
            eff_size <- fixef(m1)$cond[2]
            return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pval, model="betabinomial_glmmTMB_MAP_with_prior", phi=sigma(m1), effect_size=eff_size))
        }
    }

    if (is.null(m1) || is.na(logLik(m1)) || is.null(m0) || is.na(logLik(m0))) {
        m1 <- tryCatch(glm(cbind(y, n - y) ~ groups, data = dd, family = binomial), error = function(e) NULL)
        m0 <- tryCatch(glm(cbind(y, n - y) ~ 1, data = dd, family = binomial), error = function(e) NULL)
        if (!is.null(m1) && !is.null(m0)) {
             LR <- 2 * (logLik(m1) - logLik(m0))
             LR <- max(0, as.numeric(LR))
             pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
             eff_size <- coef(m1)[2]
             return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pval, model="binomial_glm", phi=NA, effect_size=eff_size))
        }
    }
    return(NULL)
}



# 2. Moderated glmmTMB Beta-Binomial EB
#' @export
test_model_glmmTMB_EB <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)

    # Models include 'priors' for empirical Bayes moderation
    target_log_phi = dd$z_mod[1]

    # Fallback to Binomial GLM if dispersion is missing
    if (is.na(target_log_phi)) {
        m1 <- tryCatch(glm(cbind(y, n - y) ~ groups, data = dd, family = binomial), error = function(e) NULL)
        m0 <- tryCatch(glm(cbind(y, n - y) ~ 1, data = dd, family = binomial), error = function(e) NULL)
        if (!is.null(m1) && !is.null(m0)) {
             LR <- 2 * (logLik(m1) - logLik(m0))
             pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
             eff_size <- coef(m1)[2]
             return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pval, model="binomial_glm", phi=NA, effect_size=eff_size))
        }
        return(NULL)
    }

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
#' @export
test_model_vgam_EB_init <- function(dd) {

  ## per-gene moderated LRT in VGAM
  gene  <- unique(dd$gene)
  event <- unique(dd$event)
  dd <- dd[dd$n > 0, ]
  if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
  
  phi_mod = dd$phi_mod[1]

  # Fallback to Binomial GLM if dispersion is missing
  if (is.na(phi_mod)) {
      m1 <- tryCatch(glm(cbind(y, n - y) ~ groups, data = dd, family = binomial), error = function(e) NULL)
      m0 <- tryCatch(glm(cbind(y, n - y) ~ 1, data = dd, family = binomial), error = function(e) NULL)
      if (!is.null(m1) && !is.null(m0)) {
          LR <- 2 * (logLik(m1) - logLik(m0))
          pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
          eff_size <- coef(m1)[2]
          return(data.frame(
            gene    = gene,
            event   = event,
            LRT     = as.numeric(LR),
            p.value = pval,
            phi_hat = NA,
            phi_mod = NA,
            w       = NA,
            rho     = NA,
            phi     = NA,
            model   = "binomial_glm",
            effect_size = eff_size,
            stringsAsFactors = FALSE
          ))
      }
      return(NULL)
  }

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
#' @export
test_model_wilcoxon <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    
    # Filter valid rows
    dd <- dd[dd$n > 0, ]
    
    # Ensure there are enough samples and both groups are present
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
    
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
#' @export
test_model_multinomial_vgam_EB <- function(dd) {
    gene <- unique(dd$gene); event <- unique(dd$event)
    log_prec_mod <- dd$log_prec_mod[1]
    rho_mod <- dd$rho_mod[1]

    wide_df <- dd %>% dplyr::select(sample, groups, type, count) %>% pivot_wider(names_from = type, values_from = count, values_fill = 0)
    Y <- as.matrix(wide_df[, setdiff(names(wide_df), c("sample", "groups"))])
    wide_df <- wide_df[rowSums(Y) > 0, ]
    Y <- Y[rowSums(Y) > 0, , drop = FALSE]
    
    if (nrow(Y) < 2 || length(unique(wide_df$groups)) < 2) return(NULL)

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
#' @export
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
