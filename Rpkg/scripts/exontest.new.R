library(parallel)
library(MASS)
library(matrixStats)
library(glmmTMB)
library(VGAM)
library(tidyverse)
library(grase)

# --- Data Preparation Functions ---

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

moderate_prec_log_scale <- function(prec_table) {
  z_bar <- mean(prec_table$log_prec, na.rm = TRUE)
  tau2_z <- max(var(prec_table$log_prec, na.rm = TRUE) - mean(prec_table$var_log_prec, na.rm = TRUE), 0)
  prec_table$w_prec <- tau2_z / (tau2_z + prec_table$var_log_prec)
  prec_table$log_prec_mod <- (prec_table$w_prec * prec_table$log_prec) + ((1 - prec_table$w_prec) * z_bar)
  prec_table$prec_mod <- exp(prec_table$log_prec_mod)
  return(prec_table)
}



# --- Moderated Testing Functions ---

# 1. Moderated glmmTMB Beta-Binomial with Prior
test_model_glmmTMB_with_prior <- function(dd, prior_disp) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
   
    # Models include 'priors' for empirical Bayes moderation
    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ groups, data = dd, family = glmmTMB::betabinomial(link="logit"), priors = prior_disp), 
      error = function(e) NULL
    )
    m0 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd, family = glmmTMB::betabinomial(link="logit"), priors = prior_disp),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0))) {
        LR <- 2 * (logLik(m1) - logLik(m0))
        pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
        return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pval, model="betabinomial_glmmTMB_MAP_with_prior", phi=sigma(m1)))
    }
    return(NULL)
}



# 2. Moderated glmmTMB Beta-Binomial EB
test_model_glmmTMB_EB <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)

    # Models include 'priors' for empirical Bayes moderation
    target_log_phi = dd$z_mod[1]
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
        return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pval, model="betabinomial_glmmTMB_fixed_EB", phi=sigma(m1)))
    }
    return(NULL)
}




# 3. EB-moderated VGAM beta-binomial
test_model_vgam_EB <- function(dd) {

  ## per-gene moderated LRT in VGAM
  gene  <- unique(dd$gene)
  event <- unique(dd$event)
  dd <- dd[dd$n > 0, ]
  if (nrow(dd) < 2 || length(unique(dd$groups)) < 2) return(NULL)
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

  LR <- 2 * (logLik(m1) - logLik(m0))
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
    stringsAsFactors = FALSE
  )

}



# 4. Moderated VGAM Dirichlet-Multinomial EB
test_model_multinomial_vgam_EB <- function(dd) {
    gene <- unique(dd$gene); event <- unique(dd$event)
    prec_mod <- dd$prec_mod[1]
    wide_df <- dd %>% dplyr::select(sample, groups, type, count) %>% pivot_wider(names_from = type, values_from = count, values_fill = 0)
    Y <- as.matrix(wide_df[, setdiff(names(wide_df), c("sample", "groups"))])
    wide_df <- wide_df[rowSums(Y) > 0, ]
    Y <- Y[rowSums(Y) > 0, , drop = FALSE]
    
    if (nrow(Y) < 2 || length(unique(wide_df$groups)) < 2) return(NULL)

    m1 <- tryCatch(vglm(Y ~ groups, dirmultinomial(idisp = prec_mod, zero = ncol(Y)), data = wide_df), error = function(e) NULL)
    m0 <- tryCatch(vglm(Y ~ 1, dirmultinomial(idisp = prec_mod, zero = ncol(Y)), data = wide_df), error = function(e) NULL)

    if (!is.null(m1) && !is.null(m0)) {
        LR <- 2 * (logLik(m1) - logLik(m0))
        df_diff <- df.residual(m0) - df.residual(m1)
        return(data.frame(gene=gene, event=event, LRT=as.numeric(LR), p.value=pchisq(as.numeric(LR), df_diff, lower.tail=F), model="dirmult_moderated"))
    }
    return(NULL)
}



# Updated function for Wald tests in a Multinomial/Dirichlet-Multinomial setting
test_model_multinomial_vgam_wald <- function(dd, shape0, rate0, ref_option = NULL) {
    gene <- unique(dd$gene); event <- unique(dd$event)
    prec_mod <- dd$prec_mod[1]
 
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

    # 3. Estimate precision and apply moderation (same as previous logic)
    m_init <- tryCatch(vglm(Y ~ 1, dirmultinomial, data = wide_df), error = function(e) NULL)
    if (is.null(m_init)) return(NULL)

    # 4. Fit the Full Model with moderated precision
    # 'zero = ncol(Y)' ensures precision is an intercept-only constant (constant across the group effect)
    m1 <- tryCatch(
      vglm(Y ~ groups, dirmultinomial(idisp = prec_mod, zero = ncol(Y)), data = wide_df), 
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
                        log_odds_ratio = Estimate, std_err = `Std. Error`, 
                        z_value = `z value`, p.value = `Pr(>|z|)`)
        
        return(wald_results)
    }
    return(NULL)
}



exit()
# --- Main Execution Script ---

outdir = '~/graphml.dexseq.v34/dice_exoncnts'
exoncnt_master <- paste0(outdir,'/bipartitions.internal.exoncnt.txt')
cond1 = 'B'; cond2 = 'CD8T'

if (file.exists(exoncnt_master)) {
  bipartitioncnts <- read.table(exoncnt_master, header=TRUE, row.names=NULL)
  bipartitioncnts$groups <- factor(bipartitioncnts$groups, levels = c(cond1, cond2))
}

grouped_data <- group_by_event(bipartitioncnts, 'diff', 'n')

# Global Dispersion Estimation (glmmTMB)
if (file.exists(paste0(outdir,'/phi.glmmtmb.txt'))) {
  phi_df <- read.table(paste0(outdir,'/phi.glmmtmb.txt'), header=TRUE, row.names=NULL)
} else {
  cl <- makeCluster(20, outfile='phi.glmmtmb.log')
  clusterEvalQ(cl, library(glmmTMB))
  clusterExport(cl, varlist = c("phi_estimate_glmmTMB", "grouped_data"), envir = environment())
  phi_list <- parLapply(cl, grouped_data, function(dd) {
    tryCatch({
      phi_estimate_glmmTMB(dd) 
    }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     dd$gene[1], Sys.getpid(), conditionMessage(e))
      cat(msg, file = "phi.glmmtmb.errors.log", append = TRUE)
      NULL
    })
  })
  stopCluster(cl)
  phi_df <- bind_rows(Filter(Negate(is.null), phi_list))
  write.table(phi_df, file=paste0(outdir,'/phi.glmmtmb.txt'), quote=FALSE, sep="\t") 
}


phi_trimmed <- phi_df$phi[phi_df$phi < 1e+10]

# Moderated glmmTMB with Prior 
fit_phi <- fitdistr(phi_trimmed, "gamma")
shape0  <- fit_phi$estimate["shape"]
rate0   <- fit_phi$estimate["rate"]
prior_disp <- data.frame(prior = sprintf("gamma(%g,%g)", shape0, rate0), class = "fixef_disp", coef = "")
print(prior_disp)

cl <- makeCluster(30, outfile='testglmmTMB_withprior.log')
clusterEvalQ(cl, { library(glmmTMB); library(VGAM); library(tidyverse) })
clusterExport(cl, varlist = c("test_model_glmmTMB_with_prior", "prior_disp"), envir = environment())

result_list <- parLapply(cl, grouped_data, function(dd) {
  tryCatch({
    test_model_glmmTMB_with_prior(dd, prior_disp) # Defaults to moderated glmmTMB BB
  }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     dd$gene[1], Sys.getpid(), conditionMessage(e))
      cat(msg, file = "glmmtmb_withprior.errors.log", append = TRUE)
      NULL
  })
})
stopCluster(cl)

results <- bind_rows(Filter(Negate(is.null), result_list))
write.table(results, file=paste0(outdir,'/test_glmmTMB_prior.txt'), quote=FALSE, sep="\t")


## Empirical Bayes hyperparameters: prior mean and variance
phi_table <- moderate_phi_log_scale(phi_df)
bipartitioncnts_eb <- bipartitioncnts %>%
  left_join(phi_table, by = c("gene", "event"))
grouped_data_eb <- group_by_event(bipartitioncnts_eb, 'diff', 'n')
grouped_data_eb <- Filter(function(x) !is.na(x$z_mod[1]), grouped_data_eb)

# Moderated glmmTMB with EB
cl <- makeCluster(30, outfile='testglmmTMB_withprior.log')
clusterEvalQ(cl, { library(glmmTMB); library(VGAM); library(tidyverse) })
clusterExport(cl, varlist = c("test_model_glmmTMB_EB"), envir = environment())

result_list <- parLapply(cl, grouped_data_eb, function(dd) {
  tryCatch({
    test_model_glmmTMB_EB(dd) # Defaults to moderated glmmTMB BB
  }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     dd$gene[1], Sys.getpid(), conditionMessage(e))
      cat(msg, file = "glmmtmb_withprior.errors.log", append = TRUE)
      NULL
  })
})
stopCluster(cl)

results <- bind_rows(Filter(Negate(is.null), result_list))
write.table(results, file=paste0(outdir,'/test_glmmTMB_EB.txt'), quote=FALSE, sep="\t")




# Moderated VGAM with EB
cl <- makeCluster(30, outfile='testVGAM_EB.log')
clusterEvalQ(cl, { library(glmmTMB); library(VGAM); library(tidyverse) })
clusterExport(cl, varlist = c("test_model_vgam_EB"), envir = environment())

result_list <- parLapply(cl, grouped_data_eb, function(dd) {
  tryCatch({
    test_model_vgam_EB(dd) 
  }, error = function(e) {
      msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     f, Sys.getpid(), conditionMessage(e))
      cat(msg, file = "vgam_eb.errors.log", append = TRUE)
      NULL
  })
})
stopCluster(cl)
results <- bind_rows(Filter(Negate(is.null), result_list))
write.table(results, file=paste0(outdir,'/test_VGAM_EB.txt'), quote=FALSE, sep="\t")





# DM Moderation Pipeline (NEW)
# group_by_event on counts instead of diff/n for DM
grouped_counts <- bipartitioncnts %>% group_by(gene, event) %>% group_split()
cl <- makeCluster(20); clusterEvalQ(cl, {library(VGAM); library(tidyverse)})
clusterExport(cl, "prec_estimate_vgam")
prec_list <- parLapply(cl, grouped_counts, prec_estimate_vgam)
stopCluster(cl)
prec_table <- moderate_prec_log_scale(bind_rows(Filter(Negate(is.null), prec_list)))


bipartitioncnts_eb <- bipartitioncnts %>%
  left_join(prec_table, by = c("gene", "event"))
grouped_eb <- group_by_event(bipartitioncnts_eb, 'diff', 'n')

# Run DM EB
cl <- makeCluster(30); clusterEvalQ(cl, {library(VGAM); library(tidyverse)})
clusterExport(cl, "test_model_multinomial_vgam_EB")
res_dm <- parLapply(cl, grouped_eb, function(dd) tryCatch(test_model_multinomial_vgam_EB(dd), error=function(e) NULL))
stopCluster(cl)
write.table(bind_rows(res_dm), paste0(outdir, "/test_DM_EB.txt"), sep="\t", quote=F)
                                                                                            

# Run DM Wald EB
cl <- makeCluster(30); clusterEvalQ(cl, {library(VGAM); library(tidyverse)})
clusterExport(cl, "test_model_multinomial_vgam_wald_eb")
res_wald <- parLapply(cl, grouped_data_eb, function(dd) {
  tryCatch({
    test_model_multinomial_vgam_wald_eb(dd) 
  }, error = function(e) NULL)
})
stopCluster(cl)
wald_results_df <- bind_rows(Filter(Negate(is.null), res_wald))
write.table(wald_results_df, paste0(outdir, "/test_DM_Wald_EB.txt"), sep="\t", quote=F)               

