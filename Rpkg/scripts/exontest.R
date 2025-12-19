library(parallel)
library(MASS)
library(matrixStats)
library(glmmTMB)
library(VGAM)
library(tidyverse)
library(grase)


  
group_by_event <- function(dat, col_y, col_n) {

  dat$y <- dat[[col_y]]                                     # 253105
  dat$n <- dat[[col_n]]
  dat <- dat[!is.na(dat$n),]                                # 20484
  dat <- dat %>% add_count(gene, event, name="n_samples")
  dat <- dat %>% filter(n_samples > 4)                      # 19157 
  # split data by (gene, event)                              
  grouped_data <- dat %>%                                   # 2469 gene_events
    group_by(gene, event) %>%
    group_split()

  return(grouped_data)
}



phi_estimator <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]

    if (nrow(dd) < 2) return(NULL)

    m0 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd, family = binomial(link = "logit")),
      error = function(e) NULL
    )
    m1 <- tryCatch(
      glmmTMB(cbind(y, n - y) ~ 1, data = dd, family = betabinomial(link = "logit")),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0)) ) {
      test <- tryCatch(anova(m1, m0, test = "LRT"), error = function(e) NULL)
      if (!is.null(test) && test$`Pr(>Chisq)`[2] < 0.05) {
        return(data.frame(gene = gene, event = event, phi = sigma(m1)))
      }
    }
    return(NULL)
}







test_model <- function(dd, prior_disp) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    model_type = NULL

    if (nrow(dd) < 2) return(NULL)
    if (length(unique(dd$groups)) < 2) return(NULL)
   
    print(dd)
 
    # full: add group
    m1 <- tryCatch (
      glmmTMB(
        cbind(y, n - y) ~ groups,
        data   = dd,
        family = betabinomial(link="logit"),
        priors = prior_disp
      ), 
      error = function(e) NULL
    )

    # null: intercept only, with φ prior
    m0 <- tryCatch ( 
      glmmTMB(
        cbind(y, n - y) ~ 1,
        data   = dd,
        family = betabinomial(link="logit"),
        priors = prior_disp
      ),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0)) ) {
      if (all(!is.na(summary(m1)$coefficients$cond[, "Std. Error"])) && all(!is.na(summary(m0)$coefficients$cond[, "Std. Error"]))) {
        model_type <- "betabinomial"
      }
      else {
        return (NULL)
      }
    } else {
      # full: add group
      m1 <- tryCatch (
        glmmTMB(
          cbind(y, n - y) ~ groups,
          data   = dd,
          family = binomial(link="logit")
        ), 
        error = function(e) NULL
      )

      # null: intercept only, with φ prior
      m0 <- tryCatch ( 
        glmmTMB(
          cbind(y, n - y) ~ 1,
          data   = dd,
          family = binomial(link="logit")
        ),
        error = function(e) NULL
      )

      if (!is.null(m1) && !is.na(logLik(m1)) && !is.null(m0) && !is.na(logLik(m0)) ) {
        if (all(!is.na(summary(m1)$coefficients$cond[, "Std. Error"])) && all(!is.na(summary(m0)$coefficients$cond[, "Std. Error"]))) {
          model_type <- "binomial"
        }
        else {
          return (NULL)
        }
      }
      else {
        return (NULL)
      }
    }

    # LRT statistic & p‐value
    LR   <- 2 * (logLik(m1) - logLik(m0))
    pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
    
    # extract MAP φ and group slope
    m0_beta0_hat <- fixef(m0)$cond[1]
    m1_beta0_hat <- fixef(m1)$cond[1]
    beta1_hat <- fixef(m1)$cond[2]
    phi_hat  <- sigma(m1)
    
    return (data.frame(
        gene = gene,
        event = event,
        LRT = as.numeric(LR),
        p.value = pval,
        model = model_type,
        m0_beta0 = m0_beta0_hat,
        m1_beta0 = m1_beta0_hat,
        beta1 = beta1_hat,
        phi = phi_hat
      )
    )
}

test_model_vgam <- function(dd) {
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    dd <- dd[dd$n > 0, ]
    model_type = NULL

    if (nrow(dd) < 2) return(NULL)
    if (length(unique(dd$groups)) < 2) return(NULL)
   
    # full: add group
    m1 <- tryCatch (
      vglm(cbind(y, n - y) ~ groups, family = betabinomial, data = dd), 
      error = function(e) NULL
    )

    # null: intercept only
    m0 <- tryCatch ( 
      vglm(cbind(y, n - y) ~ 1, family = betabinomial, data = dd),
      error = function(e) NULL
    )

    if (!is.null(m1) && !is.null(m0)) {
        model_type <- "betabinomial"
    } else {
      # full: add group
      m1 <- tryCatch (
        vglm(cbind(y, n - y) ~ groups, family = binomialff, data = dd), 
        error = function(e) NULL
      )

      # null: intercept only
      m0 <- tryCatch ( 
        vglm(cbind(y, n - y) ~ 1, family = binomialff, data = dd),
        error = function(e) NULL
      )

      if (!is.null(m1) && !is.null(m0)) {
          model_type <- "binomial"
      }
      else {
        return (NULL)
      }
    }

    # LRT statistic & p‐value
    LR   <- 2 * (logLik(m1) - logLik(m0))
    pval <- pchisq(as.numeric(LR), df = 1, lower.tail = FALSE)
    
    # extract coefficients
    cf0 <- coef(m0)
    cf1 <- coef(m1)
    
    if (model_type == "betabinomial") {
        # VGAM betabinomial: (Intercept):1 is logit(mu), (Intercept):2 is logit(rho)
        m0_beta0_hat <- cf0["(Intercept):1"]
        m1_beta0_hat <- cf1["(Intercept):1"]
        
        # Find group coefficient (e.g. groupsCD8T:1)
        beta1_idx <- grep("groups.*:1", names(cf1))
        beta1_hat <- if(length(beta1_idx)) cf1[beta1_idx] else NA
        
        # Calculate phi from rho: phi = (1/rho) - 1
        rho_logit <- cf1["(Intercept):2"]
        rho <- 1 / (1 + exp(-rho_logit))
        phi_hat <- (1/rho) - 1
    } else {
        # Binomial
        m0_beta0_hat <- cf0["(Intercept)"]
        m1_beta0_hat <- cf1["(Intercept)"]
        beta1_idx <- grep("groups", names(cf1))
        beta1_hat <- if(length(beta1_idx)) cf1[beta1_idx] else NA
        phi_hat <- NA
    }
    
    return (data.frame(
        gene = gene,
        event = event,
        LRT = as.numeric(LR),
        p.value = pval,
        model = model_type,
        m0_beta0 = as.numeric(m0_beta0_hat),
        m1_beta0 = as.numeric(m1_beta0_hat),
        beta1 = as.numeric(beta1_hat),
        phi = as.numeric(phi_hat)
      )
    )
}

test_model_multinomial_vgam <- function(dd) {
    # dd is expected to be in long format with columns: gene, event, sample, groups, type, count
    gene  <- unique(dd$gene)
    event <- unique(dd$event)
    
    # Pivot to wide format to get count matrix Y
    # Columns: sample, groups, type1, type2, ...
    wide_df <- dd %>%
      dplyr::select(sample, groups, type, count) %>%
      pivot_wider(names_from = type, values_from = count, values_fill = 0)
    
    # Identify count columns (all columns except sample and groups)
    count_cols <- setdiff(names(wide_df), c("sample", "groups"))
    Y <- as.matrix(wide_df[, count_cols])
    
    # Filter samples with 0 total counts
    keep <- rowSums(Y) > 0
    Y <- Y[keep, , drop = FALSE]
    wide_df <- wide_df[keep, ]
    
    if (nrow(Y) < 2 || length(unique(wide_df$groups)) < 2) return(NULL)
    
    model_type <- NULL
    
    # Try Dirichlet-Multinomial
    m1 <- tryCatch(
      vglm(Y ~ groups, dirmultinomial, data = wide_df),
      error = function(e) NULL
    )
    m0 <- tryCatch(
      vglm(Y ~ 1, dirmultinomial, data = wide_df),
      error = function(e) NULL
    )
    
    if (!is.null(m1) && !is.null(m0)) {
      model_type <- "dirmultinomial"
    } else {
      # Fallback to Multinomial
      m1 <- tryCatch(
        vglm(Y ~ groups, multinomial, data = wide_df),
        error = function(e) NULL
      )
      m0 <- tryCatch(
        vglm(Y ~ 1, multinomial, data = wide_df),
        error = function(e) NULL
      )
      
      if (!is.null(m1) && !is.null(m0)) {
        model_type <- "multinomial"
      } else {
        return(NULL)
      }
    }
    
    # Likelihood Ratio Test
    LR <- 2 * (logLik(m1) - logLik(m0))
    df_diff <- df.residual(m0) - df.residual(m1)
    pval <- pchisq(as.numeric(LR), df = df_diff, lower.tail = FALSE)
    
    return(data.frame(
      gene = gene,
      event = event,
      LRT = as.numeric(LR),
      p.value = pval,
      model = model_type,
      df = df_diff
    ))
}

outdir = '~/graphml.dexseq.v34/dice_exoncnts'
# exoncnt_master is generated by concatenating all internal .exoncnt.txt files and prefixing the header.
# cat ENSG*.exoncnt.txt | grep -v ref_ex_part > internal.exoncnt.txt
exoncnt_master <- paste0(outdir,'/internal.exoncnt.txt')

cond1 = 'B'
cond2 = 'CD8T'

if (file.exists(exoncnt_master)) {
  bipartitioncnts <- read.table(exoncnt_master, header=TRUE, row.names=NULL)
  bipartitioncnts$groups <- factor(bipartitioncnts$groups, levels = c(cond1, cond2))
}


grouped_data <- group_by_event(bipartitioncnts, 'diff', 'n')

if (file.exists(paste0(outdir,'/phi.txt'))) {
  phi_df <- read.table(paste0(outdir,'/phi.txt'), header=TRUE, row.names=NULL)
} else {

  cl <- makeCluster(20)  # or use a fixed number like makeCluster(4)
  clusterEvalQ(cl, library(glmmTMB))
  clusterExport(cl, varlist = c("phi_estimator", "grouped_data"), envir = environment())

  phi_list <- parLapply(cl, grouped_data, phi_estimator)
  #for (i in seq_along(grouped_data)) {
  # phi_list[[i]] <- phi_estimator(grouped_data[[i]])
  #}
  
  stopCluster(cl)
  phi_df <- bind_rows(Filter(Negate(is.null), phi_list))
  write.table(phi_df, file=paste0(outdir,'/phi.txt'), quote=FALSE, sep="\t") 
}

phi_trimmed <- phi_df$phi[phi_df$phi < 1e+10]
fit_phi <- fitdistr(phi_trimmed, "gamma")
shape0  <- fit_phi$estimate["shape"]
rate0   <- fit_phi$estimate["rate"]

prior_disp <- data.frame(
  prior = sprintf("gamma(%g,%g)", shape0, rate0),
  class = "fixef_disp",
  coef  = ""            # blank all dispersion parameters
)
print(prior_disp)

grouped_data <- grouped_data %>% keep (~ all(.x$n_samples == 8)) # 2469->2138
#grouped_data2 <- grouped_data2 %>% keep (~ all(.x$n_samples == 8)) # 30624->25970

cl <- makeCluster(20)  # or use a fixed number like makeCluster(4)
clusterEvalQ(cl, library(glmmTMB))
clusterExport(cl, varlist = c("test_model", "grouped_data", "prior_disp"), envir = environment())

result_list <- parLapply(cl, grouped_data, function(dd) {
  test_model(
    dd        = dd,
    prior_disp      = prior_disp
  )
})

stopCluster(cl)
results <- bind_rows(Filter(Negate(is.null), result_list))


combined_df <- bind_rows(grouped_data)

wide_df <- combined_df %>% dplyr::select(gene, event, source, sink, ref_ex_part, setdiff, sample, ref, diff)  %>%
  pivot_wider(
    names_from = sample,
    values_from = c(diff, ref)
  )
combined_results <- left_join(wide_df, results, by = c("gene", "event"))
combined_results$baseMean <- rowMeans(combined_results[, grep("diff_sample", names(combined_results))])
combined_results$pvalue <- combined_results$p.value 
combined_results <- pvalueAdjustment(combined_results, independentFiltering=TRUE, alpha = 0.1, pAdjustMethod = 'BH')
#combined_results$FDR <- p.adjust(combined_results$p.value, method = "BH")
write.table(combined_results, file=paste0(outdir,'/test.txt'), quote=FALSE, sep="\t")



