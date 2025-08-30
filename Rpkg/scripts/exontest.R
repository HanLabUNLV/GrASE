library(parallel)
library(MASS)
library(matrixStats)
library(glmmTMB)
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





estimate_phi <- function (grouped_data) {

  cl <- makeCluster(10)  # or use a fixed number like makeCluster(4)
  clusterEvalQ(cl, library(glmmTMB))
  clusterExport(cl, varlist = c("phi_estimator", "grouped_data"), envir = environment())

  phi_list <- parLapply(cl, grouped_data, phi_estimator)
  #for (i in seq_along(grouped_data)) {
  # phi_list[[i]] <- phi_estimator(grouped_data[[i]])
  #}
  
  stopCluster(cl)
  phi_df <- bind_rows(Filter(Negate(is.null), phi_list))

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

run_test_model <- function (grouped_data, prior_disp, cl_num) {

  cl <- makeCluster(cl_num)  # or use a fixed number like makeCluster(4)
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

}



outdir = '~/graphml.dexseq.v34/exoncnts.filtered'
exoncnt_master <- paste0(outdir,'/exoncnt.nonTSS.txt')

cond1 = 'B'
cond2 = 'CD8T'

if (file.exists(exoncnt_master)) {
  focalexoncnts <- read.table(exoncnt_master, header=TRUE, row.names=NULL)
  focalexoncnts$groups <- factor(focalexoncnts$groups, levels = c(cond1, cond2))
}


grouped_data <- group_by_event(focalexoncnts, 'diff', 'n')

if (file.exists(paste0(outdir,'/phi.txt'))) {
  phi_df <- read.table(paste0(outdir,'/phi.txt'), header=TRUE, row.names=NULL)
} else {
  phi_df <- estimate_phi(grouped_data)
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
results <- run_test_model(grouped_data, prior_disp, 20)
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



