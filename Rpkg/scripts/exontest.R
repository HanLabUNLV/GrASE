library(parallel)
library(MASS)
library(matrixStats)
library(glmmTMB)
library(VGAM)
library(tidyverse)
library(grase)
library(optparse)

# --- Main Execution Script ---

indir = '~/GrASE_simulation/bipartition.test'
model = 'glmmTMB_prior'
masterfile = 'bipartition.exoncnt.combined.txt'
phifile = 'phi.glmmtmb.txt'
cond1 = 'group1'; cond2 = 'group2'
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-f", "--file"), type="character",  
              help="exoncnt dataset file name", metavar="character"),
  make_option(c("-d", "--dir"), type="character", 
              help="input/output dir name", metavar="character"),
  make_option(c("-p", "--phi"), type="character", 
              help="phi estimate filename", metavar="character"),
  make_option(c("-m", "--model"), type="character", 
              help="model name", metavar="character"),
  make_option(c("--cond1"), type="character", 
              help="condition 1 name (Group 1 - Group 2)", metavar="character"),
  make_option(c("--cond2"), type="character", 
              help="condition 2 name (Reference Group)", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)
if (!is.null(opt$dir)) {
  indir = opt$dir
}
if (!is.null(opt$file)) {
  masterfile = opt$file
}
if (!is.null(opt$phi)) {
  phifile = opt$phi
}
if (!is.null(opt$model)) {
  model = opt$model 
}
if (!is.null(opt$cond1)) {
  cond1 = opt$cond1 
}
if (!is.null(opt$cond2)) {
  cond2 = opt$cond2 
}
outdir = indir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
exoncnt_master <- paste0(indir,'/', masterfile)


if (file.exists(exoncnt_master)) {
  splitcnts <- read.table(exoncnt_master, header=TRUE, row.names=NULL)
  # Set levels such that cond2 is reference, so coefficient is cond1 - cond2
  splitcnts$groups <- factor(splitcnts$groups, levels = c(cond2, cond1))
} else {
  print('exoncnt dataset file does not exist.')
  print('please combine the exoncounts into one file and specify the filename.')
  stop()
}

grouped_data <- group_by_event(splitcnts, 'diff', 'n')
#grouped_data <- grouped_data[1:10]

# Global Precision Estimation 
if (model == 'glmmTMB_prior' || model == 'glmmTMB_fixedEB' || model == 'VGAM_MLE_EB_init') {

  if (file.exists(paste0(outdir,'/', phifile))) {
    phi_df <- read.table(paste0(outdir,'/', phifile), header=TRUE, row.names=NULL)
   } else {
    #cl <- makeCluster(30, outfile='phi.glmmtmb.log')
    #clusterEvalQ(cl, library(glmmTMB))
    #clusterExport(cl, varlist = c("phi_estimate_glmmTMB", "grouped_data"), envir = environment())
    phi_list <- mclapply(grouped_data, function(dd) {
    #phi_list <- parLapply(cl, grouped_data, function(dd) {
      tryCatch({
        phi_estimate_glmmTMB(dd) 
      }, error = function(e) {
        msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       dd$gene[1], Sys.getpid(), conditionMessage(e))
        cat(msg, file = "phi.glmmtmb.errors.log", append = TRUE)
        NULL
      })
    #})
    }, mc.cores = 32)
    #stopCluster(cl)
    phi_df <- bind_rows(Filter(Negate(is.null), phi_list))
    write.table(phi_df, file=paste0(outdir,'/', phifile), quote=FALSE, sep="\t", row.names = FALSE) 
  }
} else if (model == 'DM_EB' || model == 'DM_Wald_EB') {
  # DM Moderation Pipeline (NEW)
  # group_by_event on counts instead of diff/n for DM
  prec_file = 'prec_dm.txt'
  if (file.exists(paste0(outdir,'/', prec_file))) {
    prec_table <- read.table(paste0(outdir,'/', prec_file), header=TRUE, row.names=NULL)
  } else {
    grouped_counts <- splitcnts %>% group_by(gene, event) %>% group_split()
    #cl <- makeCluster(30); clusterEvalQ(cl, {library(VGAM); library(tidyverse)})
    #clusterExport(cl, "prec_estimate_vgam")
    prec_list <- mclapply(grouped_counts, prec_estimate_vgam, mc.cores=32) 
    #prec_list <- parLapply(cl, grouped_counts, prec_estimate_vgam)
    #stopCluster(cl)
    prec_table <- moderate_prec_log_scale(bind_rows(Filter(Negate(is.null), prec_list)))
    write.table(prec_table, paste0(outdir,'/', prec_file), sep="\t", quote=F, row.names = FALSE)
  }
}

# Per Gene Model Estimation 
if (model == 'glmmTMB_prior') {
  phi_trimmed <- phi_df[phi_df$phi < 1e+10 & phi_df$phi > 0,'phi']

  # fit normal to log(phi)
  log_phi_vals <- log(phi_trimmed)
  fit_logphi <- fitdistr(log_phi_vals, "normal")
  mean_logphi <- fit_logphi$estimate["mean"]
  sd_logphi   <- fit_logphi$estimate["sd"]
  prior_disp <- data.frame(prior = sprintf("normal(%g,%g)", mean_logphi, sd_logphi), 
    class = "fixef_disp", 
    coef = "", # Explicitly target the intercept 
    stringsAsFactors = FALSE
  )
  print(prior_disp)


  #cl <- makeCluster(40, outfile='testglmmTMB_withprior.log')
  #clusterEvalQ(cl, { library(glmmTMB); library(VGAM); library(tidyverse) })
  #clusterExport(cl, varlist = c("test_model_glmmTMB_with_prior", "prior_disp"), envir = environment())

  result_list <- mclapply(grouped_data, function(dd) {
  #result_list <- parLapply(cl, grouped_data, function(dd) {
    tryCatch({
      # Suppress glmmTMB timing messages from failed optimizations
      suppressMessages(test_model_glmmTMB_with_prior(dd, prior_disp))
    }, error = function(e) {
        msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       dd$gene[1], Sys.getpid(), conditionMessage(e))
        cat(msg, file = "glmmtmb_withprior.errors.log", append = TRUE)
        NULL
    })
  #})
  }, mc.cores = 32)
  #stopCluster(cl)

  results <- bind_rows(Filter(Negate(is.null), result_list))
  
  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=FALSE, theta=NULL, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }
  
  write.table(results, file=paste0(outdir,'/test_glmmTMB_MAP_prior.txt'), quote=FALSE, sep="\t", row.names = FALSE)

} else if (model == 'glmmTMB_fixedEB') {

  ## Empirical Bayes hyperparameters: prior mean and variance
  phi_table <- moderate_phi_log_scale(phi_df)
  splitcnts_eb <- splitcnts %>%
    left_join(phi_table, by = c("gene", "event"))
  grouped_data_eb <- group_by_event(splitcnts_eb, 'diff', 'n')
  #grouped_data_eb <- grouped_data_eb[1:10]

  # Moderated glmmTMB with EB
  #cl <- makeCluster(40, outfile='testglmmTMB_EB.log')
  #clusterEvalQ(cl, { library(glmmTMB); library(VGAM); library(tidyverse) })
  #clusterExport(cl, varlist = c("test_model_glmmTMB_EB"), envir = environment())

  result_list <- mclapply(grouped_data_eb, function(dd) {
  #result_list <- parLapply(cl, grouped_data_eb, function(dd) {
    tryCatch({
      test_model_glmmTMB_EB(dd) # Defaults to moderated glmmTMB BB
    }, error = function(e) {
        msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       dd$gene[1], Sys.getpid(), conditionMessage(e))
        cat(msg, file = "glmmtmb_EB.errors.log", append = TRUE)
        NULL
    })
  #})
  }, mc.cores = 32)
  #stopCluster(cl)

  results <- bind_rows(Filter(Negate(is.null), result_list))
  
  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=FALSE, theta=NULL, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }

  write.table(results, file=paste0(outdir,'/test_glmmTMB_fixed_EB.txt'), quote=FALSE, sep="\t", row.names=FALSE)

} else if (model == 'VGAM_MLE_EB_init') {

  ## Empirical Bayes hyperparameters: prior mean and variance
  phi_table <- moderate_phi_log_scale(phi_df)
  splitcnts_eb <- splitcnts %>%
    left_join(phi_table, by = c("gene", "event"))
  grouped_data_eb <- group_by_event(splitcnts_eb, 'diff', 'n')
  #grouped_data_eb <- grouped_data_eb[1:10]


  # Moderated VGAM with EB
  #cl <- makeCluster(40, outfile='testVGAM_EB.log')
  #clusterEvalQ(cl, { library(glmmTMB); library(VGAM); library(tidyverse) })
  #clusterExport(cl, varlist = c("test_model_vgam_EB_init"), envir = environment())

  result_list <- mclapply(grouped_data_eb, function(dd) {
  #result_list <- parLapply(cl, grouped_data_eb, function(dd) {
    tryCatch({
      test_model_vgam_EB_init(dd) 
    }, error = function(e) {
        msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       dd$gene[1], Sys.getpid(), conditionMessage(e))
        cat(msg, file = "vgam_eb.errors.log", append = TRUE)
        NULL
    })
  #})
  }, mc.cores = 32)
  #stopCluster(cl)
  results <- bind_rows(Filter(Negate(is.null), result_list))
  
  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=FALSE, theta=NULL, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }
  
  write.table(results, file=paste0(outdir,'/test_VGAM_MLE_EB_init.txt'), quote=FALSE, sep="\t", row.names = FALSE)


} else if (model == 'DM_EB') {

  splitcnts_eb <- splitcnts %>%
    left_join(prec_table, by = c("gene", "event"))
  grouped_eb <- splitcnts_eb %>% group_by(gene, event) %>% group_split()
  #grouped_eb <- group_by_event(splitcnts_eb, 'diff', 'n')

  # Run DM EB
  #cl <- makeCluster(40); clusterEvalQ(cl, {library(VGAM); library(tidyverse)})
  #clusterExport(cl, "test_model_multinomial_vgam_EB")
  res_dm <- mclapply(grouped_eb, function(dd) {
  #res_dm <- parLapply(cl, grouped_eb, function(dd) {
    tryCatch(
      test_model_multinomial_vgam_EB(dd), error=function(e) NULL
    )
  }, mc.cores = 32)
  #stopCluster(cl)
  results <- bind_rows(res_dm)
  
  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=FALSE, theta=NULL, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }

  write.table(results, paste0(outdir, "/test_DM_EB.txt"), sep="\t", quote=F, row.names = FALSE)

} else if (model == 'DM_Wald_EB') {                                                                                            

  # DM Moderation Pipeline (Same as DM_EB)
  prec_file = 'prec_dm.txt'
  if (file.exists(paste0(outdir,'/', prec_file))) {
    prec_table <- read.table(paste0(outdir,'/', prec_file), header=TRUE, row.names=NULL)
  } else {
    grouped_counts <- splitcnts %>% group_by(gene, event) %>% group_split()
    #cl <- makeCluster(30); clusterEvalQ(cl, {library(VGAM); library(tidyverse)})
    #clusterExport(cl, "prec_estimate_vgam")
    prec_list <- mclapply(grouped_counts, prec_estimate_vgam, mc.cores = 32)
    #prec_list <- parLapply(cl, grouped_counts, prec_estimate_vgam)
    #stopCluster(cl)
    prec_table <- moderate_prec_log_scale(bind_rows(Filter(Negate(is.null), prec_list)))
    write.table(prec_table, paste0(outdir,'/', prec_file), sep="\t", quote=F, row.names = FALSE)
  }

  splitcnts_eb <- splitcnts %>%
    left_join(prec_table, by = c("gene", "event"))
  grouped_eb <- splitcnts_eb %>% group_by(gene, event) %>% group_split()

  # Run DM Wald EB
  #cl <- makeCluster(40); clusterEvalQ(cl, {library(VGAM); library(tidyverse)})
  #clusterExport(cl, "test_model_multinomial_vgam_wald_EB")
  res_dm <- mclapply(grouped_eb, function(dd) {
  #res_wald <- parLapply(cl, grouped_eb, function(dd) {
    tryCatch({
      test_model_multinomial_vgam_wald_EB(dd) 
    }, error = function(e) NULL)
  #})
  }, mc.cores = 32)
  #stopCluster(cl)
  wald_results_df <- bind_rows(Filter(Negate(is.null), res_dm))
  
  if (exists("pvalueAdjustment") && nrow(wald_results_df) > 0) {
      if (is.null(wald_results_df$pvalue)) wald_results_df$pvalue <- wald_results_df$p.value
      wald_results_df <- pvalueAdjustment(wald_results_df, independentFiltering=FALSE, theta=NULL, alpha=0.05, pAdjustMethod="BH")
      # wald_results_df$pvalue <- NULL # Keep if needed or remove. Wald test might have p.value or Pr(>Chi)
  }

  write.table(wald_results_df, file=paste0(outdir, "/test_DM_Wald_EB.txt"), sep="\t", quote=FALSE, row.names = FALSE)
} else if (model == 'wilcoxon') {

  result_list <- mclapply(grouped_data, function(dd) {
    tryCatch({
      test_model_wilcoxon(dd) 
    }, error = function(e) {
        msg <- sprintf("[%s] ERROR in %s (PID %d): %s\n",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       dd$gene[1], Sys.getpid(), conditionMessage(e))
        cat(msg, file = "wilcoxon.errors.log", append = TRUE)
        NULL
    })
  }, mc.cores = 32)

  results <- bind_rows(Filter(Negate(is.null), result_list))
  
  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=FALSE, theta=NULL, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }

  write.table(results, file=paste0(outdir,'/test_wilcoxon.txt'), quote=FALSE, sep="\t", row.names = FALSE)

} else {
  print ("model misspecified")
  print ("choose between glmmTMB_prior, glmmTMB_fixedEB, VGAM_MLE_EB_init, DM_EB, DM_Wald_EB, wilcoxon")

}
