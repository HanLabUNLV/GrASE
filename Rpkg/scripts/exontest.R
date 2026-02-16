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
countdir = '~/GrASE_simulation/bipartition.internal.counts'
split = 'bipartition'
model = 'glmmTMB_prior'
masterfile = 'bipartition.internal.exoncnt.combined.txt'
phifile = 'phi.glmmtmb.txt'
cond1 = 'group1'; cond2 = 'group2'
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-f", "--file"), type="character",  
              help="file name of combined exoncnt dataset", metavar="character"),
  make_option(c("-d", "--dir"), type="character", 
              help="input/output dir that contains the exoncnt file", metavar="character"),
  make_option(c("-c", "--countdir"), type="character", 
              help="count dir that contains the count and split files by gene", metavar="character"),
  make_option(c("-s", "--splittype"), type="character", 
              help="split type: bipartition or multinomial or n_choose_2", metavar="character"),
  make_option(c("-p", "--phi"), type="character", 
              help="filename for phi estimates", metavar="character"),
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
if (!is.null(opt$file)) {
  masterfile = opt$file
}
if (!is.null(opt$dir)) {
  indir = opt$dir
}
if (!is.null(opt$countdir)) {
  countdir = opt$countdir
}
if (!is.null(opt$splittype)) {
  split = opt$splittype
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
out_prefix <- sub("^([^.]+\\.[^.]+).*", "\\1", masterfile)
out_resultfile=paste0(outdir,'/test_', out_prefix,'_', model, '.txt')


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
  phi_trimmed <- phi_df[!is.na(phi_df$phi) & phi_df$phi < 1e+10 & phi_df$phi > 0,'phi']

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
  
  write.table(results, file=out_resultfile, quote=FALSE, sep="\t", row.names = FALSE)

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

  write.table(results, file=out_resultfile, quote=FALSE, sep="\t", row.names=FALSE)

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
  
  write.table(results, file=out_resultfile, quote=FALSE, sep="\t", row.names = FALSE)


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

  write.table(results, file=out_resultfile, sep="\t", quote=F, row.names = FALSE)

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

  write.table(wald_results_df, file=out_resultfile, sep="\t", quote=FALSE, row.names = FALSE)
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

  write.table(results, file=out_resultfile, quote=FALSE, sep="\t", row.names = FALSE)

} else {
  print ("model misspecified")
  print ("choose between glmmTMB_prior, glmmTMB_fixedEB, VGAM_MLE_EB_init, DM_EB, DM_Wald_EB, wilcoxon")

}





# now annotate the test results 
# Construct output filename from input filename
if (grepl("\\.txt$", out_resultfile)) {
  out_result_annotated <- sub("\\.txt$", ".annotated.txt", out_resultfile)
} else {
  out_result_annotated <- paste0(out_resultfile, ".annotated.txt")
}

# Read test results first to get the list of genes
message(paste("Reading test results from", out_resultfile, "..."))
tests <- read.table(out_resultfile, header = TRUE, stringsAsFactors = FALSE)
tests <- as_tibble(tests)
message(paste("Loaded", nrow(tests), "test results."))

# Get unique genes
unique_genes <- unique(tests$gene)
message(paste("Processing", length(unique_genes), "unique genes."))

# Function to safely read a specific gene file
read_gene_split <- function(gene_id) {
  # Construct filename based on gene ID
  f <- file.path(countdir, paste0(gene_id, ".", split, ".txt"))
  
  if (!file.exists(f)) {
    # Try resolving tilde expansion if needed (Sys.glob might help if path has ~)
    f_expanded <- Sys.glob(f)
    if (length(f_expanded) > 0) f <- f_expanded[1]
  }

  if (file.exists(f)) {
     tryCatch({
      df <- read.table(f, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
      # Ensure gene column is character to match tests
      df$gene <- as.character(df$gene)
      if("source" %in% names(df)) df$source <- as.character(df$source)
      if("sink" %in% names(df)) df$sink <- as.character(df$sink)
      return(as_tibble(df))
    }, error = function(e) {
      warning(paste("Failed to read", f, ":", e$message))
      return(NULL)
    })
  } else {
    return(NULL)
  }
}



# Iterate over unique genes and read their corresponding files
n_cores <- 32
message(paste("Using", n_cores, "cores for reading split files."))
split_list <- mclapply(unique_genes, read_gene_split, mc.cores = n_cores)
# Filter out NULLs
split_list <- split_list[!sapply(split_list, is.null)]
if (length(split_list) == 0) {
  stop("No matching split files found for any genes in the test file.")
}
splits <- bind_rows(split_list)
message(paste("Loaded split data for", length(unique(splits$gene)), "genes."))


# Combine data
# Using left_join to keep all tests, filling NAs for missing annotations
message("Merging data...")
merged_data <- left_join(tests, bipartitions, by = c("gene", "event"))
message(paste("Merged dataset has", nrow(merged_data), "rows."))

# Write output
write.table(merged_data, out_result_annotated, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Successfully wrote annotated tests to", out_result_annotated))
      
