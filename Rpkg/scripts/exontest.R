library(parallel)
library(MASS)
library(matrixStats)
library(glmmTMB)
library(tidyverse)
library(grase)
library(optparse)

# --- Main Execution Script ---

outdir = '~/GrASE_simulation/bipartition.test'
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
  make_option(c("-o", "--outdir"), type="character", 
              help="output directory that contains the test results", metavar="character"),
  make_option(c("-c", "--countdir"), type="character", 
              help="count dir that contains the count and split files by gene", metavar="character"),
  make_option(c("-s", "--splittype"), type="character", 
              help="split type: bipartition or multinomial or n_choose_2", metavar="character"),
  make_option(c("-p", "--phi"), type="character",
              help="filename for phi estimates", metavar="character"),
  make_option(c("--prec"), type="character", default=NULL,
              help="filename for precision estimates (multinomial models)", metavar="character"),
  make_option(c("-m", "--model"), type="character", 
              help="model name", metavar="character"),
  make_option(c("--cond1"), type="character",
              help="condition 1 name (Group 1 - Group 2)", metavar="character"),
  make_option(c("--cond2"), type="character",
              help="condition 2 name (Reference Group)", metavar="character"),
  make_option(c("--use_phi_loess"), action="store_true", default=FALSE,
              help="use loess phi trend (log(phi) ~ log(baseMean)) as EB shrinkage target instead of global mean"),
  make_option(c("--independent_filtering"), action="store_true", default=FALSE,
              help="use DESeq2-style independent filtering (filter on baseMean) before BH FDR correction"),
  make_option(c("--pseudocount"), type="integer", default=0L,
              help="pseudocount added to diff (ref reads) and n (total) before testing; enables detection of all-or-nothing switches when ref coverage is zero [default: 0]"),
  make_option(c("--prec_subsample_n"), type="integer", default=NULL,
              help="number of events to subsample for DM precision estimation (multinomial models); NULL uses all events [default: NULL]", metavar="integer"),
  make_option(c("--use_prec_loess"), action="store_true", default=FALSE,
              help="use loess prec trend (log(prec) ~ log(baseMean)) as EB shrinkage target instead of global mean (multinomial models)"),
  make_option(c("--padj_threshold"), type="double", default=0.05,
              help="padj threshold for the significant column [default: 0.05]", metavar="double"),
  make_option(c("--delta"), type="double", default=0,
              help="lfc_diff_net threshold for the significant column; events with lfc_diff_net <= delta are not significant [default: 0]", metavar="double")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)
if (!is.null(opt$file)) {
  masterfile = opt$file
}
if (!is.null(opt$outdir)) {
  outdir = opt$outdir
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
if (!is.null(opt$prec)) {
  prec_file = opt$prec
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
phi_trend        <- isTRUE(opt$use_phi_loess)
indep_filter     <- isTRUE(opt$independent_filtering)
pseudocount      <- as.integer(opt$pseudocount)
prec_subsample_n <- if (!is.null(opt$prec_subsample_n)) as.integer(opt$prec_subsample_n) else NULL
prec_trend       <- isTRUE(opt$use_prec_loess)
padj_thr         <- as.double(opt$padj_threshold)
delta            <- as.double(opt$delta)

# Expand tilde in paths
outdir <- path.expand(outdir)
countdir <- path.expand(countdir)

# Add a 'significant' column: padj < padj_thr AND (lfc_diff_net > delta OR lfc_diff_net absent/NA)
add_significant <- function(res, padj_thr, delta) {
  has_lfc <- "lfc_diff_net" %in% names(res)
  lfc_ok <- if (has_lfc) {
    # Events are significant if lfc_diff_net > delta.
    # If lfc_diff_net is NA, it's not filtered out.
    (res$lfc_diff_net > delta) | is.na(res$lfc_diff_net)
  } else {
    TRUE # No LFC info, don't filter
  }
  res$significant <- !is.na(res$padj) & res$padj < padj_thr & lfc_ok
  res
}
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
exoncnt_master <- paste0(countdir,'/', masterfile)
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

# baseMean: total per-event coverage across both diffs and ref
if (split != "multinomial") {
  baseMean_df <- splitcnts %>%
    mutate(n_total = rowSums(cbind(diff1, diff2, ref), na.rm = TRUE)) %>%
    group_by(gene, event) %>%
    summarise(baseMean = mean(n_total, na.rm = TRUE), .groups = "drop")

  make_comparison <- function(sc, diff_col, ref_col, comp_name) {
    sc %>%
      filter(!is.na(.data[[diff_col]]) & !is.na(.data[[ref_col]])) %>%
      mutate(diff = .data[[diff_col]],
             ref  = .data[[ref_col]],
             n    = .data[[diff_col]] + .data[[ref_col]],
             comparison = comp_name)
  }

  sc_d1r  <- make_comparison(splitcnts, "diff1", "ref",   "diff1_vs_ref")
  sc_d2r  <- make_comparison(splitcnts, "diff2", "ref",   "diff2_vs_ref")

  # pseudocount: same logic as before, applied to each comparison's ref where ref == 0
  if (pseudocount > 0L) {
    apply_pseudo <- function(sc) {
      zero_ref <- sc$ref == 0L
      sc$ref[zero_ref] <- pseudocount
      sc$n[zero_ref]   <- sc$n[zero_ref] + pseudocount
      sc
    }
    sc_d1r  <- apply_pseudo(sc_d1r)
    sc_d2r  <- apply_pseudo(sc_d2r)
  }

  lfc_d1r  <- compute_lfc_summary(sc_d1r)  %>% mutate(comparison = "diff1_vs_ref")
  lfc_d2r  <- compute_lfc_summary(sc_d2r)  %>% mutate(comparison = "diff2_vs_ref")
  lfc_summary_all <- bind_rows(lfc_d1r, lfc_d2r)
  comparisons_list <- list(sc_d1r, sc_d2r)
}
# Global Precision Estimation
if (model == 'glmmTMB_prior' || model == 'glmmTMB_fixedEB') {

  if (is.null(opt$phi)) stop("--phi is required for model '", model, "'")
  if (!is.null(opt$prec)) warning("--prec is not used for model '", model, "' and will be ignored")

  if (file.exists(paste0(outdir, "/", phifile))) {
    phi_df <- read.table(paste0(outdir,'/', phifile), header=TRUE, row.names=NULL)
  } else {
    # Estimate phi for each comparison separately
    estimate_phi_for_comparison <- function(sc, comp_name) {
      if (nrow(sc) == 0) return(NULL)
      grouped_data <- group_by_event(sc, 'diff', 'n')
      phi_progress_file <- paste0(outdir, "/phi_progress_", comp_name, ".log")
      error_log_file <- paste0("phi.glmmtmb.errors.", comp_name, ".log")

      phi_list <- mclapply(grouped_data, function(dd) {
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
      }, mc.cores = 32)

      phi_df_comp <- bind_rows(Filter(Negate(is.null), phi_list))
      if (nrow(phi_df_comp) > 0) {
        phi_df_comp$comparison <- comp_name
      }
      return(phi_df_comp)
    }

    message("Estimating phi for diff1 vs ref...")
    phi_d1r <- estimate_phi_for_comparison(sc_d1r, "diff1_vs_ref")
    message("Estimating phi for diff2 vs ref...")
    phi_d2r <- estimate_phi_for_comparison(sc_d2r, "diff2_vs_ref")

    phi_df <- bind_rows(phi_d1r, phi_d2r)
    rm(phi_d1r, phi_d2r)

    write.table(phi_df, file=paste0(outdir,'/', phifile), quote=FALSE, sep="\t", row.names = FALSE)
  }
} else if (model == 'multinomial_plugin_dm_EB') {
  # Multinomial DM Moderation Pipeline
  # group_by_event on counts instead of diff/n for multinomial
  if (is.null(opt$prec)) stop("--prec is required for model '", model, "'")
  if (!is.null(opt$phi)) warning("--phi is not used for model '", model, "' and will be ignored")

  # baseMean per event: mean total count per sample across all types
  baseMean_df_all <- splitcnts %>%
    group_by(gene, event, sample) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    group_by(gene, event) %>%
    summarise(baseMean = mean(total), .groups = "drop")

  if (file.exists(paste0(outdir,'/', prec_file))) {
    prec_raw <- read.table(paste0(outdir,'/', prec_file), header=TRUE, row.names=NULL)
    if (prec_trend) {
      prec_table <- moderate_prec_trend(prec_raw, baseMean_df_all)
    } else {
      prec_table <- moderate_prec_log_scale(prec_raw)
    }
  } else {
    grouped_counts <- splitcnts %>% group_by(gene, event) %>% group_split()
    n_events <- length(grouped_counts)

    # Subsample events for faster precision estimation if requested
    if (!is.null(prec_subsample_n) && prec_subsample_n < n_events) {
      message(sprintf("Subsampling %d of %d events for prec estimation.", prec_subsample_n, n_events))
      subsample_idx <- sample(n_events, prec_subsample_n)
      grouped_counts_prec <- grouped_counts[subsample_idx]
    } else {
      grouped_counts_prec <- grouped_counts
    }

    progress_file <- paste0(outdir, "/prec_progress.log")
    prec_list <- mclapply(grouped_counts_prec, function(dd) {
      result <- tryCatch(
        prec_estimate_plugin_dm(dd),
        error = function(e) {
          msg <- sprintf("[%s] ERROR gene=%s event=%s (PID %d): %s\n",
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                         dd$gene[1], dd$event[1], Sys.getpid(), conditionMessage(e))
          cat(msg, file = "multinomial.prec.errors.log", append = TRUE)
          NULL
        }
      )
      cat(sprintf("%s\t%s\n", dd$gene[1], dd$event[1]),
          file = progress_file, append = TRUE)
      result
    }, mc.cores=32)
    prec_raw <- bind_rows(Filter(Negate(is.null), prec_list))
    rm(prec_list)
    write.table(prec_raw, paste0(outdir,'/', prec_file), sep="\t", quote=F, row.names = FALSE)

    if (prec_trend) {
      prec_table <- moderate_prec_trend(prec_raw, baseMean_df_all)
    } else {
      prec_table <- moderate_prec_log_scale(prec_raw)
    }
  }
  prec_mod_file <- sub("\\.([^.]+)$", ".moderated.\\1", prec_file)
  write.table(prec_table, paste0(outdir, '/', prec_mod_file), sep="\t", quote=FALSE, row.names=FALSE)
  rm(baseMean_df_all)
}

# Helper: run one comparison through group_by_event + mclapply
run_one_comparison <- function(sc, test_fn, err_log, ...) {
  if (nrow(sc) == 0) return(NULL)
  gd <- group_by_event(sc, 'diff', 'n')
  res_list <- mclapply(gd, function(dd) {
    tryCatch(test_fn(dd, ...), error = function(e) {
        msg <- sprintf("[%s] ERROR gene=%s event=%s (PID %d): %s\n",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       dd$gene[1], dd$event[1], Sys.getpid(), conditionMessage(e))
        cat(msg, file = err_log, append = TRUE)
        NULL
    })
  }, mc.cores = 32)
  res <- bind_rows(Filter(Negate(is.null), res_list))
  if (nrow(res) > 0) res$comparison <- sc$comparison[1]
  res
}

# Per Gene Model Estimation
if (model == 'glmmTMB_prior') {
  # This model runs on each of the 2 comparisons (diff1_vs_ref, diff2_vs_ref).
  # Phi estimation (above) already ran on sc_d1r.
  phi_trimmed <- phi_df[!is.na(phi_df$phi) & phi_df$phi < 1e+10 & phi_df$phi > 0,'phi']
  if (phi_trend) {
    message("Using loess trend for glmmTMB_prior model's prior parameters.")
    baseMean_d1r  <- sc_d1r  %>% group_by(gene, event) %>% summarise(baseMean = mean(n, na.rm = TRUE), .groups = "drop") %>% mutate(comparison = "diff1_vs_ref")
    baseMean_d2r  <- sc_d2r  %>% group_by(gene, event) %>% summarise(baseMean = mean(n, na.rm = TRUE), .groups = "drop") %>% mutate(comparison = "diff2_vs_ref")
    baseMean_df_all <- bind_rows(baseMean_d1r, baseMean_d2r)

    # fit normal to log(phi)
    log_phi_vals <- log(phi_trimmed)
    median_logphi <- median(log_phi_vals, na.rm = TRUE)
    fit_logphi <- fitdistr(log_phi_vals, "normal")
    mean_logphi <- fit_logphi$estimate["mean"]
    sd_logphi   <- fit_logphi$estimate["sd"]
    prior_disp <- data.frame(prior = sprintf("normal(%g,%g)", mean_logphi, sd_logphi),
      class = "fixef_disp",
      coef = "", # Explicitly target the intercept
      stringsAsFactors = FALSE
    )
    print(prior_disp)
      # Get trended phi values. Assumes moderate_phi_trend returns z_trend.
      phi_table_trend <- moderate_phi_trend(phi_df, baseMean_df_all)

    result_list_all <- lapply(comparisons_list, run_one_comparison,
                            test_fn = test_model_glmmTMB_with_prior,
                            err_log = "glmmtmb_withprior.errors.log",
                            prior_disp, z_bar = median_logphi)
    # Calculate residual SD from trend, per comparison
    phi_df_with_trend <- phi_df %>%
      left_join(phi_table_trend, by = c("gene", "event", "comparison")) %>%
      filter(!is.na(phi) & phi > 0 & !is.na(z_trend))

    prior_params <- phi_df_with_trend %>%
      mutate(resid = log(phi) - z_trend) %>%
      group_by(comparison) %>%
      summarise(sd_resid = sd(resid, na.rm = TRUE), .groups = "drop")

    # Global fallback for SD
    global_sd_resid <- sd(phi_df_with_trend$resid, na.rm = TRUE)
    if (is.na(global_sd_resid)) {
        phi_trimmed <- phi_df[!is.na(phi_df$phi) & phi_df$phi < 1e+10 & phi_df$phi > 0, 'phi']
        global_sd_resid <- sd(log(phi_trimmed), na.rm = TRUE)
    }
    prior_params <- prior_params %>%
        mutate(sd_resid = ifelse(is.na(sd_resid), global_sd_resid, sd_resid))

    # Global fallbacks for mean and median
    median_logphi <- median(log(phi_df$phi), na.rm = TRUE)
    global_mean_logphi <- mean(log(phi_df$phi), na.rm = TRUE)

    # Join prior parameters to each comparison dataset
    join_prior_params <- function(sc, trend_table, params_table) {
      sc %>%
        left_join(trend_table %>% select(gene, event, comparison, z_trend), by = c("gene", "event", "comparison")) %>%
        left_join(params_table, by = "comparison")
    }
    sc_d1r_prior  <- join_prior_params(sc_d1r, phi_table_trend, prior_params)
    sc_d2r_prior  <- join_prior_params(sc_d2r, phi_table_trend, prior_params)

    # Impute missing prior parameters
    impute_prior_params <- function(df, comp_name) {
      if (nrow(df) == 0) return(df)
      fallback_sd <- prior_params$sd_resid[prior_params$comparison == comp_name]
      if (length(fallback_sd) == 0 || is.na(fallback_sd)) fallback_sd <- global_sd_resid

      df$z_trend[is.na(df$z_trend)] <- global_mean_logphi
      df$sd_resid[is.na(df$sd_resid)] <- fallback_sd
      df
    }
    sc_d1r_prior  <- impute_prior_params(sc_d1r_prior, "diff1_vs_ref")
    sc_d2r_prior  <- impute_prior_params(sc_d2r_prior, "diff2_vs_ref")

    comparisons_list_prior <- list(sc_d1r_prior, sc_d2r_prior)

    # Wrapper to construct prior inside mclapply
    test_fn_wrapper_trend <- function(dd) {
      mean_logphi <- dd$z_trend[1]
      sd_logphi   <- dd$sd_resid[1]
      prior_disp <- data.frame(prior = sprintf("normal(%g,%g)", mean_logphi, sd_logphi),
                               class = "fixef_disp", coef = "",
                               stringsAsFactors = FALSE)
      test_model_glmmTMB_with_prior(dd, prior_disp, z_bar = median_logphi)
    }

    result_list_all <- lapply(comparisons_list_prior, run_one_comparison,
                              test_fn = test_fn_wrapper_trend,
                              err_log = "glmmtmb_withprior_trend.errors.log")
  } else {
    # Original behavior: global prior for all events
    phi_trimmed <- phi_df[!is.na(phi_df$phi) & phi_df$phi < 1e+10 & phi_df$phi > 0,'phi']
    log_phi_vals <- log(phi_trimmed)
    median_logphi <- median(log_phi_vals, na.rm = TRUE)
    fit_logphi <- fitdistr(log_phi_vals, "normal")
    mean_logphi <- fit_logphi$estimate["mean"]
    sd_logphi   <- fit_logphi$estimate["sd"]
    prior_disp <- data.frame(prior = sprintf("normal(%g,%g)", mean_logphi, sd_logphi),
      class = "fixef_disp",
      coef = "", # Explicitly target the intercept
      stringsAsFactors = FALSE
    )
    print(prior_disp)

    result_list_all <- lapply(comparisons_list, run_one_comparison,
                            test_fn = test_model_glmmTMB_with_prior,
                            err_log = "glmmtmb_withprior.errors.log",
                            prior_disp, z_bar = median_logphi)
  }

  results <- bind_rows(result_list_all)
  results <- results %>%
    left_join(baseMean_df, by = c("gene", "event")) %>%
    left_join(lfc_summary_all, by = c("gene", "event", "comparison"))

  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=indep_filter, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }

  results <- add_significant(results, padj_thr, delta)
  write.table(results, file=out_resultfile, quote=FALSE, sep="\t", row.names = FALSE)

} else if (model == 'glmmTMB_fixedEB') {

  ## Empirical Bayes hyperparameters: prior mean and variance
  if (phi_trend) {
    ## Trend-based moderation: loess(log(phi) ~ log(baseMean)) as shrinkage target.
    ## baseMean is now comparison-specific.
    baseMean_d1r  <- sc_d1r  %>% group_by(gene, event) %>% summarise(baseMean = mean(n, na.rm = TRUE), .groups = "drop") %>% mutate(comparison = "diff1_vs_ref")
    baseMean_d2r  <- sc_d2r  %>% group_by(gene, event) %>% summarise(baseMean = mean(n, na.rm = TRUE), .groups = "drop") %>% mutate(comparison = "diff2_vs_ref")
    baseMean_df_all <- bind_rows(baseMean_d1r, baseMean_d2r)
    phi_table <- moderate_phi_trend(phi_df, baseMean_df_all)
    fallback_z <- median(phi_table$z_mod, na.rm = TRUE)
    # Join moderated phi to each comparison dataset by gene, event, AND comparison
    sc_d1r_eb  <- sc_d1r  %>% left_join(phi_table, by = c("gene", "event", "comparison"))
    sc_d2r_eb  <- sc_d2r  %>% left_join(phi_table, by = c("gene", "event", "comparison"))
    sc_d1r_eb$z_mod[is.na(sc_d1r_eb$z_mod)] <- fallback_z
    sc_d2r_eb$z_mod[is.na(sc_d2r_eb$z_mod)] <- fallback_z
  } else {
    ## Global-mean moderation (original behaviour)
    # To make moderation more robust, apply it independently to each comparison type,
    # acknowledging that their phi distributions may differ.
    phi_table <- phi_df %>%
      filter(!is.na(comparison)) %>%
      group_by(comparison) %>%
      group_modify(~ moderate_phi_log_scale(.x)) %>%
      ungroup()

    # Join moderated phi values.
    sc_d1r_eb  <- sc_d1r  %>% left_join(phi_table, by = c("gene", "event", "comparison"))
    sc_d2r_eb  <- sc_d2r  %>% left_join(phi_table, by = c("gene", "event", "comparison"))

    # Impute missing z_mod with the median of their respective comparison group.
    impute_z_mod <- function(df, p_table) {
        if (nrow(df) == 0) return(df)
        comp_name <- df$comparison[1]
        fallback_z <- median(p_table$z_mod[p_table$comparison == comp_name], na.rm = TRUE)
        if (is.na(fallback_z)) fallback_z <- median(p_table$z_mod, na.rm = TRUE) # Global fallback
        if (is.na(fallback_z)) fallback_z <- 0 # Ultimate fallback log(1)
        df$z_mod[is.na(df$z_mod)] <- fallback_z
        df
    }
    sc_d1r_eb  <- impute_z_mod(sc_d1r_eb, phi_table)
    sc_d2r_eb  <- impute_z_mod(sc_d2r_eb, phi_table)
  }
  comparisons_list_eb <- list(sc_d1r_eb, sc_d2r_eb)

  phi_mod_file <- sub("\\.([^.]+)$", ".moderated.\\1", phifile)
  write.table(phi_table, file=paste0(outdir, '/', phi_mod_file), sep="\t", quote=FALSE, row.names=FALSE)

  rm(splitcnts, phi_df, phi_table)
  if (exists("sc_d1r")) rm(sc_d1r, sc_d2r)

  result_list_all <- lapply(comparisons_list_eb, run_one_comparison,
                          test_fn = test_model_glmmTMB_EB,
                          err_log = "glmmtmb_EB.errors.log")

  results <- bind_rows(result_list_all)
  results <- results %>%
    left_join(baseMean_df, by = c("gene", "event")) %>%
    left_join(lfc_summary_all, by = c("gene", "event", "comparison"))

  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=indep_filter, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }

  results <- add_significant(results, padj_thr, delta)
  write.table(results, file=out_resultfile, quote=FALSE, sep="\t", row.names=FALSE)

} else if (model == 'multinomial_plugin_dm_EB') {

  splitcnts_eb <- splitcnts %>%
    left_join(prec_table, by = c("gene", "event"))
  # Fill NAs for events not covered by prec estimation (subsampled run, no prec_trend)
  if (any(is.na(splitcnts_eb$log_prec_mod))) {
    fallback_log  <- median(prec_table$log_prec_mod, na.rm = TRUE)
    fallback_prec <- exp(fallback_log)
    fallback_rho  <- pmax(pmin(1 / (1 + fallback_prec), 1 - 1e-6), 1e-6)
    na_rows <- is.na(splitcnts_eb$log_prec_mod)
    splitcnts_eb$log_prec_mod[na_rows] <- fallback_log
    splitcnts_eb$prec_mod[na_rows]     <- fallback_prec
    splitcnts_eb$rho_mod[na_rows]      <- fallback_rho
  }
  grouped_eb <- splitcnts_eb %>% group_by(gene, event) %>% group_split()

  test_fn <- test_model_multinomial_plugin_dm_EB
  err_log <- paste0(model, ".errors.log")

  res_dm <- mclapply(grouped_eb, function(dd) {
    tryCatch(
      test_fn(dd),
      error = function(e) {
        msg <- sprintf("[%s] ERROR gene=%s event=%s (PID %d): %s\n",
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       dd$gene[1], dd$event[1], Sys.getpid(), conditionMessage(e))
        cat(msg, file = err_log, append = TRUE)
        NULL
      }
    )
  }, mc.cores = 32)
  #stopCluster(cl)
  results <- bind_rows(res_dm)
  # lfc_summary is not computed for multinomial (no diff/ref columns)

  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=FALSE, theta=NULL, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }

  results <- add_significant(results, padj_thr, delta)
  write.table(results, file=out_resultfile, sep="\t", quote=F, row.names = FALSE)

} else if (model == 'wilcoxon') {
  if (!is.null(opt$phi)) warning("--phi is not used for model '", model, "' and will be ignored")

  result_list_all <- lapply(comparisons_list, run_one_comparison,
                          test_fn = test_model_wilcoxon,
                          err_log = "wilcoxon.errors.log")

  results <- bind_rows(result_list_all)
  results <- results %>%
    left_join(baseMean_df, by = c("gene", "event")) %>%
    left_join(lfc_summary_all, by = c("gene", "event", "comparison"))

  if (exists("pvalueAdjustment") && nrow(results) > 0) {
      results$pvalue <- results$p.value
      results <- pvalueAdjustment(results, independentFiltering=indep_filter, alpha=0.05, pAdjustMethod="BH")
      results$pvalue <- NULL
  }

  results <- add_significant(results, padj_thr, delta)
  write.table(results, file=out_resultfile, quote=FALSE, sep="\t", row.names = FALSE)

} else if (model == 'glmmTMB_no_prior') {
  if (!is.null(opt$phi)) warning("--phi is not used for model '", model, "' and will be ignored")

  result_list_all <- lapply(comparisons_list, run_one_comparison,
                          test_fn = function(dd) suppressMessages(test_model_glmmTMB_without_prior(dd)),
                          err_log = "glmmtmb_noprior.errors.log")

  results <- bind_rows(result_list_all)
  results <- results %>%
    left_join(baseMean_df, by = c("gene", "event")) %>%
    left_join(lfc_summary_all, by = c("gene", "event", "comparison"))

  if (exists("pvalueAdjustment") && nrow(results) > 0) {
    results$pvalue <- results$p.value
    results <- pvalueAdjustment(results, independentFiltering=indep_filter, alpha=0.05, pAdjustMethod="BH")
    results$pvalue <- NULL
  }

  results <- add_significant(results, padj_thr, delta)
  write.table(results, file=out_resultfile, quote=FALSE, sep="\t", row.names = FALSE)

} else {
  print ("model misspecified")
  print ("choose between glmmTMB_prior, glmmTMB_no_prior, glmmTMB_fixedEB, multinomial_plugin_dm_EB, wilcoxon")

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
merged_data <- left_join(tests, splits, by = c("gene", "event"))
message(paste("Merged dataset has", nrow(merged_data), "rows."))

# Write output
write.table(merged_data, out_result_annotated, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Successfully wrote annotated tests to", out_result_annotated))


# Each bipartition event was tested twice (_s1 and _s2 sides). 
if (split == 'bipartition' || split == 'n_choose_2') {
  if (nrow(merged_data) > 0 && 'setdiff1' %in% names(merged_data)) {

    # Min p-value combination for bipartition_both:
    # Take the minimum p-value across sides per event, union the setdiff exonic parts,
    # and re-apply BH FDR correction on the reduced hypothesis set (N not 2N).
    # Per-event: min p-value + union of setdiff exon parts
    combo_df <- merged_data %>%
      group_by(gene, event) %>%
      summarise(
        p_min         = {
          pv <- p.value[!is.na(p.value) & p.value > 0 & p.value <= 1]
          if (length(pv) == 0) NA_real_ else min(pv)
        },
        setdiff_union = paste(unique(trimws(unlist(strsplit(
                                na.omit(c(setdiff1, setdiff2)), ",")))), collapse = ","),
        .groups = "drop"
      )

    # Primary row per event: prefer the side with the lower p-value
    primary_df <- merged_data %>%
      group_by(gene, event) %>%
      slice(which.min(p.value)) %>%
      ungroup()

    min_data <- primary_df %>%
      left_join(combo_df, by = c("gene", "event")) %>%
      mutate(
        p.value = p_min,
        setdiff = setdiff_union
      ) %>%
      select(-p_min, -setdiff_union)

    if (exists("pvalueAdjustment") && nrow(min_data) > 0) {
      min_data$pvalue <- min_data$p.value
      min_data <- pvalueAdjustment(min_data, independentFiltering=indep_filter, alpha=0.05, pAdjustMethod="BH")
      min_data$pvalue <- NULL
    } else {
      min_data$padj <- p.adjust(min_data$p.value, method = "BH")
    }

    min_data <- add_significant(min_data, padj_thr, delta)
    out_mincomb <- sub("\\.annotated\\.txt$", ".mincomb.annotated.txt",
                      out_result_annotated)
    write.table(min_data, out_mincomb, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste("Min-p combined results written to", out_mincomb))

    # Fisher p-value combination for bipartition_both:
    # Combine p-values from _s1 and _s2 using Fisher's method:
    #   X = -2 * sum(log(p))  ~  chi-squared with 2*k df  (k = number of sides)
    # When only one side has a valid p-value, use that p directly (1 df = 2).
    fisher_combo_df <- merged_data %>%
      group_by(gene, event) %>%
      summarise(
        p_fisher      = {
          pv <- p.value[!is.na(p.value) & p.value > 0 & p.value <= 1]
          if (length(pv) == 0) NA_real_
          else pchisq(-2 * sum(log(pv)), df = 2 * length(pv), lower.tail = FALSE)
        },
        setdiff_union = paste(unique(trimws(unlist(strsplit(
                                na.omit(c(setdiff1, setdiff2)), ",")))), collapse = ","),
        .groups = "drop"
      )

    fisher_data <- primary_df %>%
      left_join(fisher_combo_df, by = c("gene", "event")) %>%
      mutate(
        p.value = p_fisher,
        setdiff = setdiff_union
      ) %>%
      select(-p_fisher, -setdiff_union)

    if (exists("pvalueAdjustment") && nrow(fisher_data) > 0) {
      fisher_data$pvalue <- fisher_data$p.value
      fisher_data <- pvalueAdjustment(fisher_data, independentFiltering=indep_filter, alpha=0.05, pAdjustMethod="BH")
      fisher_data$pvalue <- NULL
    } else {
      fisher_data$padj <- p.adjust(fisher_data$p.value, method = "BH")
    }

    fisher_data <- add_significant(fisher_data, padj_thr, delta)
    out_fisher <- sub("\\.annotated\\.txt$", ".fisher_combined.annotated.txt",
                      out_result_annotated)
    write.table(fisher_data, out_fisher, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste("Fisher combined results written to", out_fisher))
  }
}
      
