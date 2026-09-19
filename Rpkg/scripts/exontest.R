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
model = 'betabinom_EBmap'
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
  make_option(c("--contrasts"), type="character", default=NULL,
              help="comma-separated contrasts in trt:ref format, e.g. B:A,C:A (overrides --cond1/--cond2 for multi-group runs)", metavar="character"),
  make_option(c("--use_phi_loess"), action="store_true", default=FALSE,
              help="use loess phi trend (log(phi) ~ log(baseMean)) as EB shrinkage target instead of global median"),
  make_option(c("--no_independent_filtering"), action="store_false", default=TRUE,
              dest="independent_filtering",
              help="disable DESeq2-style independent filtering (filter on baseMean) before FDR correction (on by default)"),
  make_option(c("--padj_method"), type="character", default="nested_BH",
              help="p-value adjustment method: nested_BH (default) or any p.adjust method e.g. BH", metavar="character"),
  make_option(c("--pseudocount"), type="integer", default=0L,
              help="pseudocount added to diff (ref reads) and n (total) before testing; enables detection of all-or-nothing switches when ref coverage is zero [default: 0]"),
make_option(c("--use_prec_loess"), action="store_true", default=FALSE,
              help="use loess prec trend (log(prec) ~ log(baseMean)) as EB shrinkage target instead of global mean (multinomial models)"),
  make_option(c("--padj_threshold"), type="double", default=0.01,
              help=paste("padj threshold for the significant column AND the nested_BH gene screen;",
                         "both stages use this value [default: 0.01]"), metavar="double"),
  make_option(c("--delta"), type="double", default=0,
              help="lfc_diff_net threshold for the significant column; events with lfc_diff_net <= delta are not significant [default: 0]", metavar="double"),
  make_option(c("--min_reads"), type="double", default=10,
              help=paste("minimum read support for the tested distinct set, required in AT",
                         "LEAST ONE of the two contrasted groups. Either, not both: requiring",
                         "both would exclude on/off switches, where a path is absent in one",
                         "condition by design [default: 10]"), metavar="double"),
  make_option(c("--min_dpi_sj"), type="double", default=0.05,
              help=paste("as --min_dpi but for sides whose distinct set came from split reads",
                         "(merged exon+SJ runs). Split-read counts sit on a lower pi scale than",
                         "the exonic reference, so one absolute threshold is not scale-fair",
                         "[default: 0.05]"), metavar="double"),
  make_option(c("--min_dpi"), type="double", default=0.1,
              help="minimum |delta_pi| (path proportion effect size) for the significant column; events with |delta_pi| < min_dpi are not significant [default: 0.1]", metavar="double"),
  make_option(c("--mc_cores"), type="integer", default=32L,
              help="number of parallel worker processes for mclapply [default: %default]", metavar="integer")
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
if (!is.null(opt$contrasts)) {
  contrast_specs <- trimws(strsplit(opt$contrasts, ",")[[1]])
  ## Two contrast forms:
  ##   trt:ref      pairwise, 1 df -- unchanged, still the default
  ##   A+B+C        omnibus over K groups, K-1 df
  ## The omnibus form exists because the model is already cbind(y, n-y) ~ 0 +
  ## groups over every group in the run, so a K-1 df Wald test needs no refit.
  ## Its effect-size gates stay PAIRWISE (see is_significant): a range over K
  ## groups is an extreme-order statistic and drifts badly with K -- measured on
  ## DICE, median |dpi| range grows 7x from K=2 to K=13 and the share clearing
  ## 0.1 goes from 7.5% to 32.5%, so a fixed floor silently loosens as groups
  ## are added.
  contrasts <- lapply(contrast_specs, function(s) {
    if (grepl("+", s, fixed = TRUE)) {
      grps <- trimws(strsplit(s, "+", fixed = TRUE)[[1]])
      grps <- grps[nzchar(grps)]
      if (length(grps) < 2) stop("omnibus contrast needs >= 2 groups, got: ", s)
      if (anyDuplicated(grps)) stop("duplicate group in omnibus contrast: ", s)
      list(type = "omnibus", groups = grps,
           name = paste0("omnibus_", paste(grps, collapse = ".")))
    } else {
      parts <- trimws(strsplit(s, ":")[[1]])
      if (length(parts) != 2) stop("each pairwise contrast must be trt:ref, got: ", s)
      list(type = "pair", trt = parts[1], ref = parts[2],
           groups = c(parts[1], parts[2]),
           name = paste0(parts[1], "_vs_", parts[2]))
    }
  })
} else {
  contrasts <- list(list(type = "pair", trt = cond1, ref = cond2,
                        groups = c(cond1, cond2),
                        name = paste0(cond1, "_vs_", cond2)))
}
all_groups <- unique(unlist(lapply(contrasts, `[[`, "groups")))

## Omnibus contrasts need the contrast-matrix path, which only the glmmTMB
## beta-binomial models take. dirmult_EBplugin and wilcoxon instead SUBSET the
## data to the contrast's two groups, so an omnibus there would filter to an
## empty frame and yield nothing -- fail loudly instead.
if (any(vapply(contrasts, function(ct) identical(ct$type, "omnibus"), logical(1))) &&
    !model %in% c("betabinom_EBmap", "betabinom_EBapprox", "betabinom_MLE")) {
  stop("omnibus contrasts (A+B+C) are supported only by betabinom_EBmap, ",
       "betabinom_EBapprox and betabinom_MLE; model '", model,
       "' tests each contrast by subsetting to two groups.")
}
phi_trend        <- isTRUE(opt$use_phi_loess)
indep_filter     <- isTRUE(opt$independent_filtering)
pseudocount      <- as.integer(opt$pseudocount)
mc_cores         <- as.integer(opt$mc_cores)
prec_trend       <- isTRUE(opt$use_prec_loess)
padj_thr         <- as.double(opt$padj_threshold)
delta            <- as.double(opt$delta)
padj_method      <- as.character(opt$padj_method)
min_dpi          <- as.double(opt$min_dpi)
min_dpi_sj       <- as.double(opt$min_dpi_sj)
min_reads        <- as.double(opt$min_reads)

# Expand tilde in paths
outdir <- path.expand(outdir)
countdir <- path.expand(countdir)

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
exoncnt_master <- paste0(countdir,'/', masterfile)
out_prefix <- sub("^([^.]+\\.[^.]+).*", "\\1", masterfile)
out_resultfile=paste0(outdir,'/test_', out_prefix,'_', model, '.txt')


if (file.exists(exoncnt_master)) {
  splitcnts <- read.table(exoncnt_master, header=TRUE, row.names=NULL)
  # Include all groups mentioned in any contrast; drop samples from other groups.
  splitcnts$groups <- factor(splitcnts$groups, levels = all_groups)
  splitcnts <- splitcnts[!is.na(splitcnts$groups), ]
  ## Per-(gene,event,side,group) mean of the tested distinct count. Used below to
  ## build d_support = max over the two groups of each contrast, which gates the
  ## significance call (see --min_reads). Computed here because splitcnts is
  ## dropped later, and kept small: one row per gene/event/side/group.
  d_grp <- bind_rows(
    splitcnts %>% group_by(gene, event, groups) %>%
      summarise(m = mean(diff1, na.rm = TRUE), .groups = "drop") %>% mutate(side = "diff1"),
    splitcnts %>% group_by(gene, event, groups) %>%
      summarise(m = mean(diff2, na.rm = TRUE), .groups = "drop") %>% mutate(side = "diff2"))
  d_grp$gene <- as.character(d_grp$gene); d_grp$event <- as.character(d_grp$event)
  message(sprintf("Support table: %d (gene,event,side,group) means", nrow(d_grp)))

  ## L: a NAMED LIST of contrast matrices, columns = "groups<level>".
  ## A pairwise contrast is one row (+1 trt, -1 ref) and reproduces the former
  ## 1-df (est/se)^2 exactly. An omnibus over K groups is K-1 rows of successive
  ## differences, giving a full-rank K-1 df Wald test; wald_contrast() takes the
  ## numerical rank, so a gene missing one of the groups simply loses a df
  ## instead of failing.
  cols <- paste0("groups", all_groups)
  L <- lapply(contrasts, function(ct) {
    if (identical(ct$type, "omnibus")) {
      g <- paste0("groups", ct$groups)
      C <- matrix(0L, nrow = length(g) - 1L, ncol = length(cols),
                  dimnames = list(NULL, cols))
      for (i in seq_len(length(g) - 1L)) { C[i, g[i]] <- 1L; C[i, g[i + 1L]] <- -1L }
      C
    } else {
      C <- matrix(0L, nrow = 1L, ncol = length(cols), dimnames = list(NULL, cols))
      C[1, paste0("groups", ct$trt)] <-  1L
      C[1, paste0("groups", ct$ref)] <- -1L
      C
    }
  })
  names(L) <- vapply(contrasts, `[[`, character(1), "name")

  ## contrast name -> its group set, for add_support(). Previously the groups
  ## were recovered by splitting the contrast NAME on "_vs_", which breaks
  ## silently for a group whose own name contains "_vs_" (d_support becomes NA
  ## and the read floor never fires) and cannot express an omnibus at all.
  ctr_groups <- setNames(lapply(contrasts, `[[`, "groups"),
                         vapply(contrasts, `[[`, character(1), "name"))

  ## omnibus contrast name -> the names of the PAIRWISE contrasts it spans.
  ## is_significant() gates an omnibus row on whether any of these clears the
  ## effect-size floors, because the omnibus row itself has no single direction
  ## (delta_pi and lfc_diff_net are NA there). Only pairs that were actually
  ## requested can be used; an omnibus whose pairs were not also requested maps
  ## to an empty vector, and its rows then fail closed rather than passing on
  ## padj alone.
  pair_names <- vapply(Filter(function(ct) identical(ct$type, "pair"), contrasts),
                       `[[`, character(1), "name")
  omnibus_pairs <- lapply(Filter(function(ct) identical(ct$type, "omnibus"), contrasts),
    function(ct) {
      cand <- c(outer(ct$groups, ct$groups, function(a, b) paste0(a, "_vs_", b)))
      intersect(cand, pair_names)
    })
  names(omnibus_pairs) <- vapply(
    Filter(function(ct) identical(ct$type, "omnibus"), contrasts),
    `[[`, character(1), "name")
  if (!length(omnibus_pairs)) omnibus_pairs <- NULL
} else {
  print('exoncnt dataset file does not exist.')
  print('please combine the exoncounts into one file and specify the filename.')
  stop()
}

## ---- read support for the tested distinct set -------------------------------
## d_support = the LARGER of the two contrasted groups' mean distinct counts.
## max(), not min(): requiring both groups to clear the floor throws away exactly
## the on/off events the test exists to find (measured on the simulation GT:
## min() cost 200 true positives, max() cost 31, at the same precision).
##
## Kept separate from independent filtering, which cannot do this job: its
## statistic is baseMean = mean(diff_i + ref), the test's DENOMINATOR, and the
## reference usually dominates. A side with ref = 200 and diff = 0.3 has
## baseMean ~ 200 and sails through while carrying no information about D.
add_support <- function(res, contrast = NULL) {
  if (!exists("d_grp") || is.null(d_grp) || is.null(res) || !nrow(res)) return(res)
  if (!all(c("gene", "event", "comparison") %in% names(res))) return(res)
  ctr <- if (!is.null(contrast)) rep_len(as.character(contrast), nrow(res))
         else if ("contrast" %in% names(res)) as.character(res$contrast)
         else return(res)
  side <- sub("_vs_ref$", "", as.character(res$comparison))
  key  <- paste(d_grp$gene, d_grp$event, d_grp$side, d_grp$groups, sep = "\r")

  ## Groups come from ctr_groups (built with the contrasts), NOT from splitting
  ## the contrast name. Name-splitting silently yielded NA -- and so skipped the
  ## read floor entirely -- for any group containing "_vs_", and had no way to
  ## express a K-group omnibus. Unknown contrast names still give NA support,
  ## which leaves the row untested by the floor rather than wrongly filtered.
  grp_sets <- if (exists("ctr_groups")) ctr_groups[ctr] else vector("list", length(ctr))

  ## max over the contrast's groups: support in ANY ONE of them is enough.
  ## min() would discard all-or-nothing switches, which is what the test is for
  ## (measured: min cost 200 true positives vs 31 for max, at equal precision).
  ## Loop over CONTRASTS (a handful), not rows. A per-row vapply does one
  ## match() per row against a key vector with ~1.9M entries on the TSS/TTS
  ## arm -- order 10^12 operations, which ran for hours single-threaded with no
  ## output. Each contrast now contributes one vectorised match() per group,
  ## which is what the original pairwise code did (two matches, trt and ref).
  res$d_support <- NA_real_
  for (cn in unique(ctr)) {
    gs <- ctr_groups[[cn]]
    if (is.null(gs) || !length(gs)) next
    rows <- which(ctr == cn)
    acc  <- rep(NA_real_, length(rows))
    for (g in gs) {
      m <- d_grp$m[match(paste(res$gene[rows], res$event[rows], side[rows], g,
                               sep = "\r"), key)]
      acc <- pmax(acc, m, na.rm = TRUE)
    }
    res$d_support[rows] <- acc
  }
  res
}

## Drop unsupported sides from the test universe BEFORE p-value adjustment, by
## setting p.value to NA. The row is KEPT (padj then comes back NA), matching the
## independent-filtering convention: downstream evaluation must count it as a
## non-call still present in the universe. Deleting the rows instead would remove
## ground-truth positives from the denominator and inflate recall.
filter_unsupported <- function(res, contrast = NULL) {
  if (min_reads <= 0 || is.null(res) || !nrow(res)) return(res)
  res <- add_support(res, contrast)
  if (!"d_support" %in% names(res)) return(res)
  bad <- !is.na(res$d_support) & res$d_support < min_reads & !is.na(res$p.value)
  if (any(bad)) {
    res$p.value[bad] <- NA_real_
    message(sprintf("Support filter: %d of %d tests below %g reads in both contrasted groups -> not tested",
                    sum(bad), nrow(res), min_reads))
  }
  res
}

if (split != "multinomial") {

  make_comparison <- function(sc, diff_col, ref_col, comp_name) {
    sc %>%
      filter(!is.na(.data[[diff_col]]) & !is.na(.data[[ref_col]])) %>%
      mutate(diff = .data[[diff_col]],
             ref  = .data[[ref_col]],
             n    = .data[[diff_col]] + .data[[ref_col]],
             comparison = comp_name) %>%
      select(gene, event, groups, diff, ref, n, comparison)
  }

  sc_d1r  <- make_comparison(splitcnts, "diff1", "ref",   "diff1_vs_ref")
  sc_d2r  <- make_comparison(splitcnts, "diff2", "ref",   "diff2_vs_ref")

  # baseMean: per-comparison coverage n = diff_i + ref (raw, before pseudocount),
  # keyed by (gene, event, comparison). This is each beta-binomial test's denominator
  # and the correct independent-filtering statistic; a two-diff event's diff1_vs_ref
  # and diff2_vs_ref rows therefore get distinct baseMean values.
  baseMean_df <- bind_rows(
    sc_d1r %>% group_by(gene, event, comparison) %>%
      summarise(baseMean = mean(n, na.rm = TRUE), .groups = "drop"),
    sc_d2r %>% group_by(gene, event, comparison) %>%
      summarise(baseMean = mean(n, na.rm = TRUE), .groups = "drop")
  )

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

  ## PAIRWISE contrasts only. pi_trt, pi_ref, delta_pi, lfc_diff and lfc_ref are
  ## all two-group quantities, so an omnibus contrast has none of them -- calling
  ## compute_lfc_summary() with ct$trt = NULL builds the column name
  ## "mean_diff_" and aborts. Omnibus rows therefore join to NA on these
  ## columns, which is what is_significant() expects: it gates them on their
  ## constituent pairwise contrasts instead (omnibus_pairs).
  lfc_summary_all <- bind_rows(lapply(
    Filter(function(ct) identical(ct$type, "pair"), contrasts), function(ct) {
      ctr_name <- ct$name
      bind_rows(
        compute_lfc_summary(sc_d1r, cond_ref = ct$ref, cond_trt = ct$trt) %>%
          mutate(comparison = "diff1_vs_ref", contrast = ctr_name),
        compute_lfc_summary(sc_d2r, cond_ref = ct$ref, cond_trt = ct$trt) %>%
          mutate(comparison = "diff2_vs_ref", contrast = ctr_name)
      )
  }))
  comparisons_list <- list(sc_d1r, sc_d2r)
}
# Global Precision Estimation
if (model == 'betabinom_EBmap' || model == 'betabinom_EBapprox') {

  if (is.null(opt$phi)) stop("--phi is required for model '", model, "'")
  if (!is.null(opt$prec)) warning("--prec is not used for model '", model, "' and will be ignored")
  rm(splitcnts)
  gc()

  if (file.exists(paste0(outdir, "/", phifile))) {
    phi_df <- read.table(paste0(outdir,'/', phifile), header=TRUE, row.names=NULL)
  } else {
    message("Estimating phi for diff1 vs ref...")
    phi_d1r <- estimate_phi_for_comparison(sc_d1r, "diff1_vs_ref", outdir = outdir, mc_cores = mc_cores)
    message("Estimating phi for diff2 vs ref...")
    phi_d2r <- estimate_phi_for_comparison(sc_d2r, "diff2_vs_ref", outdir = outdir, mc_cores = mc_cores)

    phi_df <- bind_rows(phi_d1r, phi_d2r)
    rm(phi_d1r, phi_d2r)

    write.table(phi_df, file=paste0(outdir,'/', phifile), quote=FALSE, sep="\t", row.names = FALSE)
  }
} else if (model == 'dirmult_EBplugin') {
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

    progress_file <- paste0(outdir, "/prec_progress.log")
    prec_list <- mclapply(grouped_counts, function(dd) {
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
    }, mc.cores = mc_cores)
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

# Per Gene Model Estimation
if (model == 'betabinom_EBmap') {
  # Prior parameters from MLE phi estimates.
  # Restrict to events with an identifiable beta-binomial dispersion when that
  # flag is available: boundary/non-identifiable events (phi at the binomial
  # limit) form a spurious second mode that drags the fitted mean up and
  # inflates the sd, weakening the prior.  They are still tested (moderated),
  # just excluded from defining the prior.
  phi_prior_df <- phi_df[!is.na(phi_df$phi) & phi_df$phi < 1e+10 & phi_df$phi > 0, ]
  if ("identifiable" %in% names(phi_prior_df)) {
    phi_ident <- phi_prior_df[phi_prior_df$identifiable %in% TRUE, ]
    if (nrow(phi_ident) >= 10) {
      phi_prior_df <- phi_ident
    } else {
      warning("Fewer than 10 identifiable phi estimates; using all events for the prior.")
    }
  }
  phi_trimmed  <- phi_prior_df$phi
  log_phi_vals <- log(phi_trimmed)
  fit_logphi   <- fitdistr(log_phi_vals, "normal")
  mean_logphi  <- fit_logphi$estimate["mean"]
  sd_logphi    <- fit_logphi$estimate["sd"]

  phi_map_file <- sub("\\.([^.]+)$", ".map.moderated.\\1", phifile)
  phi_map_full_path <- paste0(outdir, '/', phi_map_file)
  if (file.exists(phi_map_full_path)) {
    message("Loading existing MAP phi moderation from: ", phi_map_file)
    phi_map_df <- read.table(phi_map_full_path, header = TRUE, row.names = NULL, sep = "\t")
  }

  if (phi_trend) {
    message("Using loess trend for betabinom_EBmap prior parameters.")
    baseMean_d1r <- sc_d1r %>%
      group_by(gene, event) %>%
      summarise(baseMean = mean(n, na.rm = TRUE), .groups = "drop") %>%
      mutate(comparison = "diff1_vs_ref")
    baseMean_d2r <- sc_d2r %>%
      group_by(gene, event) %>%
      summarise(baseMean = mean(n, na.rm = TRUE), .groups = "drop") %>%
      mutate(comparison = "diff2_vs_ref")
    baseMean_df_all <- bind_rows(baseMean_d1r, baseMean_d2r)

    # EBapprox trend moderation (fallback for MAP failures)
    phi_table_trend <- moderate_phi_trend(phi_df, baseMean_df_all)

    # Residual SD around trend, per comparison.
    # Restrict to identifiable events (when flagged): boundary events sit far
    # above the trend and would inflate the residual SD used as the prior width.
    phi_df_with_trend <- phi_df %>%
      left_join(phi_table_trend %>% select(gene, event, comparison, z_trend),
                by = c("gene", "event", "comparison")) %>%
      filter(!is.na(phi) & phi > 0 & !is.na(z_trend))
    if ("identifiable" %in% names(phi_df_with_trend) &&
        sum(phi_df_with_trend$identifiable %in% TRUE) >= 10) {
      phi_df_with_trend <- phi_df_with_trend %>% filter(identifiable %in% TRUE)
    }
    phi_df_with_trend <- phi_df_with_trend %>%
      mutate(resid = log(phi) - z_trend)

    prior_params <- phi_df_with_trend %>%
      group_by(comparison) %>%
      summarise(sd_resid = sd(resid, na.rm = TRUE), .groups = "drop")

    global_sd_resid <- sd(phi_df_with_trend$resid, na.rm = TRUE)
    if (is.na(global_sd_resid)) global_sd_resid <- sd_logphi
    prior_params <- prior_params %>%
      mutate(sd_resid = ifelse(is.na(sd_resid), global_sd_resid, sd_resid))

    # Augment sc with per-event z_trend and sd_resid
    join_trend_prior <- function(sc) {
      sc %>%
        left_join(phi_table_trend %>% select(gene, event, comparison, z_trend),
                  by = c("gene", "event", "comparison")) %>%
        left_join(prior_params, by = "comparison") %>%
        mutate(z_trend  = ifelse(is.na(z_trend),  mean_logphi,     z_trend),
               sd_resid = ifelse(is.na(sd_resid), global_sd_resid, sd_resid))
    }
    sc_d1r_aug <- join_trend_prior(sc_d1r)
    sc_d2r_aug <- join_trend_prior(sc_d2r)

    if (!exists("phi_map_df")) {
      prior_fn_trend <- function(dd) {
        pd <- data.frame(prior = sprintf("normal(%g,%g)", dd$z_trend[1], dd$sd_resid[1]),
                         class = "fixef_disp", coef = "1", stringsAsFactors = FALSE)
        phi_map_glmmTMB(dd, pd)
      }

      phi_map_d1r <- run_map_for_comparison(sc_d1r_aug, "diff1_vs_ref",
                                            prior_fn_trend, "phi_map.errors.log",
                                            checkpoint_prefix = paste0(outdir, "/ckpt_phimap_", out_prefix, "_d1r"),
                                            mc_cores = mc_cores)
      phi_map_d2r <- run_map_for_comparison(sc_d2r_aug, "diff2_vs_ref",
                                            prior_fn_trend, "phi_map.errors.log",
                                            checkpoint_prefix = paste0(outdir, "/ckpt_phimap_", out_prefix, "_d2r"),
                                            mc_cores = mc_cores)
      phi_map_df  <- bind_rows(phi_map_d1r, phi_map_d2r)
    }

    # Fallback: EBapprox trend z_mod where MAP failed
    phi_fallback <- phi_table_trend %>% select(gene, event, comparison, z_mod)
    fallback_z   <- median(phi_table_trend$z_mod, na.rm = TRUE)

    has_map <- nrow(phi_map_df) > 0 && "z_mod" %in% names(phi_map_df)
    if (!has_map) {
      message("WARNING: no MAP phi estimates obtained; using EBapprox trend fallback. Check phi_map.errors.log.")
    } else {
      message(sprintf("MAP phi trend: %d events with MAP estimate.", nrow(phi_map_df)))
    }
    build_sc_eb_trend <- function(sc_aug) {
      df <- sc_aug %>%
        left_join(phi_fallback, by = c("gene", "event", "comparison"))
      if (has_map) {
        df <- df %>%
          left_join(phi_map_df %>% rename(z_mod_map = z_mod),
                    by = c("gene", "event", "comparison")) %>%
          mutate(z_mod = dplyr::coalesce(z_mod_map, z_mod)) %>%
          select(-z_mod_map)
      }
      df$z_mod[is.na(df$z_mod)] <- fallback_z
      df
    }
    sc_d1r_eb <- build_sc_eb_trend(sc_d1r_aug)
    sc_d2r_eb <- build_sc_eb_trend(sc_d2r_aug)

  } else {
    # EBapprox z_mod as fallback for MAP failures
    phi_approx <- phi_df %>%
      filter(!is.na(comparison)) %>%
      group_by(comparison) %>%
      group_modify(~ moderate_phi_log_scale(.x)) %>%
      ungroup() %>%
      select(gene, event, comparison, z_mod)

    if (!exists("phi_map_df")) {
      # Global prior for all events
      prior_disp <- data.frame(
        prior = sprintf("normal(%g,%g)", mean_logphi, sd_logphi),
        class = "fixef_disp", coef = "1", stringsAsFactors = FALSE
      )
      print(prior_disp)

      prior_fn_global <- function(dd) phi_map_glmmTMB(dd, prior_disp)

      phi_map_d1r <- run_map_for_comparison(sc_d1r, "diff1_vs_ref",
                                            prior_fn_global, "phi_map.errors.log",
                                            checkpoint_prefix = paste0(outdir, "/ckpt_phimap_", out_prefix, "_d1r"),
                                            mc_cores = mc_cores)
      phi_map_d2r <- run_map_for_comparison(sc_d2r, "diff2_vs_ref",
                                            prior_fn_global, "phi_map.errors.log",
                                            checkpoint_prefix = paste0(outdir, "/ckpt_phimap_", out_prefix, "_d2r"),
                                            mc_cores = mc_cores)
      phi_map_df  <- bind_rows(phi_map_d1r, phi_map_d2r)
    }

    # MAP preferred; EBapprox fallback for events where MAP failed
    if (nrow(phi_map_df) > 0 && "z_mod" %in% names(phi_map_df)) {
      message(sprintf("MAP phi: %d events with MAP estimate.", nrow(phi_map_df)))
      phi_final <- phi_approx %>%
        left_join(phi_map_df %>% rename(z_mod_map = z_mod),
                  by = c("gene", "event", "comparison")) %>%
        mutate(z_mod = dplyr::coalesce(z_mod_map, z_mod)) %>%
        select(-z_mod_map)
    } else {
      message("WARNING: no MAP phi estimates obtained; using EBapprox fallback. Check phi_map.errors.log.")
      phi_final <- phi_approx
    }

    sc_d1r_eb <- sc_d1r %>%
      left_join(phi_final, by = c("gene", "event", "comparison")) %>%
      impute_z_mod(phi_final)
    sc_d2r_eb <- sc_d2r %>%
      left_join(phi_final, by = c("gene", "event", "comparison")) %>%
      impute_z_mod(phi_final)
  }

  if (!file.exists(phi_map_full_path)) {
    write.table(phi_map_df, file = phi_map_full_path,
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  sc_both_eb <- bind_rows(sc_d1r_eb, sc_d2r_eb)
  rm(sc_d1r_eb, sc_d2r_eb, phi_df, phi_map_df)
  if (exists("sc_d1r")) rm(sc_d1r, sc_d2r)
  gc()

  results <- run_one_comparison(sc_both_eb,
                                test_fn = test_model_glmmTMB_EB,
                                err_log = "glmmtmb_EBmap.errors.log",
                                model_label = "betabinom_EBmap",
                                L = L,
                                checkpoint_prefix = paste0(outdir, "/ckpt_lrt_", out_prefix, "_EBmap"),
                                mc_cores = mc_cores)
  rm(sc_both_eb)

  results <- results %>%
    left_join(baseMean_df, by = c("gene", "event", "comparison")) %>%
    left_join(lfc_summary_all, by = c("gene", "event", "comparison", "contrast"))

  if (nrow(results) > 0) {
    results <- results %>%
      group_by(contrast) %>%
      group_modify(~ {
        r <- filter_unsupported(.x, .y$contrast); r$pvalue <- r$p.value
        r <- adjust_pvalues(r, independentFiltering = indep_filter, alpha = padj_thr, method = padj_method)
        r$pvalue <- NULL; r
      }) %>%
      ungroup()
  }

  results <- add_significant(add_support(results), padj_thr, delta, min_dpi, min_dpi_sj, min_reads,
                               omnibus_pairs = omnibus_pairs)
  write.table(results, file = out_resultfile, quote = FALSE, sep = "\t", row.names = FALSE)

} else if (model == 'betabinom_EBapprox') {

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
    sc_d1r_eb  <- impute_z_mod(sc_d1r_eb, phi_table)
    sc_d2r_eb  <- impute_z_mod(sc_d2r_eb, phi_table)
  }
  phi_mod_file <- sub("\\.([^.]+)$", ".approx.moderated.\\1", phifile)
  write.table(phi_table, file=paste0(outdir, '/', phi_mod_file), sep="\t", quote=FALSE, row.names=FALSE)

  sc_both_eb <- bind_rows(sc_d1r_eb, sc_d2r_eb)
  rm(phi_df, phi_table, sc_d1r_eb, sc_d2r_eb)
  if (exists("sc_d1r")) rm(sc_d1r, sc_d2r)
  gc()

  results <- run_one_comparison(sc_both_eb,
                                test_fn = test_model_glmmTMB_EB,
                                err_log = "glmmtmb_EB.errors.log",
                                model_label = "betabinom_EBapprox",
                                L = L,
                                checkpoint_prefix = paste0(outdir, "/ckpt_lrt_", out_prefix, "_EBapprox"),
                                mc_cores = mc_cores)
  rm(sc_both_eb)
  results <- results %>%
    left_join(baseMean_df, by = c("gene", "event", "comparison")) %>%
    left_join(lfc_summary_all, by = c("gene", "event", "comparison", "contrast"))

  if (nrow(results) > 0) {
    results <- results %>%
      group_by(contrast) %>%
      group_modify(~ {
        r <- filter_unsupported(.x, .y$contrast); r$pvalue <- r$p.value
        r <- adjust_pvalues(r, independentFiltering = indep_filter, alpha = padj_thr, method = padj_method)
        r$pvalue <- NULL; r
      }) %>%
      ungroup()
  }

  results <- add_significant(add_support(results), padj_thr, delta, min_dpi, min_dpi_sj, min_reads,
                               omnibus_pairs = omnibus_pairs)
  write.table(results, file = out_resultfile, quote = FALSE, sep = "\t", row.names = FALSE)

} else if (model == 'dirmult_EBplugin') {

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
  test_fn_dm <- test_model_multinomial_plugin_dm_EB

  contrast_results <- lapply(seq_along(contrasts), function(ki) {
    ctr_ref  <- contrasts[[ki]]$ref
    ctr_trt  <- contrasts[[ki]]$trt
    ctr_name <- paste0(ctr_trt, "_vs_", ctr_ref)
    ctr_data <- splitcnts_eb %>%
      filter(groups %in% c(ctr_ref, ctr_trt)) %>%
      mutate(groups = factor(as.character(groups), levels = c(ctr_ref, ctr_trt)))
    grouped_ctr <- ctr_data %>% group_by(gene, event) %>% group_split()
    err_log <- paste0(model, ".", ctr_name, ".errors.log")
    res_dm <- mclapply(grouped_ctr, function(dd) {
      tryCatch(test_fn_dm(dd), error = function(e) {
        cat(sprintf("[%s] ERROR gene=%s event=%s (PID %d): %s\n",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    dd$gene[1], dd$event[1], Sys.getpid(), conditionMessage(e)),
            file = err_log, append = TRUE)
        NULL
      })
    }, mc.cores = mc_cores)
    res <- bind_rows(Filter(Negate(is.null), res_dm))
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res$contrast <- ctr_name
    res$pvalue <- res$p.value
    res <- adjust_pvalues(filter_unsupported(res), independentFiltering = FALSE, alpha = padj_thr, method = padj_method)
    res$pvalue <- NULL
    add_significant(add_support(res), padj_thr, delta, min_dpi, min_dpi_sj, min_reads,
                               omnibus_pairs = omnibus_pairs)
  })
  results <- bind_rows(contrast_results)
  write.table(results, file = out_resultfile, sep = "\t", quote = FALSE, row.names = FALSE)

} else if (model == 'wilcoxon') {
  if (!is.null(opt$phi)) warning("--phi is not used for model '", model, "' and will be ignored")

  contrast_results <- lapply(seq_along(contrasts), function(ki) {
    ctr_ref  <- contrasts[[ki]]$ref
    ctr_trt  <- contrasts[[ki]]$trt
    ctr_name <- paste0(ctr_trt, "_vs_", ctr_ref)
    sc_ctr <- bind_rows(comparisons_list) %>%
      filter(groups %in% c(ctr_ref, ctr_trt)) %>%
      mutate(groups = factor(as.character(groups), levels = c(ctr_ref, ctr_trt)))
    res <- run_one_comparison(sc_ctr, test_fn = test_model_wilcoxon,
                              err_log = paste0("wilcoxon.", ctr_name, ".log"),
                              checkpoint_prefix = paste0(outdir, "/ckpt_lrt_", out_prefix, "_wilcoxon_", ctr_name),
                              mc_cores = mc_cores)
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res$contrast <- ctr_name
    res <- res %>%
      left_join(baseMean_df, by = c("gene", "event", "comparison")) %>%
      left_join(lfc_summary_all, by = c("gene", "event", "comparison", "contrast"))
    res$pvalue <- res$p.value
    res <- adjust_pvalues(filter_unsupported(res), independentFiltering = indep_filter, alpha = padj_thr, method = padj_method)
    res$pvalue <- NULL
    add_significant(add_support(res), padj_thr, delta, min_dpi, min_dpi_sj, min_reads,
                               omnibus_pairs = omnibus_pairs)
  })
  results <- bind_rows(contrast_results)
  write.table(results, file = out_resultfile, quote = FALSE, sep = "\t", row.names = FALSE)

} else if (model == 'betabinom_MLE') {
  if (!is.null(opt$phi)) warning("--phi is not used for model '", model, "' and will be ignored")

  results <- run_one_comparison(bind_rows(comparisons_list),
                                test_fn = function(dd, ...) suppressMessages(test_model_glmmTMB_without_prior(dd, ...)),
                                err_log = "glmmtmb_noprior.errors.log",
                                L = L,
                                checkpoint_prefix = paste0(outdir, "/ckpt_lrt_", out_prefix, "_MLE"),
                                mc_cores = mc_cores)
  results <- results %>%
    left_join(baseMean_df, by = c("gene", "event", "comparison")) %>%
    left_join(lfc_summary_all, by = c("gene", "event", "comparison", "contrast"))

  if (nrow(results) > 0) {
    results <- results %>%
      group_by(contrast) %>%
      group_modify(~ {
        r <- filter_unsupported(.x, .y$contrast); r$pvalue <- r$p.value
        r <- adjust_pvalues(r, independentFiltering = indep_filter, alpha = padj_thr, method = padj_method)
        r$pvalue <- NULL; r
      }) %>%
      ungroup()
  }

  results <- add_significant(add_support(results), padj_thr, delta, min_dpi, min_dpi_sj, min_reads,
                               omnibus_pairs = omnibus_pairs)
  write.table(results, file = out_resultfile, quote = FALSE, sep = "\t", row.names = FALSE)

} else {
  print ("model misspecified")
  print ("choose between betabinom_EBmap, betabinom_EBapprox, betabinom_MLE, dirmult_EBplugin, wilcoxon")

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
message(paste("Using", mc_cores, "cores for reading split files."))
split_list <- mclapply(unique_genes, read_gene_split, mc.cores = mc_cores)
# Filter out NULLs
split_list <- split_list[!sapply(split_list, is.null)]
if (length(split_list) == 0) {
  stop("No matching split files found for any genes in the test file.")
}
splits <- bind_rows(split_list)
rm(split_list)
gc()
message(paste("Loaded split data for", length(unique(splits$gene)), "genes."))


# Combine data
# Using left_join to keep all tests, filling NAs for missing annotations
message("Merging data...")
merged_data <- left_join(tests, splits, by = c("gene", "event"))
message(paste("Merged dataset has", nrow(merged_data), "rows."))

# `significant` was set before the annotation merge, where setdiff1/setdiff2 did
# not exist yet, so every row got the exonic min_dpi. Now that the merge has
# attached them, recompute so junction-sourced sides use min_dpi_sj.
if (nrow(merged_data) > 0 && all(c("setdiff1", "setdiff2") %in% names(merged_data))) {
  merged_data <- add_significant(add_support(merged_data), padj_thr, delta, min_dpi, min_dpi_sj, min_reads,
                               omnibus_pairs = omnibus_pairs)
}

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
      group_by(gene, event, contrast) %>%
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
      group_by(gene, event, contrast) %>%
      slice(which.min(p.value)) %>%
      ungroup()

    min_data <- primary_df %>%
      left_join(combo_df, by = c("gene", "event", "contrast")) %>%
      mutate(p.value = p_min) %>%
      select(-p_min)

    if (nrow(min_data) > 0) {
      min_data$pvalue <- min_data$p.value
      min_data <- adjust_pvalues(filter_unsupported(min_data), independentFiltering = indep_filter, alpha = padj_thr, method = padj_method)
      min_data$pvalue <- NULL
    }

    min_data <- add_significant(add_support(min_data), padj_thr, delta, min_dpi, min_dpi_sj, min_reads,
                               omnibus_pairs = omnibus_pairs)
    out_mincomb <- sub("\\.annotated\\.txt$", ".mincomb.annotated.txt",
                      out_result_annotated)
    write.table(min_data, out_mincomb, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste("Min-p combined results written to", out_mincomb))

    # Fisher p-value combination for bipartition_both:
    # Combine p-values from _s1 and _s2 using Fisher's method:
    #   X = -2 * sum(log(p))  ~  chi-squared with 2*k df  (k = number of sides)
    # When only one side has a valid p-value, use that p directly (1 df = 2).
    fisher_combo_df <- merged_data %>%
      group_by(gene, event, contrast) %>%
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
      left_join(fisher_combo_df, by = c("gene", "event", "contrast")) %>%
      mutate(p.value = p_fisher) %>%
      select(-p_fisher)

    if (nrow(fisher_data) > 0) {
      fisher_data$pvalue <- fisher_data$p.value
      fisher_data <- adjust_pvalues(filter_unsupported(fisher_data), independentFiltering = indep_filter, alpha = padj_thr, method = padj_method)
      fisher_data$pvalue <- NULL
    }

    fisher_data <- add_significant(add_support(fisher_data), padj_thr, delta, min_dpi, min_dpi_sj, min_reads,
                               omnibus_pairs = omnibus_pairs)
    out_fisher <- sub("\\.annotated\\.txt$", ".fisher_combined.annotated.txt",
                      out_result_annotated)
    write.table(fisher_data, out_fisher, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste("Fisher combined results written to", out_fisher))
  }
}
      
