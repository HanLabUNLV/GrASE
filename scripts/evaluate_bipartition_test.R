#!/usr/bin/env Rscript

# scripts/evaluate_bipartition_test.R
#
# Evaluates GrASE bipartition statistical test results against simulation
# ground truth.
#
# Input test results file columns (tab-separated):
#   gene, event, LRT, p.value, model, phi, effect_size, padj,
#   source, sink, ref_ex_part, setdiff1, setdiff2,
#   transcripts1, transcripts2, path1, path2,
#   ref_mean, diff1_mean, diff2_mean, diff_mean, which, setdiff
#
# Ground truth (gt_dir, from infer_diff_exons_gt.R):
#   gene.exonic_parts_fc.txt   exonic_part, fold_change, group, transcripts,
#                               changed_tx, alt_tx, source, sink
#
# Level 1 - Exonic part detection:
#   Detected positive = exonic parts in `setdiff` column with padj < threshold
#   GT positive       = exonic parts labelled "changed" or "constant"
#   GT negative       = exonic parts labelled "negative"
#   Reports TP/FP/FN/TN, precision, recall, F1 at multiple padj thresholds
#
#   Two evaluations:
#     Full (total space)   : universe = all GT genes; untested = implicit neg
#     Restricted           : universe = only GrASE-testable exons (setdiff)
#
# Usage:
#   Rscript evaluate_bipartition_test.R <test_results> <gt_dir> <out_dir> <simulate_rda>
#
#   <test_results> may be a single file path or a comma-separated list of paths.
#   When multiple files are given they are combined: for Level 1 the minimum
#   padj across files is used per (gene, setdiff).
#
# Outputs:
#   level1_per_gene_padj<thr>.txt                  per-gene metrics (full GT)
#   level1_restricted_per_gene_padj<thr>.txt       per-gene metrics (restricted)
#   level1_summary_by_simtype.txt                  P/R/F1 by sim_type (full GT)
#   level1_restricted_summary_by_simtype.txt       same, restricted

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript evaluate_bipartition_test.R <test_results> <gt_dir> <out_dir> <simulate_rda>\n")
  quit(save = "no", status = 1)
}

test_files  <- trimws(unlist(strsplit(args[1], ",")))
gt_dir      <- args[2]
out_dir     <- args[3]
sim_rda     <- args[4]

padj_thresholds <- c(0.01, 0.05, 0.1, 0.2)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── helpers ───────────────────────────────────────────────────────────────────

parse_exparts <- function(x) {
  if (is.na(x) || trimws(x) == "" || x == "NA") return(character(0))
  trimws(unlist(strsplit(x, ",")))
}

f1_safe <- function(p, r) {
  if (is.na(p) || is.na(r) || (p + r) == 0) NA_real_
  else 2 * p * r / (p + r)
}

# ── load test results ─────────────────────────────────────────────────────────

cat(sprintf("Loading %d test file(s)...\n", length(test_files)))
tests_list <- lapply(test_files, function(f) {
  cat(sprintf("  %s\n", basename(f)))
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             quote = "", comment.char = "", na.strings = c("NA", ""),
             colClasses = c(source = "character", sink = "character"))
})
tests <- bind_rows(tests_list)
cat(sprintf("  %d total rows, %d unique genes across all files\n",
            nrow(tests), length(unique(tests$gene))))

# When multiple files are combined, take the minimum padj per (gene, setdiff)
if (length(test_files) > 1) {
  tests_l1 <- tests %>%
    group_by(gene, setdiff) %>%
    summarise(padj = min(padj, na.rm = TRUE), .groups = "drop")
} else {
  tests_l1 <- tests %>% select(gene, setdiff, padj)
}

# Build per-gene set of exonic parts GrASE explicitly tested (setdiff only)
grase_testable_by_gene <- lapply(
  split(tests_l1$setdiff, tests_l1$gene),
  function(x) unique(x[!is.na(x) & nchar(x) > 0])
)
cat(sprintf("  GrASE: %d genes with testable exonic parts (setdiff only)\n",
            length(grase_testable_by_gene)))

# ── load simulation gene type labels ─────────────────────────────────────────

cat("Loading simulate.rda...\n")
dge_genes <- dte_genes <- dtu_genes <- character(0)
if (file.exists(sim_rda)) {
  load(sim_rda)
  if (exists("dge.genes")) dge_genes <- dge.genes
  if (exists("dte.genes")) dte_genes <- dte.genes
  if (exists("dtu.genes")) dtu_genes <- dtu.genes
  cat(sprintf("  DGE: %d  DTE: %d  DTU: %d genes\n",
              length(dge_genes), length(dte_genes), length(dtu_genes)))
} else {
  warning("simulate.rda not found; sim_type will be 'Unknown' for all genes")
}

get_sim_type <- function(gene) {
  if (gene %in% dge_genes) return("DGE")
  if (gene %in% dte_genes) return("DTE")
  if (gene %in% dtu_genes) return("DTU")
  return("Background")
}

# ── load ground truth ─────────────────────────────────────────────────────────

cat("Loading ground truth exonic_parts_fc files...\n")
gt_fc_files <- list.files(gt_dir, pattern = "\\.exonic_parts_fc\\.txt$",
                           full.names = TRUE)
gt_genes_available <- sub("\\.exonic_parts_fc\\.txt$", "",
                          basename(gt_fc_files))

test_genes <- unique(tests$gene)
eval_genes <- gt_genes_available   # all GT genes (total-space universe)
cat(sprintf("  GT genes: %d  |  Test genes: %d  |  GT-only (untested): %d\n",
            length(gt_genes_available), length(test_genes),
            length(setdiff(gt_genes_available, test_genes))))

gt_all <- lapply(gt_fc_files, function(f) {
  tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = "",
                      colClasses = c(exonic_part = "character",
                                     fold_change = "character",
                                     source      = "character",
                                     sink        = "character")),
           error = function(e) NULL)
})
gt_all <- bind_rows(gt_all[!sapply(gt_all, is.null)])
cat(sprintf("  Loaded %d exonic part records across %d genes\n",
            nrow(gt_all), length(unique(gt_all$gene))))

gt_all$sim_type <- sapply(gt_all$gene, get_sim_type)

# ── Level 1: exonic part detection ───────────────────────────────────────────

cat("\nRunning Level 1 (exonic part detection)...\n")

add_all_row <- function(summary_df) {
  all_rows <- summary_df %>%
    group_by(padj_thr) %>%
    summarise(
      sim_type        = "ALL",
      n_genes         = sum(n_genes),
      total_TP        = sum(total_TP),
      total_FP        = sum(total_FP),
      total_FN        = sum(total_FN),
      total_TN        = sum(total_TN),
      micro_precision = round(total_TP / max(total_TP + total_FP, 1), 4),
      micro_recall    = round(total_TP / max(total_TP + total_FN, 1), 4),
      micro_f1        = round(f1_safe(micro_precision, micro_recall), 4),
      macro_precision = NA_real_,
      macro_recall    = NA_real_,
      macro_f1        = NA_real_,
      .groups = "drop"
    )
  bind_rows(summary_df, all_rows) %>% arrange(sim_type, padj_thr)
}

run_level1 <- function(label, testable_by_gene = NULL, restrict_pos = TRUE) {
  cat(sprintf("  [%s]\n", label))
  all_thresholds <- list()

  for (thr in padj_thresholds) {
    cat(sprintf("    padj < %.2f\n", thr))

    sig_rows    <- tests_l1[!is.na(tests_l1$padj) & tests_l1$padj < thr, ]
    det_by_gene <- split(sig_rows$setdiff, sig_rows$gene)

    rows <- lapply(eval_genes, function(gene) {
      gt_gene <- gt_all[gt_all$gene == gene, ]
      if (nrow(gt_gene) == 0) return(NULL)

      testable <- if (!is.null(testable_by_gene)) testable_by_gene[[gene]] else NULL
      if (!is.null(testable_by_gene) &&
          (is.null(testable) || length(testable) == 0)) return(NULL)

      gt_pos <- unique(gt_gene$exonic_part[
        gt_gene$group %in% c("changed", "constant")])
      gt_neg <- unique(gt_gene$exonic_part[gt_gene$group == "negative"])
      if (!is.null(testable)) {
        if (restrict_pos) gt_pos <- intersect(gt_pos, testable)
        gt_neg <- intersect(gt_neg, testable)
      }

      sim_type    <- gt_gene$sim_type[1]
      det_exparts <- unique(det_by_gene[[gene]])
      det_exparts <- det_exparts[!is.na(det_exparts) & nchar(det_exparts) > 0]

      TP <- length(intersect(det_exparts, gt_pos))
      FP <- length(intersect(det_exparts, gt_neg))
      FN <- length(setdiff(gt_pos, det_exparts))
      TN <- length(setdiff(gt_neg, det_exparts))

      precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
      recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
      f1        <- f1_safe(precision, recall)

      data.frame(
        gene       = gene,
        sim_type   = sim_type,
        padj_thr   = thr,
        n_gt_pos   = length(gt_pos),
        n_gt_neg   = length(gt_neg),
        n_detected = length(det_exparts),
        TP = TP, FP = FP, FN = FN, TN = TN,
        precision  = round(precision, 4),
        recall     = round(recall,    4),
        f1         = round(f1,        4),
        stringsAsFactors = FALSE
      )
    })
    level1_df <- bind_rows(rows[!sapply(rows, is.null)])
    all_thresholds[[as.character(thr)]] <- level1_df

    write.table(level1_df,
                file.path(out_dir, sprintf("%s_per_gene_padj%.2f.txt",
                                          label, thr)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  summary_rows <- lapply(names(all_thresholds), function(thr_str) {
    df  <- all_thresholds[[thr_str]]
    thr <- as.numeric(thr_str)
    df %>%
      group_by(sim_type) %>%
      summarise(
        padj_thr        = thr,
        n_genes         = n(),
        total_TP        = sum(TP,  na.rm = TRUE),
        total_FP        = sum(FP,  na.rm = TRUE),
        total_FN        = sum(FN,  na.rm = TRUE),
        total_TN        = sum(TN,  na.rm = TRUE),
        micro_precision = round(sum(TP) / max(sum(TP) + sum(FP), 1), 4),
        micro_recall    = round(sum(TP) / max(sum(TP) + sum(FN), 1), 4),
        micro_f1        = round(f1_safe(
                            sum(TP) / max(sum(TP) + sum(FP), 1),
                            sum(TP) / max(sum(TP) + sum(FN), 1)), 4),
        macro_precision = round(mean(precision, na.rm = TRUE), 4),
        macro_recall    = round(mean(recall,    na.rm = TRUE), 4),
        macro_f1        = round(mean(f1,        na.rm = TRUE), 4),
        .groups = "drop"
      )
  })
  summary1 <- bind_rows(summary_rows) %>% arrange(sim_type, padj_thr)
  summary1  <- add_all_row(summary1)

  write.table(summary1,
              file.path(out_dir, sprintf("%s_summary_by_simtype.txt", label)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  list(thresholds = all_thresholds, summary = summary1)
}

# Full evaluation (total space): universe = all GT exons
l1_full  <- run_level1("level1")
# Coverage-restricted: universe = GrASE-testable exons (setdiff only)
l1_restr <- run_level1("level1_restricted",
                        testable_by_gene = grase_testable_by_gene,
                        restrict_pos = TRUE)

cat("\n=== Level 1 Summary (full GT, micro metrics) ===\n")
print(as.data.frame(l1_full$summary %>%
  select(sim_type, padj_thr, n_genes,
         micro_precision, micro_recall, micro_f1,
         total_TP, total_FP, total_FN)))

cat("\n=== Level 1 Summary (restricted to GrASE-testable exons) ===\n")
print(as.data.frame(l1_restr$summary %>%
  select(sim_type, padj_thr, n_genes,
         micro_precision, micro_recall, micro_f1,
         total_TP, total_FP, total_FN)))

cat(sprintf("\nResults written to: %s\n", out_dir))
