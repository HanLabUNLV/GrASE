#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-g", "--glmm"), type = "character", default = "test_glmmTMB_fixed_EB.annotated.txt",
              help = "Path to glmmTMB annotated results file [default= %default]", metavar = "character"),
  make_option(c("-s", "--saturn"), type = "character", default = "../saturn_CD4_CD4_N_STIM/all.CD4_CD4_N_STIM.sumExp_filteredbyCountMultiExon.txt",
              help = "Path to Saturn results file [default= %default]", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "saturn_glmmTMB_comparison.txt",
              help = "Output file path [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

glmm_file <- opt$glmm
saturn_file <- opt$saturn
output_file <- opt$output

# Generate PDF filename from output filename
plot_file <- sub("\\.[^.]+$", ".pdf", output_file)
if (plot_file == output_file) {
  plot_file <- paste0(output_file, ".pdf")
}

message("Comparing Saturn and glmmTMB results.")
message("glmmTMB file: ", glmm_file)
message("Saturn file:  ", saturn_file)

if (!file.exists(glmm_file)) {
  stop(paste("glmmTMB file not found:", glmm_file))
}
if (!file.exists(saturn_file)) {
  stop(paste("Saturn file not found:", saturn_file))
}

# --- Load Data ---
glmm_dt <- fread(glmm_file)
saturn_dt <- fread(saturn_file)

# --- Process glmmTMB ---
# 1. Select relevant columns including 'setdiff' which contains the matching exon IDs
# 2. Expand comma-separated 'setdiff' (values like "E006,E005" -> 2 rows)
# 3. For each (gene, setdiff) pair, keep the result with the SMALLEST p-value

message("Processing glmmTMB data...")
glmm_processed <- glmm_dt %>%
  select(gene, setdiff, p.value, event, model, padj, diff_mean, ref_mean, effect_size, ref_ex_part) %>%
  mutate(setdiff = as.character(setdiff)) %>%
  separate_rows(setdiff, sep = ",") %>%
  mutate(setdiff = trimws(setdiff)) %>%
  filter(!is.na(p.value)) %>%
  group_by(gene, setdiff) %>%
  slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(
    exon = setdiff,
    glmmTMB_p.value = p.value,
    glmmTMB_padj = padj,
    glmmTMB_event = event,
    glmmTMB_effect_size = effect_size,
    glmmtmb_refexpart = ref_ex_part
  ) %>%
  mutate(glmmtmb_setdiff = exon)

message(" - Unique gene-exon pairs in glmmTMB: ", nrow(glmm_processed))

# --- Process Saturn ---
# 1. Split 'ids' (format "gene:exon") into 'gene' and 'exon' columns
message("Processing Saturn data...")
saturn_processed <- saturn_dt %>%
  separate(ids, into = c("gene", "exon"), sep = ":", remove = FALSE) %>%
  select(gene, exon, pval, regular_FDR, estimates, se, empirical_pval) %>%
  rename(
    saturn_pval = pval,
    saturn_FDR = regular_FDR,
    saturn_estimate = estimates,
    saturn_empirical_pval = empirical_pval
  )

message(" - Unique gene-exon pairs in Saturn: ", nrow(saturn_processed))

# --- Match and Compare ---
# Inner join to find common gene-exon pairs
merged_df <- inner_join(glmm_processed, saturn_processed, by = c("gene", "exon"))

message("Merged dataset size: ", nrow(merged_df))

# --- Summarize ---
stats_summary <- list(
  total_common_pairs = nrow(merged_df),
  
  # Correlation of p-values (log scale often more useful)
  pearson_log_pval = cor(-log10(merged_df$glmmTMB_p.value + 1e-300), 
                         -log10(merged_df$saturn_pval + 1e-300), 
                         method = "pearson", use = "complete.obs"),
  spearman_pval = cor(merged_df$glmmTMB_p.value, merged_df$saturn_pval, method = "spearman", use = "complete.obs"),
  
  # Correlation of empirical p-values
  pearson_log_empirical_pval = cor(-log10(merged_df$glmmTMB_p.value + 1e-300), 
                                   -log10(merged_df$saturn_empirical_pval + 1e-300), 
                                   method = "pearson", use = "complete.obs"),
  spearman_empirical_pval = cor(merged_df$glmmTMB_p.value, merged_df$saturn_empirical_pval, method = "spearman", use = "complete.obs"),

  # Correlation of effect sizes
  pearson_effect_size = cor(merged_df$glmmTMB_effect_size, merged_df$saturn_estimate, method = "pearson", use = "complete.obs"),
  spearman_effect_size = cor(merged_df$glmmTMB_effect_size, merged_df$saturn_estimate, method = "spearman", use = "complete.obs"),

  # Significant counts (p < 0.05)
  n_sig_glmmTMB = sum(merged_df$glmmTMB_p.value < 0.05, na.rm = TRUE),
  n_sig_saturn = sum(merged_df$saturn_pval < 0.05, na.rm = TRUE),
  n_sig_both = sum(merged_df$glmmTMB_p.value < 0.05 & merged_df$saturn_pval < 0.05, na.rm = TRUE)
)

print("--- Summary Statistics ---")
print(stats_summary)

# --- Plotting ---
message("Generating plots...")
generate_plots <- function(df, output_pdf) {
  pdf(output_pdf, width = 10, height = 5)
  
  # Prepare for log-log plot (add tiny constant to avoid log(0))
  df <- df %>%
    mutate(
      log_glmm = -log10(glmmTMB_p.value + 1e-300),
      log_saturn = -log10(saturn_pval + 1e-300),
      log_saturn_empirical = -log10(saturn_empirical_pval + 1e-300)
    )
  
  # 1. log P-value correlation
  p1 <- ggplot(df, aes(x = log_glmm, y = log_saturn)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
    labs(
      title = "P-value Correlation (-log10)",
      x = "-log10(glmmTMB p-value)",
      y = "-log10(Saturn p-value)"
    ) +
    theme_minimal()
  
  # 2. log Empirical P-value correlation
  p2 <- ggplot(df, aes(x = log_glmm, y = log_saturn_empirical)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
    labs(
      title = "Empirical P-value Correlation (-log10)",
      x = "-log10(glmmTMB p-value)",
      y = "-log10(Saturn empirical p-value)"
    ) +
    theme_minimal()
  
  # 3. Effect size correlation
  p3 <- ggplot(df, aes(x = glmmTMB_effect_size, y = saturn_estimate)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
    labs(
      title = "Effect Size Correlation",
      x = "glmmTMB Effect Size",
      y = "Saturn Estimate"
    ) +
    theme_minimal()
  
  # Arrange plots side-by-side (using gridExtra if available, or just printing pages)
  # Since gridExtra is not guaranteed, printing separately to the PDF
  print(p1)
  print(p2)
  print(p3)
  
  dev.off()
}

generate_plots(merged_df, plot_file)
message("Plots saved to: ", plot_file)

# --- Save ---
write.table(merged_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Results saved to: ", output_file)
