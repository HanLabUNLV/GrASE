#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-g", "--grase"), type = "character", default = "~/GrASE_simulation/bipartition.test/test_bipartition.internal_glmmTMB_fixed_EB.annotated.txt",
              help = "Path to GrASE annotated exontest results file [default= %default]", metavar = "character"),
  make_option(c("--grase2"), type = "character", default = "~/GrASE_simulation/bipartition.test/test_bipartition.TSSTTS_glmmTMB_fixed_EB.annotated.txt",
              help = "Path to GrASE second optional annotated exontest results file [default= %default]", metavar = "character"),
  make_option(c("-s", "--saturn"), type = "character", default = "~/GrASE_simulation/saturn/saturn_group1_group2/all.group1_group2.sumExp_filteredbyCountMultiExon.txt",
              help = "Path to Saturn results file [default= %default]", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "saturn_grase_glmmTMB_fixed_EB.comparison.txt",
              help = "Output file path [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

grase_file <- opt$grase
grase2_file <- opt$grase2
saturn_file <- opt$saturn
output_file <- opt$output

# Generate PDF filename from output filename
plot_file <- sub("\\.[^.]+$", ".pdf", output_file)
if (plot_file == output_file) {
  plot_file <- paste0(output_file, ".pdf")
}

message("Comparing Saturn and grase results.")
message("grase file: ", grase_file)
message("grase2 file: ", grase2_file)
message("Saturn file:  ", saturn_file)

if (!file.exists(grase_file)) {
  stop(paste("grase file not found:", grase_file))
}
if (!file.exists(grase2_file)) {
  stop(paste("grase2 file not found:", grase2_file))
}
if (!file.exists(saturn_file)) {
  stop(paste("Saturn file not found:", saturn_file))
}

# --- Load Data ---
grase_dt <- fread(grase_file)
grase2_dt <- fread(grase2_file)
grase_dt <- rbind(grase_dt, grase2_dt, use.names = TRUE, fill = TRUE)
rm(grase2_dt)
saturn_dt <- fread(saturn_file)

# --- Process grase ---
# 1. Select relevant columns including 'setdiff' which contains the matching exon IDs
# 2. Expand comma-separated 'setdiff' (values like "E006,E005" -> 2 rows)
# 3. For each (gene, setdiff) pair, keep the result with the SMALLEST p-value

message("Processing grase data...")
grase_processed <- grase_dt %>%
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
    grase_p.value = p.value,
    grase_padj = padj,
    grase_event = event,
    grase_effect_size = effect_size,
    grase_refexpart = ref_ex_part
  ) %>%
  mutate(grase_setdiff = exon)

message(" - Unique gene-exon pairs in grase: ", nrow(grase_processed))

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
merged_df <- inner_join(grase_processed, saturn_processed, by = c("gene", "exon"))

message("Merged dataset size: ", nrow(merged_df))

# --- Summarize ---
stats_summary <- list(
  total_common_pairs = nrow(merged_df),
  
  # Correlation of p-values (log scale often more useful)
  pearson_log_pval = cor(-log10(merged_df$grase_p.value + 1e-300), 
                         -log10(merged_df$saturn_pval + 1e-300), 
                         method = "pearson", use = "complete.obs"),
  spearman_pval = cor(merged_df$grase_p.value, merged_df$saturn_pval, method = "spearman", use = "complete.obs"),
  
  # Correlation of empirical p-values
  pearson_log_empirical_pval = cor(-log10(merged_df$grase_p.value + 1e-300), 
                                   -log10(merged_df$saturn_empirical_pval + 1e-300), 
                                   method = "pearson", use = "complete.obs"),
  spearman_empirical_pval = cor(merged_df$grase_p.value, merged_df$saturn_empirical_pval, method = "spearman", use = "complete.obs"),

  # Correlation of effect sizes
  pearson_effect_size = cor(merged_df$grase_effect_size, merged_df$saturn_estimate, method = "pearson", use = "complete.obs"),
  spearman_effect_size = cor(merged_df$grase_effect_size, merged_df$saturn_estimate, method = "spearman", use = "complete.obs"),

  # Significant counts (p < 0.05)
  n_sig_grase = sum(merged_df$grase_p.value < 0.05, na.rm = TRUE),
  n_sig_saturn = sum(merged_df$saturn_pval < 0.05, na.rm = TRUE),
  n_sig_both = sum(merged_df$grase_p.value < 0.05 & merged_df$saturn_pval < 0.05, na.rm = TRUE)
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
      log_grase = -log10(grase_p.value + 1e-300),
      log_saturn = -log10(saturn_pval + 1e-300),
      log_saturn_empirical = -log10(saturn_empirical_pval + 1e-300)
    )
  
  # 1. log P-value correlation
  p1 <- ggplot(df, aes(x = log_grase, y = log_saturn)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
    labs(
      title = "P-value Correlation (-log10)",
      x = "-log10(grase p-value)",
      y = "-log10(Saturn p-value)"
    ) +
    theme_minimal()
  
  # 2. log Empirical P-value correlation
  p2 <- ggplot(df, aes(x = log_grase, y = log_saturn_empirical)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
    labs(
      title = "Empirical P-value Correlation (-log10)",
      x = "-log10(grase p-value)",
      y = "-log10(Saturn empirical p-value)"
    ) +
    theme_minimal()
  
  # 3. Effect size correlation
  p3 <- ggplot(df, aes(x = grase_effect_size, y = saturn_estimate)) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
    labs(
      title = "Effect Size Correlation",
      x = "grase Effect Size",
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
