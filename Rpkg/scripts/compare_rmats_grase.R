#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

# --- Command Line Options ---
option_list <- list(
  make_option(c("-g", "--grase"), type = "character", 
              help = "Path to GrASE annotated results [default= %default]", metavar = "character"),
  make_option(c("-r", "--rmats_dir"), type = "character", 
              help = "Directory containing rMATS results (SE.MATS.JC.txt etc) [default= %default]", metavar = "character"),
  make_option(c("-m", "--map_dir"), type = "character", 
              help = "Directory containing GrASE-rMATS mapping files [default= %default]", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", 
              help = "Output file path [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Expand paths
grase_file <- path.expand(opt$grase)
rmats_dir <- path.expand(opt$rmats_dir)
map_dir <- path.expand(opt$map_dir)
output_file <- opt$output
plot_file <- sub("\\.[^.]+$", ".pdf", output_file)
if (plot_file == output_file) plot_file <- paste0(output_file, ".pdf")

message("--- Configuration ---")
message("GrASE File:    ", grase_file)
message("rMATS Dir:     ", rmats_dir)
message("Mapping Dir:   ", map_dir)
message("Output File:   ", output_file)

# --- Validation ---
if (!file.exists(grase_file)) stop("GrASE file not found: ", grase_file)
if (!dir.exists(rmats_dir)) stop("rMATS directory not found: ", rmats_dir)
if (!dir.exists(map_dir)) stop("Mapping directory not found: ", map_dir)

# --- Load GrASE Data ---
message("\nLoading GrASE results...")
grase_dt <- fread(grase_file)
# Enforce string types for join keys
grase_dt <- grase_dt %>% 
  mutate(
    gene = as.character(gene),
    event = as.character(event)
  )

message("GrASE rows: ", nrow(grase_dt))

# --- Process Each rMATS Event Type ---
event_types <- c("SE", "MXE", "A3SS", "A5SS", "RI")
merged_list <- list()

for (etype in event_types) {
  # 1. Look for mapping file (now combined.fromGTF.*.txt)
  map_filename <- paste0("combined.fromGTF.", etype, ".txt")
  map_path <- file.path(map_dir, map_filename)
  
  # 2. Look for rMATS result file
  # Trying JC.txt (Junction Counts only) or JCEC.txt (Junction Counts + Exon Counts)
  rmats_filename <- paste0(etype, ".MATS.JC.txt")
  rmats_path <- file.path(rmats_dir, rmats_filename)
  
  if (!file.exists(map_path)) {
    message("Skipping ", etype, ": Mapping file not found (", map_path, ")")
    next
  }
  if (!file.exists(rmats_path)) {
    message("Skipping ", etype, ": rMATS file not found (", rmats_path, ")")
    next
  }
  
  message("Processing ", etype, "...")
  
  # Load Mapping
  # Format of combined.fromGTF.SE.txt:
  # ID, GeneID, geneSymbol, chr, strand, longExonStart_0base, longExonEnd, shortES, shortEE, flankingES, flankingEE, DexseqFragment, DexseqRefFrag, bipartID
  map_dt <- fread(map_path)
  if (nrow(map_dt) == 0) {
    message("  Empty mapping file.")
    next
  }
  
  # Check if 'bipartID' exists (as per user instruction)
  if (!"bipartID" %in% colnames(map_dt)) {
      message("  'bipartID' column not found in ", map_filename)
      # Check if 'DexseqFragment' is the intended column if bipartID is missing, but strict request asked for bipartID
      if ("DexseqFragment" %in% colnames(map_dt)) {
          message("  Falling back to 'DexseqFragment' as event ID")
          map_dt$bipartID <- map_dt$DexseqFragment
      } else {
          message("  Skipping due to missing ID column")
          next
      }
  }

  # Prepare mapping table
  # We join with GrASE on (GeneID <-> gene, bipartID <-> event)
  # And join with rMATS on (ID <-> ID)
  # Keep all columns from map_dt to check mapping
  map_clean <- map_dt %>%
    mutate(
      GeneID_map = as.character(GeneID),
      bipartID_map = as.character(bipartID),
      ID_map = as.character(ID)     # The rMATS ID
    ) %>%
    filter(!is.na(ID_map) & ID_map != "NA")
  
  # Load rMATS
  rmats_dt <- fread(rmats_path)
  names(rmats_dt) <- make.unique(names(rmats_dt))
  
  # rMATS ID column is usually "ID"
  rmats_clean <- rmats_dt %>%
    mutate(ID = as.character(ID)) %>%
    select(ID, PValue, FDR, IncLevelDifference) %>%
    rename(
      rmats_pval = PValue,
      rmats_fdr = FDR,
      rmats_inc_diff = IncLevelDifference
    )
  
  # Merge: Mapping + rMATS (on rMATS 'ID')
  mapped_rmats <- inner_join(map_clean, rmats_clean, by = c("ID_map" = "ID"))
  
  if (nrow(mapped_rmats) == 0) {
    message("  No overlap between mapping IDs and rMATS file IDs.")
    next
  }
  
  # Merge: (Mapping + rMATS) + GrASE (on gene + event)
  final_merged <- inner_join(grase_dt, mapped_rmats, by = c("gene" = "GeneID_map", "event" = "bipartID_map")) %>%
    mutate(rMATS_EventType = etype)
  
  message("  Matches found: ", nrow(final_merged))
  if (nrow(final_merged) > 0) {
    merged_list[[etype]] <- final_merged
  }
}

# --- Combine All Types ---
if (length(merged_list) == 0) {
  stop("No common events found across any types.")
}

full_data <- rbindlist(merged_list, fill = TRUE)
message("\nTotal matched events: ", nrow(full_data))

# --- Metrics ---
# GrASE uses effect_size, rMATS uses IncLevelDifference
# GrASE uses glmmTMB_p.value (based on earlier column renaming, or p.value if standard file)
# Check columns
p_col <- if ("glmmTMB_p.value" %in% colnames(full_data)) "glmmTMB_p.value" else "p.value"
eff_col <- if ("glmmTMB_effect_size" %in% colnames(full_data)) "glmmTMB_effect_size" else "effect_size"

if (!p_col %in% colnames(full_data)) {
  message("Warning: Could not find p-value column (glmmTMB_p.value or p.value). Using first numeric column?")
  # Fallback logic or stop? Let's just create NAs if missing
  full_data$p.value <- NA
  p_col <- "p.value"
}
if (!eff_col %in% colnames(full_data)) {
  message("Warning: Could not find effect size column (glmmTMB_effect_size or effect_size).")
  full_data$effect_size <- NA
  eff_col <- "effect_size"
}

# Create classification BEFORE correlations so we can filter
SIG_CUTOFF <- 0.05
full_data <- full_data %>%
  mutate(
    grase_p = get(p_col),
    grase_sig = !is.na(grase_p) & grase_p < SIG_CUTOFF,
    rmats_sig = !is.na(rmats_pval) & rmats_pval < SIG_CUTOFF,
    category = case_when(
      grase_sig & rmats_sig ~ "Both_Sig",
      grase_sig & !rmats_sig ~ "GrASE_Only",
      !grase_sig & rmats_sig ~ "rMATS_Only",
      TRUE ~ "Neither_Sig"
    ),
    dir_agree = sign(get(eff_col)) == sign(rmats_inc_diff)
  )

# --- Function for Stats ---
calc_stats <- function(dt, p_col, eff_col, label="Dataset") {
  list(
    label = label,
    n = nrow(dt),
    pearson_log_pval = cor(-log10(dt[[p_col]] + 1e-300), 
                           -log10(dt$rmats_pval + 1e-300), 
                           method = "pearson", use = "complete.obs"),
    spearman_pval = cor(dt[[p_col]], dt$rmats_pval, method = "spearman", use = "complete.obs"),
    pearson_effect = cor(dt[[eff_col]], dt$rmats_inc_diff, method = "pearson", use = "complete.obs")
  )
}

# 1. Total Dataset Stats
stats_total <- calc_stats(full_data, p_col, eff_col, "TOTAL")
print("--- Summary Stats (Total) ---")
print(stats_total)

# 2. Filtered Dataset Stats (Removing Discrepancies)
filtered_data <- full_data %>% 
  filter(!category %in% c("GrASE_Only", "rMATS_Only"))

stats_filtered <- calc_stats(filtered_data, p_col, eff_col, "FILTERED (Consistent)")
print("--- Summary Stats (Filtered: No Discrepancies) ---")
print(stats_filtered)

# --- Discrepancy Analysis ---
message("\n--- Discrepancy Analysis ---")

# 1. Summary of categories
message("Evaluation Categories (p < ", SIG_CUTOFF, "):")
print(table(full_data$category))

# 2. Check Direction agreement for Significant events (Both Sig)
both_sig <- full_data %>% filter(category == "Both_Sig")
message("\nDirection agreement in 'Both_Sig' events:")
print(table(both_sig$dir_agree))

# 3. Correlation of raw metrics if possible (Counts vs Means)
# GrASE usually has 'diff_mean' (counts for difference) and 'ref_mean' (counts for reference)
# rMATS has IJC_SAMPLE_1 etc in the raw file, but we only loaded PValue/FDR/IncLevelDiff. 
# We can't do exact count correlation without parsing raw columns, but we can look at effect sizes.

# 4. Save Discrepancies
# Export events where they disagree significantly
discrepancies <- full_data %>%
  filter(category %in% c("GrASE_Only", "rMATS_Only")) %>%
  select(gene, event, rMATS_EventType, ID, 
         grase_p, rmats_pval, 
         grase_eff = !!eff_col, rmats_inc = rmats_inc_diff, 
         category) %>%
  arrange(grase_p) # Sort by GrASE significance

discrep_file <- sub("\\.[^.]+$", ".discrepancies.txt", output_file)
write.table(discrepancies, discrep_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Discrepancy list saved to: ", discrep_file)
message("  (Contains events significant in one tool but not the other)")

# --- Plotting ---
message("\nGenerating plots to ", plot_file)
pdf(plot_file, width = 10, height = 5)

plot_dt <- full_data %>%
  mutate(
    log_grase = -log10(grase_p + 1e-50),
    log_rmats = -log10(rmats_pval + 1e-50),
    grase_eff = get(eff_col)
  )

# Plot 1: P-values with Category Coloring
# Filter removed to show capped values at 50
p1_dt <- plot_dt 
p1 <- ggplot(p1_dt, aes(x = log_grase, y = log_rmats, color = category)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = "P-value Comparison (Colored by Significance)",
       subtitle = paste("Significance Cutoff p <", SIG_CUTOFF, "- Total Dataset"),
       x = "GrASE -log10(p)", y = "rMATS -log10(p)")

# Plot 2: Effect Sizes Total
p2 <- ggplot(plot_dt, aes(x = grase_eff, y = rmats_inc_diff, color = dir_agree)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue", "NA" = "grey")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "Effect Size Comparison (Total)", 
       x = "GrASE Effect Size", y = "rMATS IncLevelDiff") +
  theme_minimal()

print(p1)
print(p2)

# --- Plots for Filtered Data ---
plot_filt_dt <- plot_dt %>% filter(!category %in% c("GrASE_Only", "rMATS_Only"))

p3 <- ggplot(plot_filt_dt, aes(x = log_grase, y = log_rmats, color = category)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(title = "P-value Comparison (Consistent Events Only)",
       subtitle = "Removing GrASE_Only / rMATS_Only discrepancies",
       x = "GrASE -log10(p)", y = "rMATS -log10(p)")

p4 <- ggplot(plot_filt_dt, aes(x = grase_eff, y = rmats_inc_diff, color = dir_agree)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "blue", "NA" = "grey")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(title = "Effect Size Comparison (Consistent Events Only)", 
       x = "GrASE Effect Size", y = "rMATS IncLevelDiff") +
  theme_minimal()

print(p3)
print(p4)

# Faceted versions if multiple types exist
if (length(unique(plot_dt$rMATS_EventType)) > 1) {
    print(p1 + facet_wrap(~rMATS_EventType) + ggtitle("Total Dataset by Type"))
    print(p3 + facet_wrap(~rMATS_EventType) + ggtitle("Consistent Events by Type"))
}

dev.off()

# --- Saving ---
write.table(full_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Results saved to: ", output_file)

