#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(grase)
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
              help = "Output file path [default= %default]", metavar = "character"),
  make_option("--contrast", type = "character", default = NULL,
              help = "GrASE contrast to compare; required when the GrASE results hold more than one"),
  make_option("--min_dpi", type = "double", default = 0.1,
              help = "GrASE |delta_pi| gate of is_significant() [default= %default]"),
  make_option("--min_dpi_sj", type = "double", default = 0.05,
              help = "the same for a junction-sourced side [default= %default]"),
  make_option("--min_dpsi", type = "double", default = 0.1,
              help = "rMATS |IncLevelDifference| a call needs [default= %default]"),
  make_option("--min_reads", type = "double", default = 10,
              help = "GrASE reads on the distinct set in a group [default= %default]"),
  make_option("--pval", type = "double", default = 0.05,
              help = "cutoff on the adjusted p of GrASE (is_significant) and rMATS FDR [default= %default]")
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

cts <- unique(grase_dt$contrast)
if (is.null(opt$contrast)) {
  if (length(cts) > 1) stop("GrASE results hold ", length(cts), " contrasts; pick one with --contrast: ", paste(cts, collapse = " "))
} else {
  grase_dt <- grase_dt[contrast == opt$contrast]
  if (!nrow(grase_dt)) stop("no rows for --contrast ", opt$contrast, "; the file has: ", paste(cts, collapse = " "))
}
# a GrASE side is called by grase::is_significant() at --pval: padj, lfc_diff_net, delta_pi gates, read support
grase_dt[, grase_call := is_significant(grase_dt, opt$pval, 0, opt$min_dpi, opt$min_dpi_sj, opt$min_reads)]
message("GrASE rows: ", nrow(grase_dt), " | called: ", sum(grase_dt$grase_call))

# --- Process Each rMATS Event Type ---
event_types <- c("SE", "MXE", "A3SS", "A5SS", "RI")
merged_list <- list()

for (etype in event_types) {
  # 1. Look for mapping file (now combined.fromGTF.*.txt)
  map_filename <- paste0("combined.fromGTF.", etype, ".txt")
  map_path <- file.path(map_dir, map_filename)
  
  # 2. Look for rMATS result file
  # Trying JC.txt (Junction Counts only) or JCEC.txt (Junction Counts + Exon Counts)
  rmats_filename <- paste0(etype, ".MATS.JCEC.txt")
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

# Create classification BEFORE correlations so we can filter.
# Calls: GrASE is_significant() at --pval, rMATS FDR < --pval with |IncLevelDifference| >= --min_dpsi; correlations and plots use the unadjusted ones.
SIG_CUTOFF <- opt$pval
full_data <- full_data %>%
  mutate(
    grase_p = get(p_col),
    grase_sig = grase_call,
    rmats_sig = !is.na(rmats_fdr) & rmats_fdr < SIG_CUTOFF & !is.na(rmats_inc_diff) & abs(rmats_inc_diff) >= opt$min_dpsi,
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
# Coverage: the table above holds only events present on both sides (a bipartition matched to an
# rMATS event); count what each tool called outside it.
rm_all <- rbindlist(lapply(event_types, function(et) {
  f <- file.path(rmats_dir, paste0(et, ".MATS.JCEC.txt"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, select = c("ID", "FDR", "IncLevelDifference")); x[, tag := paste0(et, "_", ID)]
  x[, call := !is.na(FDR) & FDR < opt$pval & abs(IncLevelDifference) >= opt$min_dpsi]; x
}))
matched_tag <- unique(paste0(full_data$rMATS_EventType, "_", full_data$ID))
gb <- unique(grase_dt[!is.na(p.value), .(gene, event, call = grase_call)])[, .(call = any(call)), by = .(gene, event)]
mb <- unique(full_data[, .(gene, event)])
message(sprintf("Events | GrASE bipartitions tested %d, called %d; matched to an rMATS event %d (called %d)",
                nrow(gb), sum(gb$call), nrow(merge(gb, mb)), sum(merge(gb, mb)$call)))
message(sprintf("Events | rMATS events tested %d, called %d; matched to a bipartition %d (called %d)",
                sum(!is.na(rm_all$FDR)), sum(rm_all$call, na.rm = TRUE),
                sum(rm_all$tag %in% matched_tag), sum(rm_all$tag %in% matched_tag & rm_all$call, na.rm = TRUE)))

message("\n--- Discrepancy Analysis ---")

# 1. Summary of categories
message("Evaluation Categories (padj / FDR < ", SIG_CUTOFF, "):")
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
         grase_p, grase_padj = padj, rmats_pval, rmats_fdr, 
         grase_eff = !!eff_col, rmats_inc = rmats_inc_diff, 
         category) %>%
  arrange(grase_padj) # Sort by GrASE significance

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
       subtitle = paste("Colour: padj / FDR <", SIG_CUTOFF, "- Total Dataset"),
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



# --- Exonic parts ---------------------------------------------------------------------------
# The same comparison on exonic parts, with no restriction to events both tools matched. GrASE
# parts: both distinct sets (setdiff1 and setdiff2) of each bipartition, called if either side is

# called if any event holding the part has FDR < --pval. A part is tested when its tool gave it a p-value.
message("\n--- Exonic parts ---")
gp <- rbind(grase_dt[!is.na(p.value), .(gene, event, call = grase_call, part = setdiff1)],
            grase_dt[!is.na(p.value), .(gene, event, call = grase_call, part = setdiff2)])
gp <- gp[!is.na(part) & part != "NA" & nzchar(part)]
gp <- gp[, .(part = trimws(unlist(strsplit(part, ",")))), by = .(gene, event, call)]
gp <- gp[, .(grase_called = any(call)), by = .(gene, part)]
rp <- rbindlist(lapply(event_types, function(et) {
  f <- file.path(map_dir, sprintf("combined.dexseq.%s.mapped.txt", et))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, colClasses = "character"); setnames(x, 1:3, c("gene", "part", "ids"))
  x <- x[!is.na(ids) & ids != "NA", .(tag = unlist(strsplit(ids, ","))), by = .(gene = gsub('"', "", gene), part)]
  merge(x, rm_all[, .(tag, FDR, call)], by = "tag")[!is.na(FDR)]
}))
rp <- rp[, .(rmats_called = any(call)), by = .(gene, part)]
parts <- merge(gp, rp, by = c("gene", "part"), all = TRUE)
parts[, `:=`(grase_tested = !is.na(grase_called), rmats_tested = !is.na(rmats_called))]
parts[is.na(grase_called), grase_called := FALSE]; parts[is.na(rmats_called), rmats_called := FALSE]
parts[, category := fifelse(grase_called & rmats_called, "Both_Sig", fifelse(grase_called, "GrASE_Only",
                    fifelse(rmats_called, "rMATS_Only", "Neither_Sig")))]
message(sprintf("Parts | GrASE tested %d, called %d | rMATS tested %d, called %d | tested by both %d",
                sum(parts$grase_tested), sum(parts$grase_called), sum(parts$rmats_tested), sum(parts$rmats_called),
                sum(parts$grase_tested & parts$rmats_tested)))
print(table(parts$category))
parts_file <- sub("\\.[^.]+$", ".parts.txt", output_file)
if (parts_file == output_file) parts_file <- paste0(output_file, ".parts.txt")
fwrite(parts, parts_file, sep = "\t")
message("Part-level table saved to: ", parts_file)
