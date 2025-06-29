# Description:
# It integrates and compares alternative splicing analysis results from rMATS and DEXSeq.
#
# Vocabulary:
#   rMATS - (fromGTF)
#   DEXSeq - (gff)
#   A3SS / A5SS - (overhang)
#   RI / SE - (full fragment)
#   graphml - (exon, intron, splicingGraphs)

# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(fs)

# Helper function to read and preprocess DEXSeq-to-rMATS mapping files
dex_to_mats <- function(file_path) {
  df <- readr::read_tsv(file_path, col_types = cols(.default = "c")) %>%
    filter(GeneID != "GeneID") %>%
    mutate(
      GeneID = str_trim(GeneID),
      DexseqFragment = str_trim(DexseqFragment)
    ) %>%
    arrange(GeneID, DexseqFragment)
  return(df)
}

# Helper function to read and preprocess rMATS-to-DEXSeq mapping files
mats_to_dex <- function(file_path, event_type) {
  df <- readr::read_tsv(file_path, col_types = cols(.default = "c")) %>%
    filter(GeneID != "GeneID") %>%
    mutate(ID = paste(event_type, ID, sep = "_")) %>%
    arrange(GeneID, ID)
  return(df)
}

# Function to get all necessary results files
get_results_files <- function(grase_directory, rmats_directory, dexseq_results) {
  results <- list()

  # Load rMATS results if specified
  rmats_dir <- fs::path_abs(rmats_directory)
  rmats_files <- fs::dir_ls(rmats_dir)

  JCEC <- readr::read_tsv(rmats_files[str_ends(rmats_files, "A3SS.MATS.JCEC.txt")], col_types = cols(.default = "c")) 
  colnames(JCEC)[1] <- "ID"
  results$A3SS_MATS <- JCEC %>% mutate(ID = paste("A3SS", ID, sep = "_"))
  JCEC <- readr::read_tsv(rmats_files[str_ends(rmats_files, "A5SS.MATS.JCEC.txt")], col_types = cols(.default = "c")) 
  colnames(JCEC)[1] <- "ID"
  results$A5SS_MATS <- JCEC %>% mutate(ID = paste("A5SS", ID, sep = "_"))
  JCEC <- readr::read_tsv(rmats_files[str_ends(rmats_files, "SE.MATS.JCEC.txt")], col_types = cols(.default = "c")) 
  colnames(JCEC)[1] <- "ID"
  results$SE_MATS <- JCEC %>% mutate(ID = paste("SE", ID, sep = "_"))
  JCEC <- readr::read_tsv(rmats_files[str_ends(rmats_files, "RI.MATS.JCEC.txt")], col_types = cols(.default = "c")) 
  colnames(JCEC)[1] <- "ID"
  results$RI_MATS <- JCEC %>% mutate(ID = paste("RI", ID, sep = "_"))
  
  results$output_dir <- fs::path(grase_directory, "results")
  grase_results_tmp <- fs::path(results$output_dir, "tmp")
  
  tmp_files <- fs::dir_ls(grase_results_tmp)

  # Load mapping files
  results$dex_to_A3SS <- dex_to_mats(tmp_files[str_ends(tmp_files, "A3SS.mapped.txt")])
  results$dex_to_A5SS <- dex_to_mats(tmp_files[str_ends(tmp_files, "A5SS.mapped.txt")])
  results$dex_to_SE <- dex_to_mats(tmp_files[str_ends(tmp_files, "SE.mapped.txt")])
  results$dex_to_RI <- dex_to_mats(tmp_files[str_ends(tmp_files, "RI.mapped.txt")])
  results$A3SS_to_dex <- mats_to_dex(tmp_files[str_ends(tmp_files, "fromGTF.A3SS.txt")], "A3SS")
  results$A5SS_to_dex <- mats_to_dex(tmp_files[str_ends(tmp_files, "fromGTF.A5SS.txt")], "A5SS")
  results$SE_to_dex <- mats_to_dex(tmp_files[str_ends(tmp_files, "fromGTF.SE.txt")], "SE")
  results$RI_to_dex <- mats_to_dex(tmp_files[str_ends(tmp_files, "fromGTF.RI.txt")], "RI")

  # Load DEXSeq results
  results$dexseqResults <- read.table(dexseq_results, sep="\t", header=TRUE)
    
  return(results)
}

# A safer, more idiomatic R function to filter data frames
filter_df_rMATS <- function(input_df, type, level, ...) {
  df <- input_df %>% filter(...)

  if (type == "Event") {
    if (level == "Secondary") {
      df <- df %>% 
        distinct(GeneID, ID, .keep_all = TRUE) %>%
        select(GeneID, ID, FDR, DexseqFragment, padj)
    }
    num_events <- nrow(df)
  }

  if (type == "Exon") {
    if (level == "Secondary") {
      df <- df %>% 
        distinct(groupID, featureID, .keep_all = TRUE) %>%
        select(groupID, featureID, padj, rMATS_ID, FDR)
    }
    num_events <- nrow(df)
  }

  return(list(df = df, num_events = num_events))
}


get_grase_results_rMATS <- function(grase_directory, rmats_directory, dexseq_results) {

  files <- get_results_files(grase_directory, rmats_directory, dexseq_results) 
  
  # Create output directories if they don't exist
  fs::dir_create(files$output_dir, "SplicingEvents", recurse = TRUE)
  fs::dir_create(files$output_dir, "ExonParts", recurse = TRUE)
  
  # Unpack lists for easier access
  list2env(files, envir = environment())

  # Exon Counts ###############################################################################
  dex_to_rmats <- full_join(dex_to_SE, dex_to_A5SS, by = c("GeneID", "DexseqFragment")) %>%
    full_join(., dex_to_A3SS, by = c("GeneID", "DexseqFragment")) %>%
    full_join(., dex_to_RI, by = c("GeneID", "DexseqFragment")) %>% 
    unite("rMATS_ID", c(rMATS_ID_A3SS, rMATS_ID_A5SS, rMATS_ID_SE, rMATS_ID_RI), sep = ",", na.rm = TRUE, remove = TRUE)

  colnames(dex_to_rmats)[1] <- 'groupID'
  colnames(dex_to_rmats)[2] <- 'featureID'
  dex_to_rmats_dexRes <- full_join(dexseqResults, dex_to_rmats, by = c("groupID", "featureID"))  %>%
    relocate(rMATS_ID, .after = featureID) %>%
    arrange(groupID, featureID)

  dex_to_rmats_exploded <- dex_to_rmats %>%
    separate_rows(rMATS_ID, sep = ",")

  dex_to_rmats_ex_MATS <- bind_rows(
      left_join(dex_to_rmats_exploded, A3SS_MATS, by = c("groupID" = "GeneID", "rMATS_ID" = "ID")),
      left_join(dex_to_rmats_exploded, A5SS_MATS, by = c("groupID" = "GeneID", "rMATS_ID" = "ID")),
      left_join(dex_to_rmats_exploded, SE_MATS, by = c("groupID" = "GeneID", "rMATS_ID" = "ID")),
      left_join(dex_to_rmats_exploded, RI_MATS, by = c("groupID" = "GeneID", "rMATS_ID" = "ID"))
    ) %>%
    filter(!is.na(geneSymbol)) # Corresponds to dropna(subset="ID.1")

  dex_to_rmats_ex_dexRes_MATS <- dex_to_rmats_dexRes %>%
    separate_rows(rMATS_ID, sep = ",") %>%
    left_join(dex_to_rmats_ex_MATS, by = c("groupID", "rMATS_ID", "featureID")) %>%
    mutate(across(c(FDR), as.numeric))
    
  num_exons_detected <- nrow(dex_to_rmats_dexRes)

  # rMATS Detected Exons
  res <- filter_df_rMATS(dex_to_rmats_ex_dexRes_MATS, "Exon", "Secondary", !is.na(rMATS_ID))
  exon_rmats_detected <- res$df; num_exons_rmats_detected <- res$num_events

  # rMATS Tested Exons
  res <- filter_df_rMATS(dex_to_rmats_ex_dexRes_MATS, "Exon", "Secondary", !is.na(ID))
  exon_rmats_tested <- res$df; num_exons_rmats_tested <- res$num_events
  
  # rMATS Significant Exons
  res <- filter_df_rMATS(dex_to_rmats_ex_dexRes_MATS, "Exon", "Secondary", FDR <= 0.05)
  exon_rmats_sig <- res$df; num_exons_rmats_sig <- res$num_events

  # DEXSeq Tested Exons
  res <- filter_df_rMATS(dex_to_rmats_ex_dexRes_MATS, "Exon", "Primary", !is.na(padj))
  exon_dex_tested_full <- res$df
  res_sec <- filter_df_rMATS(exon_dex_tested_full, "Exon", "Secondary", !is.na(rMATS_ID))
  exon_dex_tested_rmats_detected <- res_sec$df; num_exons_dex_tested_rmats_detected <- res_sec$num_events
  res_sec <- filter_df_rMATS(exon_dex_tested_full, "Exon", "Secondary", !is.na(ID))
  exon_dex_tested_rmats_tested <- res_sec$df; num_exons_dex_tested_rmats_tested <- res_sec$num_events
  res_sec <- filter_df_rMATS(exon_dex_tested_full, "Exon", "Secondary", FDR <= 0.05)
  exon_dex_tested_rmats_sig <- res_sec$df; num_exons_dex_tested_rmats_sig <- res_sec$num_events
  
  exon_dex_tested <- exon_dex_tested_full %>%
    distinct(groupID, featureID, .keep_all = TRUE) %>%
    select(groupID, featureID, padj, rMATS_ID, FDR)
  num_exons_dex_tested <- nrow(exon_dex_tested)

  # DEXSeq Significant Exons
  res <- filter_df_rMATS(dex_to_rmats_ex_dexRes_MATS, "Exon", "Primary", padj <= 0.05)
  exon_dex_sig_full <- res$df
  res_sec <- filter_df_rMATS(exon_dex_sig_full, "Exon", "Secondary", !is.na(rMATS_ID))
  exon_dex_sig_rmats_detected <- res_sec$df; num_exons_dex_sig_rmats_detected <- res_sec$num_events
  res_sec <- filter_df_rMATS(exon_dex_sig_full, "Exon", "Secondary", !is.na(ID))
  exon_dex_sig_rmats_tested <- res_sec$df; num_exons_dex_sig_rmats_tested <- res_sec$num_events
  res_sec <- filter_df_rMATS(exon_dex_sig_full, "Exon", "Secondary", FDR <= 0.05)
  exon_dex_sig_rmats_sig <- res_sec$df; num_exons_dex_sig_rmats_sig <- res_sec$num_events
  
  exon_dex_sig <- exon_dex_sig_full %>%
    distinct(groupID, featureID, .keep_all = TRUE) %>%
    select(groupID, featureID, padj, rMATS_ID, FDR)
  num_exons_dex_sig <- nrow(exon_dex_sig)


  # Event Counts ##################################################################################
  rmats_to_dex <- bind_rows(A3SS_to_dex, A5SS_to_dex, SE_to_dex, RI_to_dex) %>%
    arrange(GeneID, ID)

  rmats_to_dex_MATS <- bind_rows(
      left_join(rmats_to_dex, A3SS_MATS, by = c("GeneID", "ID")),
      left_join(rmats_to_dex, A5SS_MATS, by = c("GeneID", "ID")),
      left_join(rmats_to_dex, SE_MATS, by = c("GeneID", "ID")),
      left_join(rmats_to_dex, RI_MATS, by = c("GeneID", "ID"))
    ) %>%
    filter(!is.na(geneSymbol)) # Corresponds to dropna(subset="ID.1")
  names(rmats_to_dex_MATS)[names(rmats_to_dex_MATS)=='ID...2'] <- 'ID'

  rmats_to_dex_ex_dexRes <- rmats_to_dex %>%
    separate_rows(DexseqFragment, sep = ",") %>%
    left_join(rename(dexseqResults, GeneID = groupID, DexseqFragment = featureID), by = c("GeneID", "DexseqFragment")) %>%
    mutate(padj = as.numeric(padj))

  rmats_to_dex_ex_MATS_dexRes <- rmats_to_dex_MATS %>%
    separate_rows(DexseqFragment, sep = ",") %>%
    full_join(rmats_to_dex_ex_dexRes, by = c("GeneID", "ID", "DexseqFragment")) %>%
    mutate(across(c(padj, FDR), as.numeric))

  num_events_detected <- nrow(rmats_to_dex)
  
  # DEXSeq Tested Events
  res <- filter_df_rMATS(rmats_to_dex_ex_MATS_dexRes, "Event", "Secondary", !is.na(padj))
  event_dex_tested <- res$df; num_events_dex_tested <- res$num_events

  # DEXSeq Significant Events
  res <- filter_df_rMATS(rmats_to_dex_ex_MATS_dexRes, "Event", "Secondary", padj <= 0.05)
  event_dex_sig <- res$df; num_events_dex_sig <- res$num_events

  # rMATS Tested Events
  res <- filter_df_rMATS(rmats_to_dex_ex_MATS_dexRes, "Event", "Primary", !is.na(geneSymbol))
  event_rmats_tested_full <- res$df
  res_sec <- filter_df_rMATS(event_rmats_tested_full, "Event", "Secondary", !is.na(padj))
  event_rmats_tested_dex_tested <- res_sec$df; num_events_rmats_tested_dex_tested <- res_sec$num_events
  res_sec <- filter_df_rMATS(event_rmats_tested_full, "Event", "Secondary", padj <= 0.05)
  event_rmats_tested_dex_sig <- res_sec$df; num_events_rmats_tested_dex_sig <- res_sec$num_events

  event_rmats_tested <- event_rmats_tested_full %>%
    distinct(GeneID, ID, .keep_all = TRUE) %>%
    select(GeneID, ID, FDR, DexseqFragment, padj)
  num_events_rmats_tested <- nrow(event_rmats_tested)

  # rMATS Significant Events
  res <- filter_df_rMATS(rmats_to_dex_ex_MATS_dexRes, "Event", "Primary", FDR <= 0.05)
  event_rmats_sig_full <- res$df
  res_sec <- filter_df_rMATS(event_rmats_sig_full, "Event", "Secondary", !is.na(padj))
  event_rmats_sig_dex_tested <- res_sec$df; num_events_rmats_sig_dex_tested <- res_sec$num_events
  res_sec <- filter_df_rMATS(event_rmats_sig_full, "Event", "Secondary", padj <= 0.05)
  event_rmats_sig_dex_sig <- res_sec$df; num_events_rmats_sig_dex_sig <- res_sec$num_events

  event_rmats_sig <- event_rmats_sig_full %>%
    distinct(GeneID, ID, .keep_all = TRUE) %>%
    select(GeneID, ID, FDR, DexseqFragment, padj)
  num_events_rmats_sig <- nrow(event_rmats_sig)

  # Output Results #############################################################################
  summary_table <- tibble::tribble(
    ~CountType, ~Counts,
    "Total Exons Detected", num_exons_detected,
    "rMATS Detected Exons", num_exons_rmats_detected,
    "rMATS Tested Exons", num_exons_rmats_tested,
    "rMATS Sig Exons", num_exons_rmats_sig,
    "DEXSeq Tested Exons", num_exons_dex_tested,
    "DEXSeq Tested & rMATS Detected Exons", num_exons_dex_tested_rmats_detected,
    "DEXSeq Tested & rMATS Tested Exons", num_exons_dex_tested_rmats_tested,
    "DEXSeq Tested & rMATS Sig Exons", num_exons_dex_tested_rmats_sig,
    "DEXSeq Sig Exons", num_exons_dex_sig,
    "DEXSeq Sig & rMATS Detected Exons", num_exons_dex_sig_rmats_detected,
    "DEXSeq Sig & rMATS Tested Exons", num_exons_dex_sig_rmats_tested,
    "DEXSeq Sig & rMATS Sig Exons", num_exons_dex_sig_rmats_sig,
    "Total Events Detected", num_events_detected,
    "DEXSeq Tested Events", num_events_dex_tested,
    "DEXSeq Sig Events", num_events_dex_sig,
    "rMATS Tested Events", num_events_rmats_tested,
    "rMATS Tested & DEXSeq Tested Events", num_events_rmats_tested_dex_tested,
    "rMATS Tested & DEXSeq Sig Events", num_events_rmats_tested_dex_sig,
    "rMATS Sig Events", num_events_rmats_sig,
    "rMATS Sig & DEXSeq Tested Events", num_events_rmats_sig_dex_tested,
    "rMATS Sig & DEXSeq Sig Events", num_events_rmats_sig_dex_sig
  )
  
  readr::write_tsv(summary_table, fs::path(output_dir, "summary.txt"))

  # Note: The intersection calculations from the Python script are directly translated.
  # These assume certain set relationships that might require validation depending on the data.
  intersection_table <- tibble::tribble(
      ~Intersection, ~Counts,
      "DEXSeq Detected Exons Only", num_exons_detected - num_exons_rmats_detected - num_exons_dex_tested + num_exons_dex_tested_rmats_detected,
      "DEXSeq Detected & rMATS Detected Exons Only", num_exons_rmats_detected - num_exons_dex_tested_rmats_detected - num_exons_rmats_tested + num_exons_dex_tested_rmats_tested,
      "DEXSeq Detected & rMATS Tested Exons Only", num_exons_rmats_tested - num_exons_dex_tested_rmats_tested - num_exons_rmats_sig + num_exons_dex_tested_rmats_sig,
      "DEXSeq Detected & rMATS Sig Exons Only", num_exons_rmats_sig - num_exons_dex_tested_rmats_sig,
      "DEXSeq Tested Exons Only", num_exons_dex_tested - num_exons_dex_tested_rmats_detected - num_exons_dex_sig + num_exons_dex_sig_rmats_detected,
      "DEXSeq Tested & rMATS Detected Exons Only", num_exons_dex_tested_rmats_detected - num_exons_dex_tested_rmats_tested - num_exons_dex_sig_rmats_detected + num_exons_dex_sig_rmats_tested,
      "DEXSeq Tested & rMATS Tested Exons Only", num_exons_dex_tested_rmats_tested - num_exons_dex_tested_rmats_sig - num_exons_dex_sig_rmats_tested + num_exons_dex_sig_rmats_sig,
      "DEXSeq Tested & rMATS Sig Exons Only", num_exons_dex_tested_rmats_sig - num_exons_dex_sig_rmats_sig,
      "DEXSeq Sig Exons Only", num_exons_dex_sig - num_exons_dex_sig_rmats_detected,
      "DEXSeq Sig & rMATS Detected Exons Only", num_exons_dex_sig_rmats_detected - num_exons_dex_sig_rmats_tested,
      "DEXSeq Sig & rMATS Tested Exons Only", num_exons_dex_sig_rmats_tested - num_exons_dex_sig_rmats_sig,
      "DEXSeq Sig & rMATS Sig Exons", num_exons_dex_sig_rmats_sig,
      "rMATS Detected Events Only", num_events_detected - num_events_dex_tested - num_events_rmats_tested + num_events_rmats_tested_dex_tested,
      "rMATS Detected & DEXSeq Tested Events Only", num_events_dex_tested - num_events_dex_sig - num_events_rmats_tested_dex_tested + num_events_rmats_tested_dex_sig,
      "rMATS Detected & DEXSeq Sig Events Only", num_events_dex_sig - num_events_rmats_tested_dex_sig,
      "rMATS Tested Events Only", num_events_rmats_tested - num_events_rmats_tested_dex_tested - num_events_rmats_sig + num_events_rmats_sig_dex_tested,
      "rMATS Tested & DEXSeq Tested Events Only", num_events_rmats_tested_dex_tested - num_events_rmats_tested_dex_sig - num_events_rmats_sig_dex_tested + num_events_rmats_sig_dex_sig,
      "rMATS Tested & DEXSeq Sig Events Only", num_events_rmats_tested_dex_sig - num_events_rmats_sig_dex_sig,
      "rMATS Sig Events Only", num_events_rmats_sig - num_events_rmats_sig_dex_tested,
      "rMATS Sig & DEXSeq Tested Events Only", num_events_rmats_sig_dex_tested - num_events_rmats_sig_dex_sig,
      "rMATS Sig & DEXSeq Sig Events Only", num_events_rmats_sig_dex_sig
  )
  readr::write_tsv(intersection_table, fs::path(output_dir, "intersections.txt"))

  # Save event files
  readr::write_tsv(event_rmats_tested, fs::path(output_dir, "SplicingEvents", "rMATS_TestedEvents.txt"))
  readr::write_tsv(event_rmats_sig, fs::path(output_dir, "SplicingEvents", "rMATS_SigEvents.txt"))
  readr::write_tsv(event_dex_tested, fs::path(output_dir, "SplicingEvents", "DexTestedEvents.txt"))
  readr::write_tsv(event_dex_sig, fs::path(output_dir, "SplicingEvents", "DexSigEvents.txt"))
  readr::write_tsv(event_rmats_tested_dex_tested, fs::path(output_dir, "SplicingEvents", "rMATS_Tested__DexTestedEvents.txt"))
  readr::write_tsv(event_rmats_tested_dex_sig, fs::path(output_dir, "SplicingEvents", "rMATS_Tested__DexSigEvents.txt"))
  readr::write_tsv(event_rmats_sig_dex_tested, fs::path(output_dir, "SplicingEvents", "rMATS_Sig__DexTestedEvents.txt"))
  readr::write_tsv(event_rmats_sig_dex_sig, fs::path(output_dir, "SplicingEvents", "rMATS_Sig__DexSigEvents.txt"))

  # Save exon files
  readr::write_tsv(exon_dex_tested, fs::path(output_dir, "ExonParts", "DexTestedExons.txt"))
  readr::write_tsv(exon_dex_sig, fs::path(output_dir, "ExonParts", "DexSigExons.txt"))
  readr::write_tsv(exon_rmats_detected, fs::path(output_dir, "ExonParts", "rMATS_DetectedExons.txt"))
  readr::write_tsv(exon_rmats_tested, fs::path(output_dir, "ExonParts", "rMATS_TestedExons.txt"))
  readr::write_tsv(exon_rmats_sig, fs::path(output_dir, "ExonParts", "rMATS_SigExons.txt"))
  readr::write_tsv(exon_dex_tested_rmats_detected, fs::path(output_dir, "ExonParts", "DexTested__rMATS_DetectedExons.txt"))
  readr::write_tsv(exon_dex_tested_rmats_tested, fs::path(output_dir, "ExonParts", "DexTested__rMATS_TestedExons.txt"))
  readr::write_tsv(exon_dex_tested_rmats_sig, fs::path(output_dir, "ExonParts", "DexTested__rMATS_SigExons.txt"))
  readr::write_tsv(exon_dex_sig_rmats_detected, fs::path(output_dir, "ExonParts", "DexSig__rMATS_DetectedExons.txt")) # Mapped to correct df
  readr::write_tsv(exon_dex_sig_rmats_tested, fs::path(output_dir, "ExonParts", "DexSig__rMATS_TestedExons.txt"))  # Mapped to correct df
  readr::write_tsv(exon_dex_sig_rmats_sig, fs::path(output_dir, "ExonParts", "DexSig__rMATS_SigExons.txt"))    # Mapped to correct df
  
  # Save mapping files
  readr::write_tsv(dex_to_rmats_dexRes, fs::path(output_dir, "DEX_to_rMATS_Events.txt"))
  readr::write_tsv(rmats_to_dex_MATS, fs::path(output_dir, "rMATS_to_DEX_Exons.txt"))
  readr::write_tsv(dex_to_rmats_ex_dexRes_MATS, fs::path(output_dir, "Mapped.ExonsToEvents.txt"))
  readr::write_tsv(rmats_to_dex_ex_MATS_dexRes, fs::path(output_dir, "Mapped.EventsToExons.txt"))

  return(0)
}

grase_directory = "/Users/mirahan/Work/GrASE/grase_output_monaco_b_vs_cd8"
rmats_directory = "/Users/mirahan/Work/GrASE/rmats_output_monaco_b_vs_cd8"
dexseq_results = "/Users/mirahan/Work/GrASE/dexseq_output_monaco_b_vs_cd8/DEXSeq_Monaco_B_vs_CD8_results.txt"
nthread = 1 
                    

get_grase_results_rMATS(grase_directory, rmats_directory, dexseq_results)


