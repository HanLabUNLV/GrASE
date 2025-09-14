#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(readr)
library(fs)

# =========================
# CONFIG — set your paths
# =========================
output_dir         <- "~/graphml.dexseq.v34/grase_results.integrate/results"
rmats_results_path <- file.path(output_dir, "SplicingEvents", "rMATS_TestedEvents.txt")
dexseq_tmp_dir     <- file.path(output_dir, "tmp")   # where combined.dexseq.*.mapped.txt live
merged_graph_path  <- "~/graphml.dexseq.v34/dice_exoncnts.filtered/test.txt"

summary_path       <- file.path(output_dir, "summary.new.txt")
intersections_path <- file.path(output_dir, "intersections.new.txt")

dir_create(output_dir)
dir_create(file.path(output_dir, "GraphEvents"))
dir_create(file.path(output_dir, "Crosswalks"))
dir_create(file.path(output_dir, "Merged"))

# =========================
# Helpers
# =========================

split_parts <- function(x, sep = ",", na_rm = TRUE, unique_out = FALSE, sort_out = FALSE) {
  if (is.null(x) || length(x) == 0) return(character())
  if (na_rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(character())

  y <- unlist(strsplit(as.character(x), sep, fixed = TRUE), use.names = FALSE)
  y <- trimws(y)
  y <- y[nzchar(y) & y != "NA"]

  if (unique_out) {
    y <- if (sort_out) sort(unique(y)) else unique(y)  # unique preserves first-seen order
  }
  y
}

canon_join <- function(x) paste(sort(unique(x)), collapse = ",")
set_equal_chr <- function(a, b) {
  a <- unique(a); b <- unique(b)
  length(a) == length(b) && all(a %in% b)
}
set_overlap_chr <- function(a, b) length(intersect(unique(a), unique(b))) > 0

append_rows <- function(rows_dt, path) {
  if (file.exists(path)) {
    old <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE))
    new <- bind_rows(old, as_tibble(rows_dt))
  } else {
    new <- as_tibble(rows_dt)
  }
  readr::write_tsv(new, path)
}

# =========================
# Readers
# =========================

# rMATS-DEX mapping files you already produced
read_rmats_dex_mapped <- function(tmp_dir) {
  types <- c("A3SS","A5SS","SE","RI")
  out <- rbindlist(lapply(types, function(ty) {
    fn <- file.path(tmp_dir, paste0('combined.fromGTF.', ty, '.txt')) 
    if (!file.exists(fn)) return(NULL)
    dt <- data.table(read.table(fn, sep = "\t", header = TRUE, na.strings = c("", "NA")))
    colnames(dt)[which(names(dt) == "ID")] <- "rMATSEventID"
    colnames(dt)[which(names(dt) == "DexseqFragment")] <- "rmats_alt_parts"
    colnames(dt)[which(names(dt) == "DexseqRefFrag")] <- "rmats_ref_parts"
    if (!"GeneID" %in% names(dt)) stop("Expected 'GeneID' in ", basename(fn))
    if (!"rMATSEventID" %in% names(dt)) stop("Expected rMATS ID column in ", basename(fn))
    dt$rMATSEventID = paste0(ty, "_", dt$rMATSEventID)
    # DexseqFragment (ALT) may be multiple rows per event; keep as-is
    dt[, EventType := ty]
    dt[, rmats_alt_parts := lapply(rmats_alt_parts,              function(z) split_parts(z, unique_out=TRUE, sort_out=TRUE))]
    dt[, rmats_ref_parts := lapply(rmats_ref_parts, function(z) split_parts(z, unique_out=TRUE, sort_out=TRUE))]
    dt[]
  }), fill = TRUE)
  if (is.null(out) || nrow(out) == 0) {
    stop("No rMATS-DEX mapping files found under: ", tmp_dir)
  }
  out
}

# Collapse to per-(GeneID, rMATSEventID) ALT/REF sets
rmats_parts_per_event <- function(map_dt) {
  map_dt[, .(
    rmats_alt_parts = list(unique(unlist(strsplit(paste(na.omit(DexseqFragment), collapse=","), ",")))),
    rmats_ref_parts = list(unique(unlist(strsplit(paste(na.omit(rMATSRefParts), collapse=","), ",")))),
    EventTypes      = paste(unique(na.omit(EventType)), collapse = ",")
  ), by = .(GeneID, rMATSEventID)]
}

# Graph (bipartition-exon) events
read_graph_events <- function(merged_graph_path) {
  g <- data.table(read.table(merged_graph_path, sep = "\t", header = TRUE, na.strings = c("", "NA")))
  needed <- c("gene","event","ref_ex_part","setdiff","LRT","pvalue","padj")
  miss <- setdiff(needed, names(g))
  if (length(miss)) stop("merged graph file missing columns: ", paste(miss, collapse=", "))
  g[, GeneID := as.character(gene)]
  g[, GraphEvent := as.integer(event)]
  g[, graph_ref_parts := lapply(ref_ex_part, function(z) split_parts(z, unique_out=TRUE, sort_out=TRUE))]
  g[, graph_alt_parts := lapply(setdiff, function(z) split_parts(z, unique_out=TRUE, sort_out=TRUE))]
  g[]
}

# rMATS per-event results (your rMATS_TestedEvents.txt)
read_rmats_results <- function(p) {
  stopifnot(file.exists(p))
  dt <- data.table(read.table(p, sep = "\t", header = TRUE, na.strings = c("", "NA")))
  if (!"GeneID" %in% names(dt)) stop("Expected column 'GeneID' in rMATS results.")
  if ("ID" %in% names(dt)) setnames(dt, "ID", "rMATSEventID")
  if (!"rMATSEventID" %in% names(dt)) stop("Expected column 'ID' (event id).")
  if ("FDR" %in% names(dt)) dt[, FDR := suppressWarnings(as.numeric(FDR))] else dt[, FDR := NA_real_]
  if ("padj" %in% names(dt)) dt[, rMATS_padj := suppressWarnings(as.numeric(padj))]
  if (!"EventType" %in% names(dt)) dt[, EventType := sub("_.*$", "", rMATSEventID)]
  dt[]
}

# =========================
# Crosswalk: Graph <-> rMATS
# =========================
match_graph_to_rmats_one_gene <- function(ge_gene, re_gene, mode = c("strict","lenient")) {
  mode <- match.arg(mode)
  if (nrow(ge_gene) == 0L || nrow(re_gene) == 0L) return(data.table())
  out <- vector("list", nrow(ge_gene))
  for (i in seq_len(nrow(ge_gene))) {
    g_ref <- ge_gene$graph_ref_parts[[i]]
    g_alt <- ge_gene$graph_alt_parts[[i]]
    hits <- logical(nrow(re_gene))
    for (j in seq_len(nrow(re_gene))) {
      r_ref <- re_gene$rmats_ref_parts[[j]]
      r_alt <- re_gene$rmats_alt_parts[[j]]
      hits[j] <- if (mode == "strict") {
        set_equal_chr(g_ref, r_ref) && set_equal_chr(g_alt, r_alt)
      } else {
        set_overlap_chr(g_ref, r_ref) && set_overlap_chr(g_alt, r_alt)
      }
    }
    if (any(hits)) {
      out[[i]] <- data.table(
        GeneID       = ge_gene$GeneID[i],
        GraphEvent   = ge_gene$GraphEvent[i],
        rMATSEventID = re_gene$rMATSEventID[hits],
        match_mode   = mode
      )
    } else {
      out[[i]] <- data.table(
        GeneID       = ge_gene$GeneID[i],
        GraphEvent   = ge_gene$GraphEvent[i],
        rMATSEventID = NA_character_,
        match_mode   = mode
      )
    }
  }
  rbindlist(out, fill = TRUE)
}

build_graph_rmats_crosswalk <- function(merged_graph_path, tmp_dir, mode = c("strict","lenient")) {
  mode <- match.arg(mode)
  map_dt <-  read_rmats_dex_mapped(tmp_dir)
  #re <- rmats_parts_per_event(map_dt)
  re <- map_dt
  ge <- read_graph_events(merged_graph_path)
  genes <- intersect(unique(ge$GeneID), unique(re$GeneID))
  res <- vector("list", length(genes))
  for (k in seq_along(genes)) {
    g <- genes[k]
    res[[k]] <- match_graph_to_rmats_one_gene(
      ge_gene = ge[GeneID == g, .(GeneID, GraphEvent, graph_ref_parts, graph_alt_parts)],
      re_gene = re[GeneID == g, .(rMATSEventID, rmats_ref_parts, rmats_alt_parts)],
      mode = mode
    )
  }
  unique(rbindlist(res, fill = TRUE))
}

# =========================
# MAIN
# =========================
#exit()

graph_all <- read_graph_events(merged_graph_path)
rmats_res <- read_rmats_results(rmats_results_path)

# Graph tested/sig
graph_tested <- graph_all[!is.na(p.value) | !is.na(padj)]
graph_sig    <- graph_all[!is.na(padj) & padj <= 0.05]
fwrite(graph_tested, file.path(output_dir, "GraphEvents", "Graph_TestedEvents.txt"), sep="\t")
fwrite(graph_sig,    file.path(output_dir, "GraphEvents", "Graph_SigEvents.txt"),    sep="\t")

xwalk_strict  <- build_graph_rmats_crosswalk(merged_graph_path, dexseq_tmp_dir, "strict")
xwalk_lenient <- build_graph_rmats_crosswalk(merged_graph_path, dexseq_tmp_dir, "lenient")
fwrite(xwalk_strict,  file.path(output_dir, "Crosswalks", "Graph_to_rMATS.strict.tsv"),  sep="\t")
fwrite(xwalk_lenient, file.path(output_dir, "Crosswalks", "Graph_to_rMATS.lenient.tsv"), sep="\t")

# Prefer 1:1 matches for final event-level merge
xwalk_1to1 <- xwalk_strict[!is.na(rMATSEventID)][, .N, by=.(GeneID, GraphEvent)][N==1][
  xwalk_strict, on=.(GeneID, GraphEvent)][, N := NULL][]

# Carry key graph annotations through
graph_keep <- graph_all[, .(gene, GeneID, event = GraphEvent, LRT, p.value, padj,
                            source, sink, ref_ex_part = graph_ref_parts, setdiff = graph_alt_parts)]

event_level_merged <- merge(
  graph_keep,
  xwalk_1to1, by.x = c("GeneID","event"), by.y = c("GeneID","GraphEvent"), all.x = TRUE
)

# Attach rMATS stats
rmats_keep <- rmats_res[, .(GeneID, rMATSEventID, rMATS_FDR = FDR)]
event_level_merged <- merge(event_level_merged, rmats_keep, by = c("GeneID","rMATSEventID"), all.x = TRUE)

# Canonicalize the set columns for readability
event_level_merged[, `:=`(
  ref_ex_part  = vapply(ref_ex_part, canon_join, "", USE.NAMES = FALSE),
  setdiff     = vapply(setdiff, canon_join, "", USE.NAMES = FALSE)
)]

fwrite(event_level_merged, file.path(output_dir, "Merged", "Graph_rMATS_event_level_merged.tsv"), sep="\t")

# Counts
num_graph_tested   <- nrow(graph_tested)
num_graph_sig      <- nrow(graph_sig)
num_xwalk_strict   <- sum(!is.na(xwalk_strict$rMATSEventID))
num_xwalk_lenient  <- sum(!is.na(xwalk_lenient$rMATSEventID))
num_rmats_tested   <- rmats_res %>% nrow()
num_rmats_sig      <- rmats_res %>% filter(!is.na(FDR) & FDR <= 0.05) %>% nrow()

# Event-level sig overlap via strict crosswalk
rmats_sig_ids <- rmats_res %>% filter(!is.na(FDR) & FDR <= 0.05) %>% pull(rMATSEventID) %>% unique()
graph_sig_ids <- event_level_merged %>%
  filter(!is.na(padj) & padj <= 0.05, !is.na(rMATSEventID)) %>%
  pull(rMATSEventID) %>% unique()
inter_event_sig_overlap <- length(intersect(graph_sig_ids, rmats_sig_ids))

summary_rows <- data.table(
  CountType = c(
    "Graph/FocalExon Tested Events",
    "Graph/FocalExon Sig Events (padj<=0.05)",
    "Graph:rMATS Crosswalk (strict) matches",
    "Graph:rMATS Crosswalk (lenient) matches",
    "rMATS Tested Events",
    "rMATS Sig Events (FDR<=0.05)",
    "Event-level Sig Overlap (via strict crosswalk)"
  ),
  Counts = c(
    num_graph_tested,
    num_graph_sig,
    num_xwalk_strict,
    num_xwalk_lenient,
    num_rmats_tested,
    num_rmats_sig,
    inter_event_sig_overlap
  )
)
append_rows(summary_rows, summary_path)

# Gene-level context
genes_graph_sig <- unique(graph_sig$GeneID)
genes_rmats_sig <- unique(rmats_res$GeneID[!is.na(rmats_res$FDR) & rmats_res$FDR <= 0.05])
inter_gene_sig_overlap <- length(intersect(genes_graph_sig, genes_rmats_sig))
intersections_rows <- data.table(
  Intersection = c(
    "Genes with Graph Sig (padj<=0.05)",
    "Genes with rMATS Sig (FDR<=0.05)",
    "Gene-level Sig Overlap"
  ),
  Counts = c(length(genes_graph_sig), length(genes_rmats_sig), inter_gene_sig_overlap)
)
append_rows(intersections_rows, intersections_path)


