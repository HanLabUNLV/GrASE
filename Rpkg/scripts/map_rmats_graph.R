#!/usr/bin/env Rscript

# scripts/map_rmats_graph.R
#
# Graph-based mapping of rMATS events to DEXSeq exonic parts, without
# bipartition splits.
#
# For each rMATS event the splice graph defines a local bubble with two
# boundary vertices (divergence point and reconvergence point).  Only
# transcripts incident to BOTH boundary vertices are local to the event --
# transcripts that do not pass through that region at all are excluded.
# Local transcripts are then split into path1 (inclusion/long-exon form,
# those that use exonic edges within the alternative region) and path2
# (exclusion/short-exon form, remaining local transcripts).  The setdiff
# of their edge sets -- edges where any path1 transcript is TRUE and all
# path2 transcripts are FALSE -- gives the DEXSeq fragments specific to the
# alternative isoform.
#
# Key graph structure note:
#   The graphml has two classes of edges:
#     ex/in   Regular splice-graph edges: have from_pos/to_pos attributes,
#             per-transcript booleans, but empty dexseq_fragment.
#     ex_part DEXSeq exonic-part edges: have dexseq_fragment values and
#             per-transcript booleans, but from_pos/to_pos stored as "NA".
#             Their positions must be read via their source/target vertex positions.
#
# Usage:
#   Rscript map_rmats_graph.R \
#     <graphml_dir>  directory with per-gene .graphml files
#     <rmats_dir>    directory with fromGTF.*.txt (and *.MATS.JCEC.txt)
#     <out_dir>      output directory for combined.fromGTF.*.txt files
#
# Output format identical to map_rmats_coords.R (ID, GeneID, DexseqFragment).

suppressPackageStartupMessages({
  library(igraph)
  library(data.table)
  library(grase)
})

# Source shared graph utilities from the GrASE package R file
.script_dir <- dirname(normalizePath(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]),
  mustWork = FALSE
))


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat(paste(
    "Usage: Rscript map_rmats_graph.R",
    "  <graphml_dir>  directory with per-gene .graphml files",
    "  <rmats_dir>    directory with fromGTF.*.txt files",
    "  <out_dir>      output directory",
    sep = "\n"
  ), "\n")
  quit(save = "no", status = 1)
}

graphml_dir <- args[1]
rmats_dir   <- args[2]
out_dir     <- args[3]
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# per-gene and per-event mapping 
# Thin wrappers around shared functions from Rpkg/R/rMATS.R

map_event <- function(ge, etype, strand, row) {
  reg <- rmats_event_boundaries(etype, strand, row)
  if (is.null(reg)) return(character(0))
  graph_setdiff_frags(ge, reg$b1, reg$b2, reg$alt_low, reg$alt_hi)
}

# load gene list 

event_types <- c("SE", "A3SS", "A5SS", "RI")

cat("Scanning rMATS fromGTF files for gene list...\n")
all_genes <- character(0)
for (etype in event_types) {
  fp <- file.path(rmats_dir, sprintf("fromGTF.%s.txt", etype))
  if (!file.exists(fp)) next
  dt <- fread(fp, select = "GeneID", quote = "\"")
  dt[, GeneID := gsub("\"", "", GeneID)]
  all_genes <- union(all_genes, unique(dt$GeneID))
}
cat(sprintf("  %d unique genes across all event types\n", length(all_genes)))

# load and precompute graphs 

cat("Loading and precomputing graphs...\n")
gene_data <- list()
n_loaded  <- 0L
for (gene in all_genes) {
  gpath <- file.path(graphml_dir, paste0(gene, ".graphml"))
  if (!file.exists(gpath)) next
  g <- tryCatch(
    read_graph(gpath, format = "graphml"),
    error = function(e) NULL
  )
  if (is.null(g)) next
  gene_data[[gene]] <- precompute_gene_graph(g)
  n_loaded <- n_loaded + 1L
}
cat(sprintf("  Precomputed %d genes (%.1f%%)\n",
            n_loaded, 100 * n_loaded / length(all_genes)))

# map each event type 

for (etype in event_types) {
  cat(sprintf("\nMapping %s events...\n", etype))

  fp <- file.path(rmats_dir, sprintf("fromGTF.%s.txt", etype))
  if (!file.exists(fp)) {
    message("  Not found, skipping")
    next
  }

  ev <- fread(fp, quote = "\"")
  ev[, GeneID := gsub("\"", "", GeneID)]
  if ("geneSymbol" %in% names(ev)) {
    ev[, geneSymbol := gsub("\"", "", geneSymbol)]
  }
  cat(sprintf("  %d events across %d genes\n", nrow(ev), uniqueN(ev$GeneID)))

  genes_ok <- intersect(ev$GeneID, names(gene_data))
  cat(sprintf("  %d genes have graph (%.1f%%)\n",
              length(genes_ok), 100 * length(genes_ok) / uniqueN(ev$GeneID)))

  result_rows <- vector("list", nrow(ev))

  for (gene in genes_ok) {
    ge   <- gene_data[[gene]]
    ev_g <- ev[GeneID == gene]

    for (i in seq_len(nrow(ev_g))) {
      row   <- ev_g[i]
      frags <- map_event(ge, etype, row$strand, row)
      if (length(frags) == 0L) next
      idx <- which(ev$ID == row$ID)[1L]
      result_rows[[idx]] <- data.table(
        ID             = as.character(row$ID),
        GeneID         = gene,
        DexseqFragment = paste(frags, collapse = ",")
      )
    }
  }

  keep   <- !vapply(result_rows,
                    function(x) is.null(x) || nrow(x) == 0L,
                    logical(1))
  out_dt <- rbindlist(result_rows[keep])

  n_with_fdr <- tryCatch({
    jf <- file.path(rmats_dir, sprintf("%s.MATS.JCEC.txt", etype))
    if (file.exists(jf)) {
      j <- fread(jf, select = "ID", quote = "\"")
      length(intersect(out_dt$ID, j$ID))
    } else NA_integer_
  }, error = function(e) NA_integer_)

  cat(sprintf("  Mapped: %d / %d events (%.1f%%)\n",
              nrow(out_dt), nrow(ev), 100 * nrow(out_dt) / nrow(ev)))
  if (!is.na(n_with_fdr)) {
    cat(sprintf("  Of those, %d also have MATS.JCEC FDR values\n",
                n_with_fdr))
  }

  out_file <- file.path(out_dir, sprintf("combined.fromGTF.%s.txt", etype))
  fwrite(out_dt, out_file, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("  Written -> %s\n", out_file))
}

cat("\nDone.\n")
