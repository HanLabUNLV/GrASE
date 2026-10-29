#!/usr/bin/env Rscript
## Map rMATS events onto GrASE bipartitions through the splicing graph.
##
## Each rMATS event has two forms (inclusion and skipping; long and short;
## retained and spliced). Both are placed on the splicing graph and the event
## is matched to the bipartition whose two sides hold the two forms; see
## match_rmats_bipartitions() in R/rMATS.R for the rule and tie-breaking.
##
## Writes, per event type, into --outdir:
##   combined.fromGTF.<type>.txt        every fromGTF row, plus DexseqFragment
##       (the matched bipartition's distinct parts), DexseqRefFrag (its
##       reference parts) and bipartID (its event number); NA when unmatched
##   combined.dexseq.<type>.mapped.txt  one row per exonic part of each gene:
##       the events whose matched bipartition has the part in a distinct set
##       (rMATS_ID_<type>) or in its reference (rMATS_ID_<type>_ref)
## compare_rmats_grase.R joins the first on (GeneID, bipartID) = (gene, event);
## compare_rmats_dexseq.R reads both.
##
## Usage:
##   Rscript map_rmats_bipartitions.R --graphdir=<graphml dir> \
##     --splitdir=<dir of per-gene bipartition files> [--class=internal|terminal] \
##     --rmatsdir=<dir with fromGTF.<type>.txt> --outdir=<dir> [--mc_cores=8]
## --splitdir holds <gene>.bipartition.txt (as in a count directory) or, with
## --class, <gene>.bipartition.<class>.txt (as in bipartition.filtered). The
## event number is the file's `event` column when present, else the row number.
suppressPackageStartupMessages({library(optparse); library(data.table); library(igraph); library(parallel)})
o <- parse_args(OptionParser(option_list = list(
  make_option("--graphdir", type = "character"),
  make_option("--splitdir", type = "character"),
  make_option("--class",    type = "character", default = NULL),
  make_option("--rmatsdir", type = "character"),
  make_option("--outdir",   type = "character"),
  make_option("--mc_cores", type = "integer", default = 8L))))
for (r in c("graphdir", "splitdir", "rmatsdir", "outdir"))
  if (is.null(o[[r]])) stop("--", r, " is required", call. = FALSE)
dir.create(o$outdir, showWarnings = FALSE, recursive = TRUE)
suffix <- if (is.null(o$class)) ".bipartition.txt" else paste0(".bipartition.", o$class, ".txt")
etypes <- c("SE", "A3SS", "A5SS", "RI")

ev <- list()
for (et in etypes) {
  fp <- file.path(o$rmatsdir, sprintf("fromGTF.%s.txt", et))
  if (!file.exists(fp)) { message("no ", fp, ", skipping ", et); next }
  ev[[et]] <- fread(fp, quote = "\"", showProgress = FALSE)
  ev[[et]][, `:=`(ID = as.character(ID), GeneID = gsub("\"", "", GeneID))]
  if ("geneSymbol" %in% names(ev[[et]])) ev[[et]][, geneSymbol := gsub("\"", "", geneSymbol)]
}
genes <- unique(unlist(lapply(ev, `[[`, "GeneID")))

one_gene <- function(gene) {
  gf <- file.path(o$graphdir, paste0(gene, ".graphml"))
  bf <- file.path(o$splitdir, paste0(gene, suffix))
  if (!file.exists(gf)) return(list(status = "no_graph"))
  if (!file.exists(bf)) return(list(status = "no_bipartitions"))
  bip <- tryCatch(read.table(bf, header = TRUE, sep = "\t", colClasses = "character",
                             quote = "", comment.char = ""), error = function(e) data.frame())
  if (!nrow(bip)) return(list(status = "no_bipartitions"))
  evno <- if ("event" %in% names(bip)) bip$event else as.character(seq_len(nrow(bip)))
  ge <- grase::precompute_gene_graph(read_graph(gf, format = "graphml"))
  frags <- ge$dex_frag[!is.na(ge$dex_frag) & nzchar(ge$dex_frag)]
  frags <- sort(unique(sprintf("E%03d", as.integer(frags))))
  res <- list(status = "ok")
  for (et in names(ev)) {
    e <- ev[[et]][GeneID == gene]
    m <- grase:::match_rmats_bipartitions(ge, bip, et, e)
    k <- !is.na(m)
    tag <- paste0(et, "_", e$ID)
    dist <- ifelse(k, paste(bip$setdiff1[m], bip$setdiff2[m], sep = ","), NA)
    dist <- vapply(dist, function(s) if (is.na(s)) NA_character_ else {
      p <- trimws(strsplit(s, ",")[[1]]); paste(unique(p[nzchar(p) & p != "NA"]), collapse = ",") }, "")
    refp <- ifelse(k, bip$ref_ex_part[m], NA)
    hit <- function(lst) vapply(frags, function(f) {
      w <- vapply(lst, function(s) !is.na(s) && f %in% strsplit(s, ",")[[1]], logical(1))
      if (any(w)) paste(tag[w], collapse = ",") else NA_character_ }, "")
    res[[et]] <- list(
      rows = data.table(ID = e$ID, DexseqFragment = unname(dist), DexseqRefFrag = refp,
                        bipartID = ifelse(k, evno[m], NA)),
      mapped = setnames(data.table(GeneID = gene, DexseqFragment = frags,
                                   a = unname(hit(dist)), b = unname(hit(refp))),
                        c("a", "b"), paste0("rMATS_ID_", et, c("", "_ref"))))
  }
  res
}
res <- mclapply(genes, function(g) tryCatch(one_gene(g),
                error = function(e) list(status = paste("error:", conditionMessage(e)))),
                mc.cores = o$mc_cores, mc.preschedule = FALSE)
names(res) <- genes
st <- vapply(res, `[[`, "", "status")
cat("genes:", length(genes), "|", paste(names(table(st)), table(st), sep = " ", collapse = ", "), "\n")
if (any(startsWith(st, "error"))) cat("first errors:", head(unique(st[startsWith(st, "error")]), 3), sep = "\n  ")

for (et in names(ev)) {
  rows <- rbindlist(lapply(res, function(r) r[[et]]$rows))
  out <- merge(ev[[et]], rows, by = "ID", all.x = TRUE, sort = FALSE)
  fwrite(out, file.path(o$outdir, sprintf("combined.fromGTF.%s.txt", et)), sep = "\t", quote = FALSE, na = "NA")
  fwrite(rbindlist(lapply(res, function(r) r[[et]]$mapped)),
         file.path(o$outdir, sprintf("combined.dexseq.%s.mapped.txt", et)), sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("%-5s events %6d | matched %6d (%5.1f%%)\n", et, nrow(out),
              sum(!is.na(out$bipartID)), 100 * mean(!is.na(out$bipartID))))
}
