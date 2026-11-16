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
##       reference parts), bipartID (its event number) and inclSide (1 or 2, the
##       side whose path holds the inclusion form); NA when unmatched
##   combined.dexseq.<type>.mapped.txt  one row per exonic part of each gene:
##       the events whose matched bipartition has the part in a distinct set
##       (rMATS_ID_<type>) or in its reference (rMATS_ID_<type>_ref)
## compare_rmats_grase.R joins the first on (GeneID, bipartID) = (gene, event);
## compare_rmats_dexseq.R reads both.
##
## Several rMATS runs (e.g. the contrasts of one study, all on the same annotation) can be
## mapped in one call: --rmatsdir and --outdir take comma-separated lists of the same length.
## rMATS detects events separately per run, but the lists overlap heavily, and an event's
## match depends only on its coordinates and its gene. So the distinct events (same gene and
## coordinates) are mapped once and each run's output is written under its own event IDs,
## identical to mapping the runs one by one.
##
## Usage:
##   Rscript map_rmats_bipartitions.R --graphdir=<graphml dir> \
##     --splitdir=<dir of per-gene bipartition files> [--class=internal|terminal] \
##     --rmatsdir=<dir with fromGTF.<type>.txt>[,<dir>...] --outdir=<dir>[,<dir>...] [--mc_cores=8]
## --splitdir holds <gene>.bipartition.txt (as in a count directory) or, with
## --class, <gene>.bipartition.<class>.txt (as in bipartition.filtered). The
## event number is the file's `event` column when present, else the row number.
suppressPackageStartupMessages({library(optparse); library(data.table); library(igraph); library(parallel)})
o <- parse_args(OptionParser(option_list = list(
  make_option("--graphdir", type = "character"),
  make_option("--splitdir", type = "character"),
  make_option("--class",    type = "character", default = NULL),
  make_option("--rmatsdir", type = "character", help = "rMATS output dir(s), comma separated"),
  make_option("--outdir",   type = "character", help = "output dir(s), one per --rmatsdir"),
  make_option("--mc_cores", type = "integer", default = 8L))))
for (r in c("graphdir", "splitdir", "rmatsdir", "outdir"))
  if (is.null(o[[r]])) stop("--", r, " is required", call. = FALSE)
rdirs <- strsplit(o$rmatsdir, ",")[[1]]; odirs <- strsplit(o$outdir, ",")[[1]]
if (length(rdirs) != length(odirs)) stop("--rmatsdir and --outdir need the same number of directories", call. = FALSE)
for (d in odirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
suffix <- if (is.null(o$class)) ".bipartition.txt" else paste0(".bipartition.", o$class, ".txt")
etypes <- c("SE", "A3SS", "A5SS", "RI")

## events of each run; `key` = gene and coordinates (every column but ID and geneSymbol)
read_run <- function(dir) {
  ev <- list()
  for (et in etypes) {
    fp <- file.path(dir, sprintf("fromGTF.%s.txt", et))
    if (!file.exists(fp)) { message("no ", fp, ", skipping ", et); next }
    e <- fread(fp, quote = "\"", showProgress = FALSE)
    e[, `:=`(ID = as.character(ID), GeneID = gsub("\"", "", GeneID))]
    if ("geneSymbol" %in% names(e)) e[, geneSymbol := gsub("\"", "", geneSymbol)]
    kc <- setdiff(names(e), c("ID", "geneSymbol"))
    e[, key := do.call(paste, c(unname(as.list(e[, kc, with = FALSE])), sep = "\t"))]
    ev[[et]] <- e
  }
  ev
}
runs <- lapply(rdirs, read_run)
## distinct events over the runs, numbered uid
U <- list()
for (et in unique(unlist(lapply(runs, names)))) {
  u <- unique(rbindlist(lapply(runs, `[[`, et), use.names = TRUE), by = "key")
  u[, uid := seq_len(.N)]
  for (i in seq_along(runs)) if (!is.null(runs[[i]][[et]])) runs[[i]][[et]][, uid := u$uid[match(key, u$key)]]
  U[[et]] <- u
  if (length(runs) > 1)
    message(sprintf("%s: %d events in %d runs -> %d distinct", et,
                    sum(vapply(runs, function(r) NROW(r[[et]]), 0)), length(runs), nrow(u)))
}
Ug <- lapply(U, function(u) split(u, by = "GeneID"))
genes <- unique(unlist(lapply(U, `[[`, "GeneID")))

## match the distinct events of one gene
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
  res <- list(status = "ok", frags = frags)
  for (et in names(U)) {
    e <- Ug[[et]][[gene]]; if (is.null(e)) e <- U[[et]][0]
    m <- grase:::match_rmats_bipartitions(ge, bip, et, e)
    k <- !is.na(m)
    dist <- ifelse(k, paste(bip$setdiff1[m], bip$setdiff2[m], sep = ","), NA)
    dist <- vapply(dist, function(s) if (is.na(s)) NA_character_ else {
      p <- trimws(strsplit(s, ",")[[1]]); paste(unique(p[nzchar(p) & p != "NA"]), collapse = ",") }, "")
    res[[et]] <- data.table(uid = e$uid, dist = unname(dist), refp = ifelse(k, bip$ref_ex_part[m], NA),
                            bipartID = ifelse(k, evno[m], NA), inclSide = attr(m, "incl_side"))
  }
  res
}
## progress on stderr about every 5% of genes (genes are dispatched in order, one worker each, so the
## index of a finishing gene tracks the progress closely); the time per percent gives the remaining time
t0 <- Sys.time(); n_genes <- length(genes); step <- max(1L, n_genes %/% 20L)
message(sprintf("[%s] mapping %d genes on %d cores", format(t0, "%H:%M:%S"), n_genes, o$mc_cores))
res <- mclapply(seq_along(genes), function(i) {
  r <- tryCatch(one_gene(genes[i]), error = function(e) list(status = paste("error:", conditionMessage(e))))
  if (i %% step == 0L) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("[%s] %d / %d genes (%.0f%%), %.1f min elapsed, ~%.0f min left",
                    format(Sys.time(), "%H:%M:%S"), i, n_genes, 100 * i / n_genes, el, el * (n_genes - i) / i))
  }
  r
}, mc.cores = o$mc_cores, mc.preschedule = FALSE)
names(res) <- genes

## each run's output, under its own event IDs
for (ri in seq_along(runs)) {
  ev <- runs[[ri]]; od <- odirs[ri]
  if (length(runs) > 1) cat("==", rdirs[ri], "->", od, "\n")
  gr <- unique(unlist(lapply(ev, `[[`, "GeneID")))
  st <- vapply(res[gr], `[[`, "", "status")
  cat("genes:", length(gr), "|", paste(names(table(st)), table(st), sep = " ", collapse = ", "), "\n")
  if (any(startsWith(st, "error"))) cat("first errors:", head(unique(st[startsWith(st, "error")]), 3), sep = "\n  ")
  ok <- gr[st == "ok"]
  for (et in names(ev)) {
    Eg <- split(ev[[et]], by = "GeneID")
    per_gene <- mclapply(ok, function(gene) {
      r <- res[[gene]]; u <- r[[et]]; frags <- r$frags
      e <- Eg[[gene]]; if (is.null(e)) e <- ev[[et]][0]
      i <- match(e$uid, u$uid)
      dist <- u$dist[i]; refp <- u$refp[i]
      tag <- paste0(et, "_", e$ID)
      hit <- function(lst) vapply(frags, function(f) {
        w <- vapply(lst, function(s) !is.na(s) && f %in% strsplit(s, ",")[[1]], logical(1))
        if (any(w)) paste(tag[w], collapse = ",") else NA_character_ }, "")
      list(rows = data.table(ID = e$ID, DexseqFragment = dist, DexseqRefFrag = refp,
                             bipartID = u$bipartID[i], inclSide = u$inclSide[i]),
           mapped = setnames(data.table(GeneID = gene, DexseqFragment = frags,
                                        a = unname(hit(dist)), b = unname(hit(refp))),
                             c("a", "b"), paste0("rMATS_ID_", et, c("", "_ref"))))
    }, mc.cores = o$mc_cores)
    rows <- rbindlist(lapply(per_gene, `[[`, "rows"))
    out <- merge(ev[[et]][, !c("key", "uid")], rows, by = "ID", all.x = TRUE, sort = FALSE)
    fwrite(out, file.path(od, sprintf("combined.fromGTF.%s.txt", et)), sep = "\t", quote = FALSE, na = "NA")
    fwrite(rbindlist(lapply(per_gene, `[[`, "mapped")),
           file.path(od, sprintf("combined.dexseq.%s.mapped.txt", et)), sep = "\t", quote = FALSE, na = "NA")
    cat(sprintf("%-5s events %6d | matched %6d (%5.1f%%)\n", et, nrow(out),
                sum(!is.na(out$bipartID)), 100 * mean(!is.na(out$bipartID))))
  }
}
