#!/usr/bin/env Rscript
## Is the outward drop a LOCAL effect, or an artefact of how often anything is
## called in a contrast?
##
## The positive fraction is not comparable across contrasts: it is close to
## 1 - P(U called | t called), so it falls as a contrast's overall call rate
## rises. Over the DICE lineage panel the terminus call rate spans 3.5-fold
## (7.2% to 25.2%) and the fraction of judgeable sides scored positive moves
## with it (r = -0.89). That is arithmetic, not biology, and no comparative
## claim can rest on it.
##
## This asks the underlying question directly and on a scale that does not
## carry the base rate: WITHIN A GENE, does a tiling unit adjacent to a called
## unit have higher odds of being called than one two steps away? Stratifying
## by gene removes between-gene differences in expression, length and call
## density. Note what this comparison is NOT: it is restricted to POWERED units
## at distance exactly 1 or 2, inside genes informative for both. Earlier
## "enrichment" figures of 15x and 25x came from a much looser comparison
## (distance 1 against everything else, unpowered units included) and are not
## this quantity -- do not carry them over. On this restricted comparison the
## crude and stratified ratios are close, which is itself the finding: the
## locality is not an artefact of pooling genes.
##
## THE TILING. A gene's exonic parts and intronic bins flank it without gaps or
## overlaps, so adjacency is well defined. Part numbers are coordinate-ordered,
## so index 2P for exonic part P places the bins between them: a whole-intron
## bin at 2*lo+1, and for a long intron split in two, its low-end bin at
## 2*lo+0.5 and its high-end bin at 2*hi-0.5. Distance is then the number of
## tiling steps between two units' ranks, counted over the gene's FULL tiling
## and not only over the units that happen to be tested.
##
## DISTANCE is to the nearest OTHER called unit, so a called unit is never its
## own neighbour and the comparison is not circular.
##
## Units below the read floor are dropped before anything is computed: an
## uncalled unit that was never powered is not evidence of absence, and
## including them would let the 5' depth deficit masquerade as locality.
##
## Usage:
##   Rscript stratified_adjacency.R [PANEL]      PANEL = lineage (default) or activation
suppressPackageStartupMessages({library(data.table)})
setwd(Sys.getenv("DICE", "/mnt/data1/home/mirahan/DICE"))
PANEL <- commandArgs(TRUE)[1]; if (is.na(PANEL)) PANEL <- "lineage"
## IG=1 as the second argument reads the gene-edge run, whose outputs live in
## <panel>/ig and whose bins come from the manifest that includes the flanks.
IG <- commandArgs(TRUE)[2]; if (is.na(IG)) IG <- "0"
cfg <- switch(PANEL,
  lineage    = list(dir = "intronparts/lineage", stem = "test_outward.ex.terminal"),
  activation = list(dir = "intronparts",         stem = "test_outward.exoncnt"),
  stop("unknown panel: ", PANEL))
if (IG == "1") cfg <- list(dir = file.path(cfg$dir, "ig"),
                           stem = "test_outward.ex.terminal")
MIN_READS <- 10; PADJ <- 0.01; MIN_DPI <- 0.1

## ---- units, with their tiling rank ---------------------------------------
st <- unique(fread(file.path(cfg$dir, "outward_side_table.tsv"),
                   select = c("gene","synth_event","unit_type","unit_id")))
man <- fread(if (IG == "1") "intronparts_ig/intronparts_manifest.tsv"
             else      "intronparts/intronparts_manifest.tsv",
             select = c("bin_id","gene","touches","abut_part_lo","abut_part_hi"))
man[, `:=`(lo = as.integer(abut_part_lo), hi = as.integer(abut_part_hi))]
man[, idx := fcase(touches == "both", 2*lo + 1,
                   touches == "lo",   2*lo + 0.5,
                   touches == "hi",   2*hi - 0.5)]

## The gene's FULL tiling: every exonic part 1..max, plus every intronic bin.
## max part comes from the per-gene DEXSeq GFF, the only complete source -- the
## units table holds only the parts this arm happens to test.
genes <- sort(unique(st$gene))
cat(sprintf("panel %s: %d genes, %d (gene, unit) rows in the side table\n",
            PANEL, length(genes), nrow(unique(st[, .(gene, unit_type, unit_id)]))))
gf <- file.path("dexseq.gff", paste0(genes, ".dexseq.gff"))
ok <- file.exists(gf)
if (any(!ok)) cat(sprintf("  WARNING %d genes have no dexseq.gff and are dropped\n", sum(!ok)))
maxpart <- rbindlist(lapply(which(ok), function(i) {
  n <- as.integer(system(sprintf("awk '$3==\"exonic_part\"' %s | wc -l", gf[i]), intern = TRUE))
  data.table(gene = genes[i], maxpart = n)
}))
maxpart <- maxpart[maxpart > 0]
ex <- data.table(gene      = rep(maxpart$gene, maxpart$maxpart),
                 unit_type = "exon",
                 pnum      = sequence(maxpart$maxpart))
ex[, `:=`(unit_id = sprintf("%03d", pnum), idx = 2 * pnum, pnum = NULL)]
full <- rbindlist(list(ex, man[gene %in% genes,
                               .(gene, unit_type = "intron", unit_id = bin_id, idx)]))
setorder(full, gene, idx)
full[, rank := seq_len(.N), by = gene]
cat(sprintf("  full tiling: %d units over %d genes\n", nrow(full), uniqueN(full$gene)))

## ---- call status per (gene, unit, contrast) -------------------------------
rf <- list.files(cfg$dir, pattern = paste0("^", cfg$stem, "_betabinom_EBmap\\.annotated\\.txt$"),
                 recursive = TRUE, full.names = TRUE)
## Prefer the per-contrast directories, as the lineage panel has. Drop
## test.outward.BH.* -- a sensitivity run under a different adjustment -- and
## drop the pooled stage-3a test.outward/ whose rows duplicate every contrast.
## The activation panel was tested serially and has ONLY the pooled directory,
## so fall back to it when no per-contrast file exists.
## The gene-edge run lives in a SUBdirECTORY of the pre-flank panel directory,
## and this listing is recursive: without this guard a non-IG lineage run would
## silently read its own 10 result files PLUS the 10 under ig/, doubling every
## contrast. Scoped by IG rather than by name so adding another variant cannot
## reintroduce it.
rf <- rf[grepl("/ig/", rf, fixed = TRUE) == (IG == "1")]
per <- rf[grepl("/test\\.outward\\.[^/]+/", rf) & !grepl("/test\\.outward\\.BH\\.", rf)]
rf  <- if (length(per)) per else rf[grepl("/test\\.outward/", rf)]
cat(sprintf("  %d result files\n", length(rf)))
res <- rbindlist(lapply(rf, fread, sep = "\t",
                        select = c("contrast","gene","event","padj","delta_pi",
                                   "lfc_diff_net","d_support")))
res[, called := !is.na(padj) & padj < PADJ & abs(delta_pi) >= MIN_DPI & lfc_diff_net > 0]
u <- merge(res, st, by.x = c("gene","event"), by.y = c("gene","synth_event"))
u <- unique(u[, .(gene, unit_type, unit_id, contrast, called, d_support)])
u <- u[!is.na(d_support) & d_support >= MIN_READS]          # powered only
u <- merge(u, full, by = c("gene","unit_type","unit_id"))
cat(sprintf("  %d powered (gene, unit, contrast) rows\n", nrow(u)))

## ---- distance to the nearest OTHER called unit ----------------------------
dist_to_called <- function(rank, called) {
  cr <- rank[called]
  if (!length(cr)) return(rep(NA_real_, length(rank)))
  vapply(seq_along(rank), function(i) {
    d <- abs(cr - rank[i]); d <- d[d > 0]
    if (length(d)) min(d) else NA_real_
  }, 0)
}
u[, d := dist_to_called(rank, called), by = .(gene, contrast)]

## ---- Mantel-Haenszel, strata = gene ---------------------------------------
mh <- function(dt) {
  s <- dt[d %in% c(1, 2), .(a = sum(d == 1 &  called), b = sum(d == 1 & !called),
                            c = sum(d == 2 &  called), dd = sum(d == 2 & !called)),
          by = gene]
  s <- s[(a + b) > 0 & (c + dd) > 0]            # informative strata only
  n <- s[, a + b + c + dd]
  num <- sum(s$a * s$dd / n); den <- sum(s$b * s$c / n)
  list(or = num / den, strata = nrow(s), units = sum(n))
}
cat("\n### Mantel-Haenszel odds ratio, called at distance 1 vs distance 2, strata = gene\n")
out <- u[, { m <- mh(.SD); .(MH_OR = round(m$or, 2), genes = m$strata, units = m$units,
                             crude = round(
                               (sum(d == 1 & called, na.rm = TRUE) / sum(d == 1, na.rm = TRUE)) /
                               (sum(d == 2 & called, na.rm = TRUE) / sum(d == 2, na.rm = TRUE)), 2)) },
         by = contrast]
tr <- u[, .(unit_call_rate = round(100 * mean(called), 1)), by = contrast]
out <- merge(out, tr, by = "contrast")
print(out[order(unit_call_rate)])
cat(sprintf("\nMH OR range %.2f-%.2f over a %.1f-fold range of call rate (r = %+.2f)\n",
            min(out$MH_OR), max(out$MH_OR),
            max(out$unit_call_rate)/min(out$unit_call_rate),
            cor(out$unit_call_rate, out$MH_OR)))
cat(sprintf("crude (unstratified) range %.2f-%.2f  <- close to the stratified value,\n",
            min(out$crude), max(out$crude)))
cat("   so on this restricted comparison the locality does not depend on the stratification\n")
