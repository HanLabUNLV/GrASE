#!/usr/bin/env Rscript
## Enumerate the OUTWARD test units for the TERMINAL-event boundary arm.
##
## WHAT THE ARM CLAIMS. A bipartition side asserts that its transcription starts
## (TSS) or ends (TTS) differ in usage from the other side's. A side has one
## start per ROUTE, not one per side -- 57.7% of sides have more than one, mean
## 2.53, max 14 -- so the claim is tested once per terminus, plus once outward of
## D's own outer edge.
##
## WHY NOT A SINGLE BOUNDARY. Taking one position per side -- which is what the
## removed edge-coverage and flank-window tests did -- lands on the OUTERMOST
## terminus only 50.8% of the time and on some middle terminus 48.5% of the
## time. So every boundary result from those tests interrogated an effectively
## arbitrary one of the side's starts. Enumerating them removes the choice.
##
## TERMINI COME FROM EDGES, NOT VERTICES. The terminal exonic part of a route is
## the `ex_part` edge SOURCED at the route's first vertex after R (sunk at the
## last vertex before L, for TTS). That yields the part ID directly, which is
## already the tiling's key -- no coordinates, no strand arithmetic. Do NOT scan
## the path for its first ex_part edge: routes traverse `ex` edges (e.g. 1->6)
## that span several parts and are not quantification units, so a scan returns
## the wrong part. There are 35 of each in PDCD4 and only 13 shared vertex pairs.
##
## NO ref_ex_part GUARD IS NEEDED HERE. A TSS bubble's paths reconverge at the
## sink, so shared sequence lies INWARD; stepping outward cannot reach it.
## Verified: 0 of 10,573 outward units are the side's reference, against 20.7%
## of inward units. The inward arm is separate and needs that guard.
##
## THE TERMINUS PART ITSELF IS ALSO A TEST UNIT (anchor_kind = "terminus_self").
## Testing only the flank OUTWARD of each start answers "does the effect continue
## past it" but never "is this start used at all". With both, the walk along a
## side alternates
##     U_n  t_n  U_n-1  t_n-1  ...  U_1  t_1  [D]
## from outermost inward, and the real start is the outermost t_i that is called
## with everything outward of it uncalled -- a localization rather than a
## bounded/unbounded flag. It also rescues sides whose outermost flank is below
## the read floor: if t_n itself is powered and uncalled, that start is unused
## regardless of what lies beyond it. 72.4% of terminus parts are not reachable
## as any anchor's outward flank, so this is not redundant; all are exonic, so
## their counts come from the existing per-part DEXSeq files with no new
## counting pass.
##
## A terminus part that lies INSIDE D is still emitted. The test is then a
## sub-test of the distinct set the exonic arm already called, which resolves
## which of the side's starts carries the signal -- but it is NOT independent
## evidence and will correlate with D's own call by construction. The reference
## stays the side's own ref_ex_part in every case.
##
## ROWS ARE KEYED BY UNIT, NOT BY TERMINUS. A side's termini are often
## consecutive parts (PDCD4 event 11 side 1: 001, 002, 003, 004), so the flanks
## outward of them collapse onto one another and onto the termini themselves.
## Emitting per terminus would test the same unit repeatedly under different
## labels. The test unit is (gene, unit, ref_ex_part) -- reference counts are
## just the counts of ref_ex_part, so sides sharing a reference part share a
## test. The side table carries the one-to-many mapping back.
##
## Usage:
##   Rscript scripts/enumerate_outward_units.R --results <annotated.txt> \
##     --split_dir <bipartition.filtered> --graphml_dir <graphml.v34> \
##     --intron_manifest <intronparts_manifest.tsv> \
##     --out <units.tsv> [--mc_cores 8]
suppressPackageStartupMessages({library(data.table); library(igraph); library(optparse)
                                library(grase)})

o <- parse_args(OptionParser(option_list = list(
  make_option("--results",         type = "character"),
  make_option("--split_dir",       type = "character"),
  make_option("--graphml_dir",     type = "character"),
  make_option("--intron_manifest", type = "character"),
  make_option("--out",             type = "character"),
  make_option("--mc_cores",        type = "integer", default = 8L),
  make_option("--all_sides",       action = "store_true", default = FALSE,
              help = "enumerate EVERY tested terminal side, not only the called
                      ones. For building a matched control set: an uncalled
                      side's termini are the same geometry, and the long-read
                      check is positional, so no counting is needed. Writes a
                      `called` column so the two can be told apart. Do not use
                      this for the production units table -- build_outward_counts.R
                      would then count every side in the annotation."))))
for (r in c("results","split_dir","graphml_dir","intron_manifest","out"))
  if (is.null(o[[r]])) stop("--", r, " is required")
pid <- function(x) sprintf("%03d", as.integer(x))

res <- fread(o$results, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene","event")))
res[, side := fifelse(comparison == "diff1_vs_ref", 1L, 2L)]
## ref_ex_part is per (gene, event, side) and does not depend on the call, so
## taking it for every tested side is the same lookup over a bigger set.
sig <- unique(res[significant == TRUE, .(gene, event, side)])
sig[, called := TRUE]
if (o$all_sides) {
  called <- unique(res[, .(gene, event, side, ref_ex_part)])
  called <- merge(called, sig, by = c("gene", "event", "side"), all.x = TRUE)
  called[is.na(called), called := FALSE]
} else {
  called <- unique(res[significant == TRUE, .(gene, event, side, ref_ex_part)])
  called[, called := TRUE]
}
setkey(called, gene, event, side)
cat(sprintf("%s sides: %d over %d genes (called %d)\n",
            if (o$all_sides) "all tested" else "called",
            nrow(called), uniqueN(called$gene), sum(called$called)))

## GEOMETRY COMES FROM THE GRAPH, not from a GFF-derived table. This script
## already reads each gene's graphml to resolve a route's terminal vertex to a
## part, so taking coordinates from a separate file meant two sources for one
## fact. The old --exon_parts pointed at tss_edge/geom_parts.tsv, a by-product
## of the GFF-based edge test, which is also why this arm still reached into
## that deleted test's output directory.
##
## An ex_part edge spans its two vertices' positions HALF-OPEN, [lo, hi), while
## the GFF writes a closed [start, end]; hence end = hi - 1L. Verified against
## geom_parts.tsv on 1,872 parts over 60 genes: start == lo and end == hi - 1 on
## every one. Strand is the sign of (to - from) along ex_part edges, median over
## the gene: correct on 80 of 80 genes. intronparts.R builds the intronic-bin
## manifest from the GFFs and the two agree exactly under that convention, so it
## needs no regeneration.
graph_geom <- function(gg) {
  el <- as.data.table(as_edgelist(gg, names = TRUE)); setnames(el, c("from", "to"))
  el$type <- E(gg)$ex_or_in; el$frag <- E(gg)$dexseq_fragment
  pos <- suppressWarnings(setNames(as.numeric(V(gg)$position), V(gg)$name))
  ep <- el[type == "ex_part" & !is.na(frag) & frag != ""]
  if (!nrow(ep)) return(NULL)
  ep[, `:=`(pf = pos[from], pt = pos[to])]
  ep <- ep[!is.na(pf) & !is.na(pt)]
  if (!nrow(ep)) return(NULL)
  d <- ep[, .(start = min(pmin(pf, pt)), end = max(pmax(pf, pt)) - 1L),
          by = .(part = as.integer(frag))]
  setorder(d, start)
  list(geo = d, strand = if (median(ep$pt - ep$pf) > 0) "+" else "-",
       s2p = setNames(ep$frag, ep$from), k2p = setNames(ep$frag, ep$to))
}

## Is there any part OUTWARD of this one, or is it the gene's edge? Emitted here
## rather than recomputed downstream: this is the only stage that holds the
## gene's full part list, and score_outward_walk.R needs the answer but has no
## other reason to load geometry. `geo` is ordered by start, so outward-and-
## decreasing means position 1 and outward-and-increasing means the last row.
gene_edge_flag <- function(a, dec, geo) {
  j <- match(as.integer(a), geo$part)
  if (dec) is.na(j) | j == 1L else is.na(j) | j == nrow(geo)
}

ip <- fread(o$intron_manifest, sep = "\t", showProgress = FALSE)
setkey(ip, gene)

## The rule for the unit immediately outward of a part -- abutting exonic part
## when the tiling is contiguous, else the intronic bin at the end facing the
## walk, else nothing at a gene edge -- is grase::outward_unit(), unit tested
## there. A missing unit is UNMEASURED SPACE, not a measured absence of signal.

one_gene <- function(g) {
  sf <- file.path(o$split_dir, sprintf("%s.bipartition.terminal.txt", g))
  gf <- file.path(o$graphml_dir, sprintf("%s.graphml", g))
  if (!file.exists(sf) || !file.exists(gf)) return(NULL)
  d <- tryCatch(fread(sf, sep = "\t", showProgress = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d[, event := as.character(.I)]
  gg <- tryCatch(read_graph(gf, format = "graphml"), error = function(e) NULL)
  if (is.null(gg)) return(NULL)
  gm <- graph_geom(gg)
  if (is.null(gm)) return(NULL)
  s2p <- gm$s2p; k2p <- gm$k2p; st <- gm$strand; geo <- gm$geo
  ib  <- ip[J(g), .(bin_id, abut_part_lo, abut_part_hi, touches), nomatch = NULL]
  rbindlist(lapply(seq_len(nrow(d)), function(i) {
    ## which graph terminal the side is anchored at: R gives a start, L an end.
    ## `kind` is that semantic label, not the terminal name -- outward_decreasing()
    ## pairs it with strand to get the walk direction.
    ## 5prime / 3prime, not start / end: the free end of an R-anchored side is
    ## the transcription start and of an L-anchored side the transcription end,
    ## and the same vocabulary is used by the bounded.<end> label downstream.
    kind <- if (identical(as.character(d$source[i]), "R")) "5prime" else
            if (identical(as.character(d$sink[i]),   "L")) "3prime" else NA
    if (is.na(kind)) return(NULL)
    dec <- outward_decreasing(kind, st)
    rbindlist(lapply(1:2, function(sd) {
      cs <- called[J(g, d$event[i], sd), nomatch = NULL]
      if (!nrow(cs)) return(NULL)
      p <- if (sd == 1) d$path1[i] else d$path2[i]
      dp <- if (sd == 1) d$setdiff1[i] else d$setdiff2[i]
      anchors <- data.table(anchor = character(0), anchor_kind = character(0))
      tp <- if (is.na(p) || p %in% c("", "NA")) character(0) else
              side_termini(p, kind, s2p, k2p)
      if (length(tp)) anchors <- rbind(anchors,
        data.table(anchor = tp, anchor_kind = "terminus"))
      if (!is.na(dp) && !dp %in% c("", "NA")) {
        dn <- as.integer(gsub("E", "", trimws(unlist(strsplit(dp, ",", fixed = TRUE)))))
        dn <- dn[!is.na(dn)]
        dr <- geo[part %in% dn]
        if (nrow(dr)) {
          oe <- if (dec) dr$part[which.min(dr$start)] else dr$part[which.max(dr$end)]
          anchors <- rbind(anchors, data.table(anchor = pid(oe), anchor_kind = "d_edge"))
        }
      }
      if (!nrow(anchors)) return(NULL)
      flanks <- outward_unit(anchors$anchor, dec, geo, ib)
      out1 <- data.table(gene = g, event = d$event[i], side = sd, kind = kind, strand = st,
                 anchor = anchors$anchor, anchor_kind = anchors$anchor_kind,
                 unit_type = flanks$unit_type, unit_id = flanks$unit_id,
                 gene_edge = gene_edge_flag(anchors$anchor, dec, geo),
                 ref_ex_part = cs$ref_ex_part[1], called = cs$called[1])
      ## the terminus parts themselves
      tn <- anchors[anchor_kind == "terminus"]
      out2 <- if (nrow(tn))
        data.table(gene = g, event = d$event[i], side = sd, kind = kind, strand = st,
                   anchor = tn$anchor, anchor_kind = "terminus_self",
                   unit_type = "exon", unit_id = tn$anchor,
                   gene_edge = gene_edge_flag(tn$anchor, dec, geo),
                   ref_ex_part = cs$ref_ex_part[1], called = cs$called[1]) else NULL
      rbind(out1, out2)
    }))
  }))
}

genes <- unique(called$gene)
cat(sprintf("genes to process: %d  (mc_cores %d)\n", length(genes), o$mc_cores))
L <- parallel::mclapply(genes, one_gene, mc.cores = o$mc_cores)
u <- rbindlist(L[!vapply(L, is.null, TRUE)], fill = TRUE)
cat(sprintf("anchor rows: %d\n", nrow(u)))
n_edge <- sum(is.na(u$unit_id)); u <- u[!is.na(unit_id)]
cat(sprintf("  at a gene edge, no outward flank: %d\n", n_edge))
cat(sprintf("  anchors: %d terminus flanks, %d terminus parts, %d D-edge\n",
            sum(u$anchor_kind == "terminus"), sum(u$anchor_kind == "terminus_self"),
            sum(u$anchor_kind == "d_edge")))
cat(sprintf("  distinct (gene, unit, ref) test units: %d  from %d anchor rows (collapse %.1f%%)\n",
            uniqueN(u[, paste(gene, unit_id, ref_ex_part)]), nrow(u),
            100 * (1 - uniqueN(u[, paste(gene, unit_id, ref_ex_part)]) / nrow(u))))
setorder(u, gene, event, side, anchor_kind, anchor)
fwrite(u, o$out, sep = "\t", quote = FALSE)
cat("wrote ", o$out, "\n", sep = "")
