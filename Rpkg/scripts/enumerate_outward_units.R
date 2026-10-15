#!/usr/bin/env Rscript
## Enumerate the OUTWARD test units for the EXTERNAL-event boundary arm.
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
## Testing only the tile OUTWARD of each start answers "does the effect continue
## past it" but never "is this start used at all". With both, the walk along a
## side alternates
##     U_n  t_n  U_n-1  t_n-1  ...  U_1  t_1  [D]
## from outermost inward, and the real start is the outermost t_i that is called
## with everything outward of it uncalled -- a localization rather than a
## bounded/unbounded flag. It also rescues sides whose outermost tile is below
## the read floor: if t_n itself is powered and uncalled, that start is unused
## regardless of what lies beyond it. 72.4% of terminus parts are not reachable
## as any anchor's outward tile, so this is not redundant; all are exonic, so
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
## consecutive parts (PDCD4 event 11 side 1: 001, 002, 003, 004), so the tiles
## outward of them collapse onto one another and onto the termini themselves.
## Emitting per terminus would test the same unit repeatedly under different
## labels. The test unit is (gene, unit, ref_ex_part) -- reference counts are
## just the counts of ref_ex_part, so sides sharing a reference part share a
## test. The side table carries the one-to-many mapping back.
##
## Usage:
##   Rscript scripts/enumerate_outward_units.R --results <annotated.txt> \
##     --split_dir <bipartition.filtered> --graphml_dir <graphml.v34> \
##     --exon_parts <geom_parts.tsv> --intron_manifest <intronparts_manifest.tsv> \
##     --out <units.tsv> [--mc_cores 8]
suppressPackageStartupMessages({library(data.table); library(igraph); library(optparse)
                                library(grase)})

o <- parse_args(OptionParser(option_list = list(
  make_option("--results",         type = "character"),
  make_option("--split_dir",       type = "character"),
  make_option("--graphml_dir",     type = "character"),
  make_option("--exon_parts",      type = "character"),
  make_option("--intron_manifest", type = "character"),
  make_option("--out",             type = "character"),
  make_option("--mc_cores",        type = "integer", default = 8L))))
for (r in c("results","split_dir","graphml_dir","exon_parts","intron_manifest","out"))
  if (is.null(o[[r]])) stop("--", r, " is required")
pid <- function(x) sprintf("%03d", as.integer(x))

res <- fread(o$results, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene","event")))
res[, side := fifelse(comparison == "diff1_vs_ref", 1L, 2L)]
called <- unique(res[significant == TRUE, .(gene, event, side, ref_ex_part)])
setkey(called, gene, event, side)
cat(sprintf("called sides: %d over %d genes\n", nrow(called), uniqueN(called$gene)))

ex <- fread(o$exon_parts, showProgress = FALSE)
ex[, `:=`(start = as.numeric(start), end = as.numeric(end), part = as.integer(part))]
setorder(ex, gene, start)
setkey(ex, gene, part)

ip <- fread(o$intron_manifest, sep = "\t", showProgress = FALSE)
setkey(ip, gene)

## The rule for the unit immediately outward of a part -- abutting exonic part
## when the tiling is contiguous, else the intronic bin at the end facing the
## walk, else nothing at a gene edge -- is grase::outward_unit(), unit tested
## there. A missing unit is UNMEASURED SPACE, not a measured absence of signal.

one_gene <- function(g) {
  sf <- file.path(o$split_dir, sprintf("%s.bipartition.external.txt", g))
  gf <- file.path(o$graphml_dir, sprintf("%s.graphml", g))
  if (!file.exists(sf) || !file.exists(gf)) return(NULL)
  d <- tryCatch(fread(sf, sep = "\t", showProgress = FALSE), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d[, event := as.character(.I)]
  gg <- tryCatch(read_graph(gf, format = "graphml"), error = function(e) NULL)
  if (is.null(gg)) return(NULL)
  el <- as.data.table(as_edgelist(gg, names = TRUE)); setnames(el, c("from","to"))
  el[, `:=`(type = E(gg)$ex_or_in, frag = E(gg)$dexseq_fragment)]
  epg <- el[type == "ex_part" & !is.na(frag) & frag != ""]
  if (!nrow(epg)) return(NULL)
  s2p <- setNames(epg$frag, epg$from); k2p <- setNames(epg$frag, epg$to)
  st <- ex[J(g), strand, nomatch = NULL][1]
  if (is.na(st)) return(NULL)
  geo <- ex[J(g), .(part, start, end), nomatch = NULL]
  ib  <- ip[J(g), .(bin_id, abut_part_lo, abut_part_hi, touches), nomatch = NULL]
  rbindlist(lapply(seq_len(nrow(d)), function(i) {
    ## which graph terminal the side is anchored at: R gives a start, L an end.
    ## `kind` is that semantic label, not the terminal name -- outward_decreasing()
    ## pairs it with strand to get the walk direction.
    kind <- if (identical(as.character(d$source[i]), "R")) "start" else
            if (identical(as.character(d$sink[i]),   "L")) "end" else NA
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
        dr <- ex[J(g, dn), nomatch = NULL]
        if (nrow(dr)) {
          oe <- if (dec) dr$part[which.min(dr$start)] else dr$part[which.max(dr$end)]
          anchors <- rbind(anchors, data.table(anchor = pid(oe), anchor_kind = "d_edge"))
        }
      }
      if (!nrow(anchors)) return(NULL)
      tiles <- outward_unit(anchors$anchor, dec, geo, ib)
      out1 <- data.table(gene = g, event = d$event[i], side = sd, kind = kind, strand = st,
                 anchor = anchors$anchor, anchor_kind = anchors$anchor_kind,
                 unit_type = tiles$unit_type, unit_id = tiles$unit_id,
                 ref_ex_part = cs$ref_ex_part[1])
      ## the terminus parts themselves
      tn <- anchors[anchor_kind == "terminus"]
      out2 <- if (nrow(tn))
        data.table(gene = g, event = d$event[i], side = sd, kind = kind, strand = st,
                   anchor = tn$anchor, anchor_kind = "terminus_self",
                   unit_type = "exon", unit_id = tn$anchor,
                   ref_ex_part = cs$ref_ex_part[1]) else NULL
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
cat(sprintf("  at a gene edge, no outward tile: %d\n", n_edge))
cat(sprintf("  anchors: %d terminus tiles, %d terminus parts, %d D-edge\n",
            sum(u$anchor_kind == "terminus"), sum(u$anchor_kind == "terminus_self"),
            sum(u$anchor_kind == "d_edge")))
cat(sprintf("  distinct (gene, unit, ref) test units: %d  from %d anchor rows (collapse %.1f%%)\n",
            uniqueN(u[, paste(gene, unit_id, ref_ex_part)]), nrow(u),
            100 * (1 - uniqueN(u[, paste(gene, unit_id, ref_ex_part)]) / nrow(u))))
setorder(u, gene, event, side, anchor_kind, anchor)
fwrite(u, o$out, sep = "\t", quote = FALSE)
cat("wrote ", o$out, "\n", sep = "")
