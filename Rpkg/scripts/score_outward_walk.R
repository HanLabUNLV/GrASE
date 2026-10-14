#!/usr/bin/env Rscript
## Score the outward walk for TSS/TTS calls: is the differential effect
## localized at an annotated transcription start / end, or does it continue
## past it?
##
## THE QUESTION. A GrASE TSS/TTS call says a side's distinct set D differs in
## usage between conditions. It does NOT say where transcription begins. This
## walks outward from the side's annotated termini and asks, with the SAME test
## and the SAME shared reference, where signal stops.
##
## UNITS. For a called side, its terminus parts t_1 .. t_n are the ex_part edges
## sourced at the first vertex after R (sunk at the last vertex before L for a
## TTS), ordered OUTERMOST FIRST. For each t_i, U_i is the tile immediately
## outward of it: the abutting exonic part when the tiling is contiguous there,
## otherwise the intronic part between them. U_i does not exist when t_i is the
## gene's outermost part, or when the intronic gap was too short to bin under
## T >= 2W. Every t_i and U_i is tested substituted for D against the side's own
## ref_ex_part, so static contributions from other isoforms cancel.
##
## THE INFERRED START is t_k, the outermost terminus that is both powered and
## called. Everything is judged relative to it.
##
## VERDICTS, in evaluation order:
##   no powered start            no t_i clears the floor; nothing is judgeable
##   no start carries signal     some t_i powered, none powered AND called
##   signal outward              U_k is powered and called: no drop
##   drop: measured              U_k is powered and uncalled -- drop OBSERVED
##   drop: gene edge             U_k does not exist, t_k is the gene's edge
##   drop: tile below floor      U_k exists but is under the floor; the drop is
##                               ASSERTED, not measured
##
## Only `drop: measured` is evidence. The other two drop rows are claims about
## space that was never measured, and they are what made every earlier headline
## figure move when a threshold changed. Report them separately, never summed
## into a single "localized" fraction.
##
## ONLY U_k CAN CONTRADICT. The boundary claim is about t_k, so a called tile
## further out is a separate event -- another start, or an effect that skips the
## intervening tile -- and does not negate a drop immediately beyond the start.
## The alternative (any U_i with i <= k contradicts) was tried and moves 1.1% of
## sides; it is strictly more conservative and was dropped as superseded.
##
## TWO FLAGS THAT MUST TRAVEL WITH THE NUMBERS.
##
## 1. `taut`. When a side's ENTIRE distinct set is one terminus part, testing
##    that terminus re-runs the test that selected the side, so `taut` marks the
##    sides where the start test is CIRCULAR. Those units are called 100.0% of
##    the time by construction, and such sides can never land in `no start
##    carries signal`. They are a third of called sides, so NEVER report a
##    figure pooled over this flag -- split on it.
##
##    `taut == FALSE` means ONLY that the start test is independent of the call.
##    It is NOT "D is larger than its terminus": two thirds of it is a
##    SINGLE-part D sitting inward of every terminus (13,398 of 20,755), with
##    just 35% having a multi-part D. Do not label that stratum by D's size.
##
##    The structurally meaningful split is `reaches` -- whether D overlaps any
##    terminus part at all. Only 50.3% of called sides do; for the other half
##    the annotated start is an extrapolation from the graph that D never
##    touches. Conditional on the question being answerable, the measured-drop
##    rate is 23.6-29.9% across EVERY class, so the classes differ in how often
##    the question is answerable, not in the answer.
##
## THE TWO DENOMINATORS, and always name which one is meant:
##   pct_of_sides       drop: measured over ALL called sides in the cell
##   pct_of_uk_powered  drop: measured over the sides where U_k cleared the
##                      floor, i.e. where the outward space was MEASURED at all.
##                      `uk_powered` TRUE partitions into `drop: measured` and
##                      `signal outward of the start`, nothing else.
##
## 2. `contrast`. Contrasts are independent tests. Pooling them mixes two
##    populations and hides the between-contrast replication rate, which for
##    the drop verdict is only about half (47.7%) against strong agreement on
##    continuation.
##
## THE FLOOR IS d_support, NOT A POOLED MEAN. d_support is the max over the
## contrast's two group means, matching add_support() in exontest.R (max() not
## min() is deliberate there). A mean pooled over samples dilutes a unit that is
## on in one group and off in the other -- exactly the on/off case max() exists
## to protect -- and inflates the untestable and below-floor categories by
## roughly 5% of units and 7 points of `no powered start`.
##
## A PER-BASE DENSITY FLOOR WAS TRIED AND REMOVED. It divided d_support by
## effective length (len + readlen - 1, because dexseq_count counts a read that
## OVERLAPS a bin, so dividing by len alone inflates short bins 2.9x at 26bp
## against 1.4x at 137bp). It gave measured drop 4.5% / 3.5% in the independent
## stratum and 20.3% / 16.5% in the circular one, i.e. the same conclusion, so
## it served as a robustness check and is not kept as a mode. DO NOT reintroduce
## density as the drop CRITERION -- the per-unit means are pooled over conditions,
## so a dens(U_k)/dens(t_k) ratio is total coverage, not differential, and is
## confounded by pi exactly as the old step statistic was.
##
## Usage:
##   Rscript scripts/score_outward_walk.R \
##     --side_table intronparts/outward_side_table.tsv \
##     --results    intronparts/test.outward/test_outward.exoncnt_betabinom_EBmap.annotated.txt \
##     --exonic     grase_activation_merged/test_bipartition.merged.TSSTTS_betabinom_EBmap.annotated.txt \
##     --group_means intronparts/outward_unit_group_means.tsv \
##     --geom       tss_edge/geom_parts.tsv \
##     --intronparts intronparts/intronparts_manifest.tsv \
##     --contrasts  'CD4_N_STIM:CD4,CD8_N_STIM:CD8' \
##     --out        intronparts/outward_walk_verdict.tsv
##
## If --group_means does not exist it is built from --combined with awk (the
## combined table runs to gigabytes; a line-by-line R loop over it is not
## viable, same reason as build_flank_counts.R).
suppressPackageStartupMessages({ library(data.table); library(optparse)
                                 library(grase) })

opt_list <- list(
  make_option("--side_table",  type = "character",
              help = "outward_side_table.tsv from build_outward_counts.R"),
  make_option("--results",     type = "character",
              help = "exontest annotated output for the outward units"),
  make_option("--exonic",      type = "character",
              help = "main TSS/TTS annotated results; supplies the called sides and setdiff"),
  make_option("--group_means", type = "character", default = NULL,
              help = "per (gene, unit, group) mean counts; built from --combined if absent"),
  make_option("--combined",    type = "character", default = NULL,
              help = "outward combined count table, only needed to build --group_means"),
  make_option("--geom",        type = "character", help = "geom_parts.tsv (exonic part coordinates)"),
  make_option("--intronparts", type = "character", help = "intronparts_manifest.tsv"),
  make_option("--contrasts",   type = "character",
              help = "trt:ref pairs, comma separated, e.g. 'CD4_N_STIM:CD4,CD8_N_STIM:CD8'"),
  make_option("--min_reads",   type = "double", default = 10,
              help = "read floor on d_support [default %default]"),
  make_option("--padj",        type = "double", default = 0.01),
  make_option("--min_dpi",     type = "double", default = 0.1),
  make_option("--out",         type = "character", help = "per-side verdict table")
)
a <- parse_args(OptionParser(option_list = opt_list))
for (req in c("side_table", "results", "exonic", "geom", "intronparts", "contrasts", "out"))
  if (is.null(a[[req]])) stop("--", req, " is required")

CTR <- lapply(strsplit(strsplit(a$contrasts, ",", fixed = TRUE)[[1]], ":", fixed = TRUE),
              function(x) { if (length(x) != 2) stop("bad contrast: ", paste(x, collapse = ":")); trimws(x) })
names(CTR) <- vapply(CTR, function(x) paste0(x[1], "_vs_", x[2]), "")
cat("contrasts:", paste(names(CTR), collapse = ", "), "\n")

## ---- exonic part geometry --------------------------------------------------
ex <- fread(a$geom, showProgress = FALSE)
ex[, `:=`(start = as.numeric(start), end = as.numeric(end), part = as.integer(part))]
setorder(ex, gene, start)
ex[, `:=`(nxt = shift(part, type = "lead"), prv = shift(part, type = "lag")), by = gene]

## ---- d_support: max over the contrast's two group means -------------------
if (is.null(a$group_means) || !file.exists(a$group_means)) {
  if (is.null(a$combined)) stop("--group_means is absent; pass --combined so it can be built")
  gmf <- if (is.null(a$group_means)) tempfile(fileext = ".tsv") else a$group_means
  cat("building group means from", a$combined, "\n")
  ## columns of the combined table: gene event source sink ref_ex_part setdiff1
  ## setdiff2 sample ref diff1 diff2 groups diff1_source diff2_source. The unit
  ## count sits in diff1 (field 10); diff2 is NA for a single-unit row.
  cmd <- sprintf(paste0("awk -F'\\t' 'NR>1 && $10!=\"NA\" && $10!=\"\" {k=$1\"\\t\"$2\"\\t\"$12; ",
                        "s[k]+=$10; n[k]++} END{print \"gene\\tsynth_event\\tgroups\\tgmean\\tgn\"; ",
                        "for(k in s) printf \"%%s\\t%%.6f\\t%%d\\n\", k, s[k]/n[k], n[k]}' %s > %s"),
                 shQuote(a$combined), shQuote(gmf))
  if (system(cmd) != 0) stop("awk group-mean aggregation failed")
  a$group_means <- gmf
}
gm <- fread(a$group_means, sep = "\t", showProgress = FALSE,
            colClasses = list(character = c("gene", "synth_event", "groups")))
miss_grp <- setdiff(unlist(CTR), unique(gm$groups))
if (length(miss_grp))
  stop("groups named in --contrasts are absent from the group means: ", paste(miss_grp, collapse = ", "))
dsup <- rbindlist(lapply(names(CTR), function(cn)
  gm[groups %in% CTR[[cn]], .(contrast = cn, d_support = max(gmean)), by = .(gene, synth_event)]))

## ---- units, ordered outermost first ---------------------------------------
st <- fread(a$side_table, sep = "\t", showProgress = FALSE,
            colClasses = list(character = c("gene", "event", "synth_event")))
u <- unique(st[anchor_kind %in% c("terminus", "terminus_self"),
               .(gene, event, side, kind, strand, anchor, anchor_kind, synth_event)])
u[, anchor_n := as.integer(anchor)]
## out_dec: outward is the direction of DECREASING coordinate, so rank 1 is the
## lowest part number there and the highest otherwise.
u[, out_dec := outward_decreasing(kind, strand)]
u[, trank := frank(fifelse(out_dec, anchor_n, -anchor_n), ties.method = "dense"),
  by = .(gene, event, side)]
u <- u[, .(contrast = names(CTR)), by = names(u)]

res <- fread(a$results, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene", "event")))
res[, called := padj < a$padj & abs(delta_pi) >= a$min_dpi & lfc_diff_net > 0]
u <- merge(u, res[, .(gene, synth_event = event, contrast, called)],
           by = c("gene", "synth_event", "contrast"), all.x = TRUE)
u <- merge(u, dsup, by = c("gene", "synth_event", "contrast"), all.x = TRUE)
u[, powered := !is.na(d_support) & d_support >= a$min_reads & !is.na(called)]
cat(sprintf("floor: d_support >= %g reads (max over the contrast's two groups)\n", a$min_reads))

## ---- called sides, and the tautology flag ---------------------------------
e2 <- fread(a$exonic, sep = "\t", showProgress = FALSE,
            colClasses = list(character = c("gene", "event")))
e2[, side := fifelse(comparison == "diff1_vs_ref", 1L, 2L)]
cl <- e2[significant == TRUE, .(gene, event, side, contrast,
                                sd = fifelse(side == 1L, setdiff1, setdiff2))]
cl <- cl[contrast %in% names(CTR)]
if (!nrow(cl)) stop("no significant sides for these contrasts in ", a$exonic)
norm_parts <- function(x) {
  p <- strsplit(gsub("[^0-9,]", "", x), ",", fixed = TRUE)
  vapply(p, function(v) { v <- v[nzchar(v)]
    if (!length(v)) NA_character_ else
      paste(sort(sprintf("%03d", as.integer(v))), collapse = ",") }, "")
}
cl[, dset := norm_parts(sd)]
cat(sprintf("called sides (side x contrast): %d\n", nrow(cl)))
u <- merge(u, cl[, .(gene, event, side, contrast, dset)],
           by = c("gene", "event", "side", "contrast"))

## ---- assemble t_i beside U_i on trank -------------------------------------
ts <- merge(u[anchor_kind == "terminus_self"], ex[, .(gene, part, prv, nxt)],
            by.x = c("gene", "anchor_n"), by.y = c("gene", "part"), all.x = TRUE)
ts[, gene_edge := is.na(fifelse(out_dec, prv, nxt))]
## eq_D: this terminus part IS the side's entire distinct set, so its test is
## the same test that selected the side.
ts[, eq_D := !is.na(dset) & dset == sprintf("%03d", anchor_n)]
## reaches: does D include THIS terminus part? Collapsed per side below. Half of
## called sides have a D that touches no terminus at all, so the start they are
## labelled with is an extrapolation D never reaches.
ts[, in_D := !is.na(dset) &
     mapply(function(d, p) p %in% strsplit(d, ",", fixed = TRUE)[[1]],
            dset, sprintf("%03d", anchor_n))]
w <- merge(ts[, .(gene, event, side, contrast, trank, gene_edge, eq_D, in_D,
                  part_called = called, part_powered = powered)],
           u[anchor_kind == "terminus",
             .(gene, event, side, contrast, trank, tile_called = called,
               tile_powered = powered, tile_present = TRUE)],
           by = c("gene", "event", "side", "contrast", "trank"), all = TRUE)
setorder(w, contrast, gene, event, side, trank)

v <- w[, {
  k <- which(part_powered & part_called)[1]
  .(n_term = .N, n_powered = sum(part_powered, na.rm = TRUE), start_rank = k,
    taut      = any(eq_D & part_powered & part_called, na.rm = TRUE),
    reaches   = any(in_D, na.rm = TRUE),
    ## start_in_D: is t_k, the terminus the verdict is anchored on, a MEMBER of
    ## D? If not, both sides of the bipartition enter on that same part, so the
    ## split is internal divergence and a drop beyond it carries no claim about
    ## alternative TSS/TTS. A positive call requires this. `reaches` is weaker:
    ## it only asks whether SOME terminus is in D.
    start_in_D = if (is.na(k)) NA else isTRUE(in_D[k]),
    uk_present= if (is.na(k)) NA else isTRUE(tile_present[k]),
    uk_powered= if (is.na(k)) NA else isTRUE(tile_powered[k]),
    uk_called = if (is.na(k)) NA else isTRUE(tile_called[k]),
    uk_edge   = if (is.na(k)) NA else isTRUE(gene_edge[k]))
}, by = .(contrast, gene, event, side)]

v[, verdict := fifelse(n_powered == 0, "no powered start",
                fifelse(is.na(start_rank), "no start carries signal",
                 fifelse(uk_powered & uk_called, "signal outward of the start",
                  fifelse(uk_powered, "drop: measured",
                   fifelse(!uk_present & uk_edge, "drop: gene edge",
                           "drop: tile below floor")))))]
v[, drop_measured := verdict == "drop: measured"]
v[, drop_any      := grepl("^drop", verdict)]
## The positive call: a drop AND the anchoring terminus is in D.
v[, tsstts := drop_any & !is.na(start_in_D) & start_in_D]

## ---- the assertion that matters -------------------------------------------
## A side that lost its units to a failed join would be scored `no powered
## start` -- indistinguishable from a real lack of depth, and that category is
## one of the reported numbers.
lost <- cl[!unique(v[, .(gene, event, side, contrast)]),
           on = .(gene, event, side, contrast)]
if (nrow(lost))
  stop(sprintf(paste0("FATAL: %d called sides have no terminus unit, e.g. %s\n",
                      "They would be scored 'no powered start', indistinguishable from\n",
                      "low depth. Check the enumeration covers every called side."),
               nrow(lost), paste(head(lost[, paste(gene, event, side)], 5), collapse = ", ")))
cat(sprintf("all %d called sides resolve to terminus units\n", nrow(cl)))

fwrite(v, a$out, sep = "\t")
cat("wrote ", a$out, "\n", sep = "")

## ---- summary: per contrast, split on taut, never pooled -------------------
show <- function(d, lab) {
  cat("\n=== ", lab, "  (n = ", nrow(d), ")\n", sep = "")
  p <- d[, .(sides = .N), by = verdict][order(-sides)]
  p[, pct := round(100 * sides / nrow(d), 1)]
  print(p)
  cat(sprintf("drop MEASURED %d (%.1f%%) | drop any %d (%.1f%%) | signal outward %d (%.1f%%)\n",
              sum(d$drop_measured), 100 * mean(d$drop_measured),
              sum(d$drop_any), 100 * mean(d$drop_any),
              sum(d$verdict == "signal outward of the start"),
              100 * mean(d$verdict == "signal outward of the start")))
  cat(sprintf("TSS/TTS (drop AND t_k in D) %d (%.1f%%)\n",
              sum(d$tsstts), 100 * mean(d$tsstts)))
}
for (cn in names(CTR)) {
  show(v[contrast == cn & taut == FALSE], paste(cn, "-- start test INDEPENDENT of the call"))
  show(v[contrast == cn & taut == TRUE ], paste(cn, "-- terminus unit IS D (start test CIRCULAR)"))
}
cat("\n=== by whether D reaches any terminus part (the structural split)\n")
## Two denominators, both named after the column they come from. uk_powered is
## the sides where U_k clears the floor, so the outward space was MEASURED; when
## it is TRUE the verdict is either `drop: measured` or `signal outward of the
## start`, and those two partition it.
rr <- v[, .(sides = .N,
            drop_measured = sum(drop_measured),
            pct_of_sides  = round(100 * mean(drop_measured), 1),
            uk_powered    = sum(uk_powered, na.rm = TRUE),
            pct_of_uk_powered = round(100 * sum(drop_measured) /
                                      sum(uk_powered, na.rm = TRUE), 1)),
        by = .(contrast, reaches, taut)]
setorder(rr, contrast, -reaches, -taut)
print(rr)
cat("\n=== NOTE no pooled figure is reported on purpose: contrasts are independent",
    "\n    tests and the taut strata behave oppositely. `taut == FALSE` means the",
    "\n    start test is independent of the call -- NOT that D is larger than its",
    "\n    terminus. Quote `drop: measured` with its denominator named.\n")
