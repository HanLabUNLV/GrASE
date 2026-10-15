#!/usr/bin/env Rscript
## Score the outward walk for TSS/TTS calls: is the differential effect
## localized at an annotated transcription start / end, or does it continue
## past it?
##
## THE QUESTION. A GrASE TSS/TTS call says a partition's distinct set D differs in
## usage between conditions. It does NOT say where transcription begins. This
## walks outward from D to the paritition's annotated termini and asks, with the SAME test
## and the SAME shared reference, where signal stops.
##
## UNITS. For a called partition, its terminus parts t_1 .. t_n are the ex_part edges
## sourced at the first vertex after R (sunk at the last vertex before L for a
## TTS), ordered OUTERMOST FIRST. For each t_i, U_i is the tile immediately
## outward of it: the abutting exonic part when the tiling is contiguous there,
## otherwise the intronic part between them. U_i does not exist when t_i is the
## gene's outermost part, or when the gene has no intronic bins at all --
## intronic_parts() returns empty for a single-part gene, for overlapping parts,
## and for a gene whose parts all abut. There is NO minimum gap length:
## `has <- ge >= gs` bins every gap of one base or more, a short intron whole
## and a long one as two end bins. 
## Every t_i and U_i is tested substituted for D against the side's own
## ref_ex_part. 
##
##
## VERDICTS, in evaluation order:
##   no start above the read floor
##                       no t_i clears the min_read_support; nothing is judgeable
##   no start carries signal
##                       some t_i clears the min_read_support, but no t_i is called
##   signal outward      U_k clears the floor and is called: no drop
##   drop: measured      U_k clears the floor and is uncalled -- drop OBSERVED
##   drop: gene edge     U_k does not exist, t_k is the gene's edge
##   drop: tile below floor
##                       U_k exists but is under the floor; the drop is
##                       ASSERTED, not measured
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
#
## WHAT MUST TRAVEL WITH THE NUMBERS.
##
## 1. The five positive kinds, as five rows. They differ in how much of the
##    outward space was quantified -- from a tile tested above the floor to no
##    tile at all -- and on the activation panel the untestable and gene-edge
##    kinds are together the majority of positives. Pooling them into one
##    figure is what made every earlier headline number move when a threshold
##    changed.
##
## 2. `contrast`. Contrasts are independent tests. Pooling them mixes two
##    populations and hides the between-contrast replication rate, which for
##    the drop verdict is only about half (47.7%) against strong agreement on
##    continuation.
##
## THE FLOOR IS d_support, NOT A POOLED MEAN. d_support is the max over the
## contrast's two group means, matching add_support() in exontest.R (max() not
## min() is deliberate there). A mean pooled over samples dilutes a unit that is
## on in one group and off in the other. 
##
## Usage:
##   Rscript scripts/score_outward_walk.R \
##     --side_table intronparts/outward_side_table.tsv \
##     --results    intronparts/test.outward/test_outward.exoncnt_betabinom_EBmap.annotated.txt \
##     --exonic     grase_activation_exsj/test_bipartition.exsj.external_betabinom_EBmap.annotated.txt \
##     --group_means intronparts/outward_unit_group_means.tsv \
##     --intronparts intronparts/intronparts_manifest.tsv \
##     --contrasts  'CD4_N_STIM:CD4,CD8_N_STIM:CD8' \
##     --out        intronparts/outward_walk_verdict.tsv
##
## If --group_means does not exist it is built from --combined with awk (the
## combined table runs to gigabytes; a line-by-line R loop over it is not
## viable).
suppressPackageStartupMessages({ library(data.table); library(optparse)
                                 library(grase) })

opt_list <- list(
  make_option("--side_table",  type = "character",
              help = "outward_side_table.tsv from build_outward_counts.R"),
  make_option("--results",     type = "character",
              help = paste("exontest annotated output for the outward units.",
                           "Comma-separated for several files, e.g. one per",
                           "contrast from a parallel run -- they are read and",
                           "row-bound here, so nothing has to concatenate them",
                           "on disk first.")),
  make_option("--exonic",      type = "character",
              help = "main TSS/TTS annotated results; supplies the called sides and setdiff"),
  make_option("--group_means", type = "character", default = NULL,
              help = "per (gene, unit, group) mean counts; built from --combined if absent"),
  make_option("--combined",    type = "character", default = NULL,
              help = "outward combined count table, only needed to build --group_means"),
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
for (req in c("side_table", "results", "exonic", "intronparts", "contrasts", "out"))
  if (is.null(a[[req]])) stop("--", req, " is required")

CTR <- lapply(strsplit(strsplit(a$contrasts, ",", fixed = TRUE)[[1]], ":", fixed = TRUE),
              function(x) { if (length(x) != 2) stop("bad contrast: ", paste(x, collapse = ":")); trimws(x) })
names(CTR) <- vapply(CTR, function(x) paste0(x[1], "_vs_", x[2]), "")
cat("contrasts:", paste(names(CTR), collapse = ", "), "\n")

## ---- no exonic part geometry is read ---------------------------------------
## No geometry is read here. gene_edge -- whether anything exists outward of a
## terminus -- now arrives in the side table from enumerate_outward_units.R,
## which is the only stage that holds a gene's full part list and already has
## the graph open. This script used to fread a GFF-derived table
## (tss_edge/geom_parts.tsv) solely to shift part numbers within a gene and
## recover that one boolean.

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
## The floor comes from grase::add_support(), the same function the exon test
## applies, rather than a second max() written here. This arm writes every unit
## into diff1 and leaves diff2 NA (verified: 0 of 26,357,294 rows), so `side` is
## constant and the awk group means above need no side key -- but the
## max-over-groups rule and the ctr_groups handling are now shared, not copied.
d_grp <- data.frame(gene   = as.character(gm$gene),
                    event  = as.character(gm$synth_event),
                    groups = as.character(gm$groups),
                    m      = as.numeric(gm$gmean),
                    side   = "diff1", stringsAsFactors = FALSE)

## ---- units, ordered outermost first ---------------------------------------
st <- fread(a$side_table, sep = "\t", showProgress = FALSE,
            colClasses = list(character = c("gene", "event", "synth_event")))
if (!"gene_edge" %in% names(st))
  stop("--side_table has no gene_edge column. It is emitted by\n",
       "  enumerate_outward_units.R and carried through build_outward_counts.R.\n",
       "  A side table built before that change must be regenerated; this script\n",
       "  no longer reads geometry and cannot recover the flag itself.",
       call. = FALSE)
u <- unique(st[anchor_kind %in% c("terminus", "terminus_self"),
               .(gene, event, side, kind, strand, anchor, anchor_kind, synth_event,
                 gene_edge)])
u[, anchor_n := as.integer(anchor)]
## out_dec: outward is the direction of DECREASING coordinate, so rank 1 is the
## lowest part number there and the highest otherwise.
u[, out_dec := outward_decreasing(kind, strand)]
u[, trank := frank(fifelse(out_dec, anchor_n, -anchor_n), ties.method = "dense"),
  by = .(gene, event, side)]
u <- u[, .(contrast = names(CTR)), by = names(u)]

## SEVERAL RESULTS FILES ARE ALLOWED. The LRT is run per contrast (see
## run_outward_lrt_parallel.sh), and this script immediately works per contrast
## anyway, so requiring one file meant concatenating ~1 GB on disk purely to
## satisfy the interface and then re-splitting it here. Read them directly.
res_files <- trimws(strsplit(as.character(a$results), ",", fixed = TRUE)[[1]])
res_files <- res_files[nzchar(res_files)]
if (!length(res_files)) stop("--results is empty", call. = FALSE)
miss_f <- res_files[!file.exists(res_files)]
if (length(miss_f))
  stop("no such results file: ", paste(miss_f, collapse = ", "), call. = FALSE)
res <- rbindlist(lapply(res_files, function(f)
         fread(f, sep = "\t", showProgress = FALSE,
               colClasses = list(character = c("gene", "event")))),
       use.names = TRUE)
if (length(res_files) > 1)
  cat(sprintf("results: %d files, %d rows, %d contrasts\n",
              length(res_files), nrow(res), uniqueN(res$contrast)))
## Two files covering the same contrast would double every row of it, which
## would not error anywhere downstream -- it would just inflate the counts.
dupc <- res[, .(files = uniqueN(.I)), by = .(contrast, gene, event, comparison)][files > 1]
if (nrow(dupc))
  stop(nrow(dupc), " (contrast, gene, event, comparison) keys appear more than once ",
       "across --results; the files overlap and the verdict counts would be inflated",
       call. = FALSE)
res[, called := padj < a$padj & abs(delta_pi) >= a$min_dpi & lfc_diff_net > 0]
## exontest.R already wrote a d_support column into --results. Recompute it with
## the shared function and CHECK, rather than trusting either copy: a silent
## disagreement here would move the floor and every verdict with it.
res_own <- if ("d_support" %in% names(res)) as.numeric(res$d_support) else NULL
res <- as.data.frame(res)
res <- grase::add_support(res, d_grp, CTR)
if (!is.null(res_own)) {
  cmp <- !is.na(res_own) & !is.na(res$d_support)
  if (any(cmp)) {
    worst <- max(abs(res_own[cmp] - res$d_support[cmp]))
    cat(sprintf("d_support vs the column in --results: %d compared, max abs diff %.3g\n",
                sum(cmp), worst))
    if (worst > 1e-3)
      stop("recomputed d_support disagrees with the column in --results by ",
           signif(worst, 4), " reads; the two are meant to be the same statistic")
  }
}
res <- as.data.table(res)

## SUPPORT FOR EVERY ENUMERATED UNIT, including the ones exontest emitted no row
## for at all. It drops a unit with no reads entirely rather than writing an NA
## row -- on the activation panel 4,244 of 68,818 counted units (6.2%), of which
## 2,745 are termini, 78.5% with a mean of exactly zero. Taking support from
## --results alone would leave those a bare NA, indistinguishable from a unit
## lost to a bad join; for a TILE that is the difference between the strongest
## boundary evidence and an unjudgeable row.
allu <- unique(u[, .(gene, event = synth_event, contrast,
                     comparison = "diff1_vs_ref")])
allu <- as.data.table(grase::add_support(as.data.frame(allu), d_grp, CTR))
u <- merge(u, res[, .(gene, synth_event = event, contrast, called)],
           by = c("gene", "synth_event", "contrast"), all.x = TRUE)
u <- merge(u, allu[, .(gene, synth_event = event, contrast, d_support)],
           by = c("gene", "synth_event", "contrast"), all.x = TRUE)
u[, tested := !is.na(called)]
u[, min_read_support := !is.na(d_support) & d_support >= a$min_reads]
## zero_count separates "measured and empty" from "never measured". A unit whose
## counts are in the table with a mean of zero IS a measurement: it says no
## reads lie there. One with no support value at all says nothing.
u[, zero_count := !is.na(d_support) & d_support == 0]
cat(sprintf("floor: d_support >= %g reads (max over the contrast's two groups)\n", a$min_reads))
cat(sprintf("units: %d enumerated, %d tested, %d with support but no test row, %d of those zero-count\n",
            nrow(u), sum(u$tested), sum(!u$tested & !is.na(u$d_support)),
            sum(!u$tested & u$zero_count)))

## ---- called sides, and D's relation to each terminus ----------------------
e2 <- fread(a$exonic, sep = "\t", showProgress = FALSE,
            colClasses = list(character = c("gene", "event")))
e2[, side := fifelse(comparison == "diff1_vs_ref", 1L, 2L)]
cl <- e2[significant == TRUE, .(gene, event, side, contrast,
                                sd = fifelse(side == 1L, setdiff1, setdiff2))]
cl <- cl[contrast %in% names(CTR)]
if (!nrow(cl)) stop("no significant sides for these contrasts in ", a$exonic)
## isTRUE() is scalar; this is the vectorised form. NA becomes FALSE, so an
## unknown is never silently treated as a yes.
isTRUE_v <- function(x) !is.na(x) & x

cat(sprintf("called sides (side x contrast): %d\n", nrow(cl)))
u <- merge(u, unique(cl[, .(gene, event, side, contrast)]),
           by = c("gene", "event", "side", "contrast"))

## ---- assemble t_i beside U_i on trank -------------------------------------
ts <- u[anchor_kind == "terminus_self"]
w <- merge(ts[, .(gene, event, side, contrast, trank, gene_edge,
                  part_called = called, part_min_read_support = min_read_support,
                  part_tested = tested, part_zero = zero_count)],
           u[anchor_kind == "terminus",
             .(gene, event, side, contrast, trank, tile_called = called,
               tile_min_read_support = min_read_support, tile_present = TRUE,
               tile_tested = tested, tile_zero = zero_count, tile_ds = d_support)],
           by = c("gene", "event", "side", "contrast", "trank"), all = TRUE)
setorder(w, contrast, gene, event, side, trank)

## ---- the verdict, PER TERMINUS --------------------------------------------
## There is no single anchor. Every terminus is its own test:
##
##   t_ok   <- part_min_read_support & part_called      signal at this terminus
##   u_dead <- !tile_min_read_support | !tile_called    none immediately outward
##   positive for the side <- any(t_ok & u_dead)
##
## An earlier version took k <- which(part_min_read_support & part_called)[1]
## and judged only that one terminus. That forced a single-start model onto
## sides with up to nine termini, let a below-floor OUTER terminus silently
## shift the anchor inward on 34.3% of sides, and needed a separate rule
## ("only U_k can contradict") whose whole job was to excuse ignoring the other
## tiles. With every terminus judged on its own there is nothing to excuse: two
## starts of one side can both be real, and a drop at either is a drop.
##
## Whether the terminus overlaps D does not enter. The test at a terminus is
## against the side's own ref_ex_part, never against the other side, so
## D-membership says nothing about whether the drop is real -- and 84.5% of
## termini lie in neither side's D, so requiring overlap would discard most of
## the data for a claim the walk does not make.
##
## THE THREE SOURCES OF NA MUST NOT COLLAPSE. `!NA | !NA` is NA, not TRUE, and
## each NA means something different:
##   no tile row at all + gene_edge   nothing outward exists  -> POSITIVE
##   tile row, never tested           the join lost it        -> NOT judgeable
##   terminus itself not tested       same                    -> NOT judgeable
## Letting the second become TRUE is the failure this script already warns
## about further down: a lost join reading as a real absence.
## The asymmetry is deliberate. A zero-count TILE is a positive (above): the
## space outward was measured and is empty. A zero-count TERMINUS is NEGATIVE:
## there is no signal there, so it is not a start and cannot anchor anything. It
## fails the floor regardless; the explicit clause makes the asymmetry visible
## rather than an accident of NA handling.
w[, t_ok := isTRUE_v(part_min_read_support) & isTRUE_v(part_called) &
            !isTRUE_v(part_zero)]
w[, u_dead := fcase(
    ## nothing outward exists: the terminus is the gene's outermost part
    !isTRUE_v(tile_present) & isTRUE_v(gene_edge),           TRUE,
    ## measured and EMPTY. exontest emits no row for a unit with no reads, so
    ## this arrives untested -- but a zero mean IS a measurement, and an empty
    ## tile outward of a used start is the strongest boundary evidence there is.
    isTRUE_v(tile_zero),                                     TRUE,
    ## a tile with no support value at all says nothing: a bad join, not a drop
    !isTRUE_v(tile_tested),                                  NA,
    ## measured: dead if under the floor, or tested and uncalled
    default = !isTRUE_v(tile_min_read_support) | !isTRUE_v(tile_called))]

## The terminus-level label, kept beside the boolean so a reader can see WHICH
## kind of drop each positive is. All five are positives; the label records how
## much was known about the tile, which is fixed by its d_support:
##   gene edge                 no tile exists outward
##   outward empty             d_support exactly 0
##   outward near-empty        0 < d_support < 1
##   outward untestable        1 <= d_support < floor; coverage, but untestable
##   outward tested, uncalled  d_support >= floor, tested, no differential
w[, term_verdict := fcase(
    !t_ok & isTRUE_v(part_zero),               "terminus zero count",
    !t_ok & !isTRUE_v(part_tested),            "terminus not tested",
    !t_ok & !isTRUE_v(part_min_read_support),  "terminus below floor",
    !t_ok,                                     "terminus above floor, not called",
    is.na(u_dead),                             "outward tile not tested",
    !u_dead,                                   "signal continues outward",
    !isTRUE_v(tile_present),                   "drop: gene edge",
    isTRUE_v(tile_zero),                       "drop: outward empty",
    !isTRUE_v(tile_min_read_support) & !is.na(tile_ds) & tile_ds < 1,
                                               "drop: outward near-empty",
    !isTRUE_v(tile_min_read_support),          "drop: outward untestable",
    default =                                  "drop: outward tested, uncalled")]
w[, term_positive := t_ok & !is.na(u_dead) & u_dead]

## ---- side level: any() over the termini -----------------------------------
## An existential, so a side with nine termini has nine chances. The per-unit
## p-values are already nested_BH-adjusted within gene, but the side-level OR
## is NOT, so quote the terminus-level figure as the primary one and treat this
## as a convenience rollup.
v <- w[, .(n_term      = .N,
           n_judgeable = sum(t_ok & !is.na(u_dead)),
           n_positive  = sum(term_positive),
           ## LOST means neither tested nor measured. A zero-count unit is
           ## untested by design -- exontest emits no row for it -- but it IS
           ## measured, so it is not lost.
           n_lost      = sum((!isTRUE_v(part_tested) & !isTRUE_v(part_zero)) |
                             (isTRUE_v(tile_present) & !isTRUE_v(tile_tested) &
                              !isTRUE_v(tile_zero))),
           positive    = any(term_positive),
           ## the strongest kind of evidence the side has, since a side with
           ## several termini can carry several kinds at once. Ordered by how
           ## much of the outward space was actually quantified.
           n_tested_uncalled = sum(term_verdict == "drop: outward tested, uncalled"),
           n_empty           = sum(term_verdict == "drop: outward empty"),
           n_near_empty      = sum(term_verdict == "drop: outward near-empty"),
           n_untestable      = sum(term_verdict == "drop: outward untestable"),
           n_gene_edge       = sum(term_verdict == "drop: gene edge"),
           n_continues       = sum(term_verdict == "signal continues outward")),
       by = .(contrast, gene, event, side)]

## ---- side-level verdict ---------------------------------------------------
## Three states, and the denominator matters: a side with no judgeable terminus
## is not a negative, it is a side the question could not be asked of.
v[, side_verdict := fcase(
    n_judgeable == 0, "not judgeable",
    positive,         "positive: a start with nothing outward",
    default =         "negative: signal continues past every start")]
## Among positives, the best evidence present. A side counts under the strongest
## kind it has, so these partition the positives.
## best_evidence is DESCRIPTIVE ONLY. The side verdict is any(term_positive) and
## nothing else; this column just says which kind of evidence the side's
## strongest positive terminus had, so a reader can see how much of the positive
## set is measured. Do not turn it into a gate.
v[, best_evidence := fcase(
    !positive,                              NA_character_,
    n_tested_uncalled > 0 | n_empty > 0,    "measured",
    n_near_empty > 0,                       "near-empty",
    n_untestable > 0,                       "untestable",
    default =                               "gene edge")]

## ---- the assertion that matters -------------------------------------------
## A side that lost its units to a failed join would be scored `no min_read_support
## start` -- indistinguishable from a real lack of depth, and that category is
## one of the reported numbers.
lost <- cl[!unique(v[, .(gene, event, side, contrast)]),
           on = .(gene, event, side, contrast)]
if (nrow(lost))
  stop(sprintf(paste0("FATAL: %d called sides have no terminus unit, e.g. %s\n",
                      "They would be scored 'no start above the read floor', indistinguishable from\n",
                      "low depth. Check the enumeration covers every called side."),
               nrow(lost), paste(head(lost[, paste(gene, event, side)], 5), collapse = ", ")))
cat(sprintf("all %d called sides resolve to terminus units\n", nrow(cl)))

## The TERMINUS table is the primary result -- it is the unit of analysis and
## carries term_verdict, the label that says which kind of drop each positive
## is. The side rollup is an any() over it and goes beside it rather than in
## place of it: writing only the rollup would discard every per-terminus label.
fwrite(w, a$out, sep = "\t")
cat("wrote ", a$out, "  (one row per terminus x contrast)\n", sep = "")
side_out <- sub("\\.tsv$", "", a$out)
side_out <- paste0(side_out, ".sides.tsv")
fwrite(v, side_out, sep = "\t")
cat("wrote ", side_out, "  (side rollup, any() over termini)\n", sep = "")

## ---- summary -------------------------------------------------------------
## The TERMINUS is the unit of analysis, so it is the primary table. The
## side-level rollup is an any() over termini and is NOT multiplicity-adjusted
## across them, so it is reported second and labelled as a rollup.
##
## TWO RATES TRAVEL WITH EVERY FRACTION, because the bounded fraction is not
## comparable across contrasts without them. Measured on the 10-contrast DICE
## lineage panel, over a 17-fold range of background call rate (1.6% to 26.6%):
##
##   t call rate  -> how many sides are judgeable at all (r = +0.83 with
##                   background). 7.2% of termini callable in CD8_vs_CD4
##                   against 25.2% in MONO.CLASSIC_vs_CD8, so 28.7% vs 54.1%
##                   of sides could be examined.
##   background   -> how often ANY measured outward unit is called. This drives
##                   P(U called | t called) (r = +0.83), and the bounded
##                   fraction is essentially its complement (r = -0.77). So a
##                   high bounded fraction in a closely related pair is mostly
##                   a low background, NOT a sharper boundary.
##
## The two move oppositely, so bounded over ALL called sides is nearly flat
## (25-36%, r = +0.28) and is the figure that does not invite the misreading.
## The local effect itself IS constant: stratified by gene, a unit at distance 1
## from a called unit has 2.1-2.7x the odds of being called against one at
## distance 2, across every contrast. Quote that for the effect, and the
## fractions only with both rates beside them.
bg <- res[!is.na(d_support) & d_support >= a$min_reads,
          .(background = mean(called, na.rm = TRUE)), by = contrast]
for (cn in names(CTR)) {
  d <- w[contrast == cn]
  cat("\n=== ", cn, " -- per terminus (n = ", nrow(d), ")\n", sep = "")
  p <- d[, .(termini = .N), by = term_verdict][order(-termini)]
  p[, pct := round(100 * termini / nrow(d), 1)]
  print(p)
  nj <- sum(d$t_ok & !is.na(d$u_dead)); np <- sum(d$term_positive)
  cat(sprintf("judgeable %d of %d (%.1f%%); positive %d (%.1f%% OF JUDGEABLE)\n",
              nj, nrow(d), 100 * nj / nrow(d), np,
              if (nj) 100 * np / nj else NA_real_))
  ## The five positive kinds are already rows of the table printed above, so
  ## there is nothing to summarise further. Report them as they stand: pooling
  ## them into one figure is what made every earlier headline number move when
  ## a threshold changed.
}

## ---- side-level summary ---------------------------------------------------
## An any() over termini, so a side with nine termini has nine chances and this
## is NOT adjusted across them. The terminus table is the primary result.
for (cn in names(CTR)) {
  d <- v[contrast == cn]
  cat("\n=== ", cn, " -- per side (n = ", nrow(d), ")\n", sep = "")
  wt <- w[contrast == cn]
  b_ <- bg[contrast == cn]$background
  cat(sprintf("  background call rate %.1f%% | t call rate %.1f%% (%d of %d termini)\n",
              100 * if (length(b_)) b_ else NA_real_,
              100 * mean(wt$t_ok), sum(wt$t_ok), nrow(wt)))
  p <- d[, .(sides = .N), by = side_verdict][order(-sides)]
  p[, pct := round(100 * sides / nrow(d), 1)]
  print(p)
  nj <- sum(d$n_judgeable > 0)
  cat(sprintf("judgeable %d of %d (%.1f%%); positive %d (%.1f%% OF JUDGEABLE SIDES)\n",
              nj, nrow(d), 100 * nj / nrow(d), sum(d$positive),
              if (nj) 100 * sum(d$positive) / nj else NA_real_))
  ## bounded over ALL called sides: the product of the two rates above, and the
  ## only one of the three that is near-flat across contrasts.
  cat(sprintf("  bounded over ALL called sides %.1f%%  <- quote this one across contrasts\n",
              100 * sum(d$positive) / nrow(d)))
  cat("  positives by their best evidence:\n")
  be <- d[positive == TRUE, .(sides = .N), by = best_evidence][order(-sides)]
  be[, pct_of_positives := round(100 * sides / sum(sides), 1)]
  print(be)
  cat(sprintf("  termini per side: median %g, max %d; positive termini per positive side: median %g, max %d\n",
              median(d$n_term), max(d$n_term),
              median(d[positive == TRUE]$n_positive),
              max(c(0L, d[positive == TRUE]$n_positive))))
}

cat("\n=== termini lost to a failed join, by contrast (NOT drops)\n")
print(w[, .(lost = sum((!isTRUE_v(part_tested) & !isTRUE_v(part_zero)) |
                           (isTRUE_v(tile_present) & !isTRUE_v(tile_tested) &
                            !isTRUE_v(tile_zero)))),
        by = contrast])

cat("\n=== NOTE no pooled figure is reported on purpose: contrasts are",
    "\n    independent tests, and the between-contrast replication rate for the",
    "\n    drop verdict is only about half. Quote the positive count against",
    "\n    JUDGEABLE termini, and say how much of it is measured.\n")
