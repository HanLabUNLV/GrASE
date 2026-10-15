#!/usr/bin/env Rscript
## Attach the outward walk's verdict to the external exontest results.
##
## Splits a called external event's sides on whether the differential effect is
## CONFINED to an annotated transcription terminus:
##
##   bounded.5prime     effect at a 5' terminus, absent immediately outward
##   bounded.3prime     the same at a 3' terminus
##   unbounded.5prime   effect present at the terminus AND still outward of
##   unbounded.3prime   every one of them
##   untested.5prime    no terminus could be judged: below the read floor, or
##   untested.3prime    uncalled, or its outward tile never tested
##   NA                 the row is not a significant external call, so the walk
##                      never looked at it
##
## WHY NOT "TSS" AND "NON-TSS". filterExternal.R puts both source == "R" and
## sink == "L" events in the external class, so 71% of bounded sides are 3'
## (5,693 against 2,367 on the 2026-10-14 activation run). Calling the bounded
## set "TSS" would mislabel seven in ten. 5prime/3prime carries the distinction
## instead of hiding it.
##
## WHY THREE STATES AND NOT TWO. `untested` is 48.5% of significant external
## rows. It is not evidence either way, and folding it into either class would
## put 15,009 rows on the wrong side of a claim. A two-way split is only
## defensible after dropping it, with the drop stated.
##
## `unbounded` describes the observation and asserts no mechanism: an effect
## extending past every annotated terminus is consistent with read-through, with
## an unannotated terminus further out, and with the event not being about
## initiation or termination at all. Do not relabel it "read-through".
##
## The verdict is any() over the side's termini -- a side with several termini is
## bounded if ANY of them has a dead tile outward. See score_outward_walk.R.
##
## OUTPUT IS ONE FILE PER CLASS: <prefix>.<class>.txt, following the pipeline's
## convention of a file per subset rather than a column to filter on.
##
## By default each file holds only that class's rows, so the six do NOT sum to
## the input: the 97.4% of rows that are not significant external calls belong
## to no class and are dropped. That is fine for anything keyed on the called
## set, and WRONG for GO enrichment, which needs the untested universe as its
## background. Pass --keep_universe to put every non-significant row into every
## class file as well -- correct for enrichment, at six times the size.
##
## Usage:
##   Rscript scripts/annotate_bounded.R \
##     --results <external ...annotated.txt> \
##     --sides   intronparts/outward_walk_verdict.sides.tsv \
##     --side_table intronparts/outward_side_table.tsv \
##     --out_prefix <...annotated>  [--keep_universe]
suppressPackageStartupMessages({library(data.table); library(optparse)})

a <- parse_args(OptionParser(option_list = list(
  make_option("--results",    type = "character"),
  make_option("--sides",      type = "character"),
  make_option("--side_table", type = "character"),
  make_option("--out_prefix", type = "character",
              help = "writes <prefix>.<class>.txt for each class"),
  make_option("--keep_universe", action = "store_true", default = FALSE,
              help = paste("also write every non-significant row into each",
                           "class file, as a shared background. Needed for GO",
                           "enrichment; six times the output size.")))))
for (r in c("results", "sides", "side_table", "out_prefix"))
  if (is.null(a[[r]])) stop("--", r, " is required", call. = FALSE)

res <- fread(a$results, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene", "event")))
sv  <- fread(a$sides, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene", "event")))
st  <- fread(a$side_table, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene", "event")))

## kind is constant per side and per event (verified: 0 of 20,868 sides mixed),
## so one lookup per side is well defined.
kd <- unique(st[, .(gene, event, side, kind)])
dup <- kd[, .N, by = .(gene, event, side)][N > 1]
if (nrow(dup)) stop(nrow(dup), " sides carry more than one `kind`; the label ",
                    "would be ambiguous", call. = FALSE)

sv <- merge(sv, kd, by = c("gene", "event", "side"), all.x = TRUE)
if (anyNA(sv$kind)) stop("side verdict rows with no `kind` in --side_table")
## `kind` already carries 5prime / 3prime (enumerate_outward_units.R), so the
## label is the value itself rather than a translation.
if (!all(sv$kind %in% c("5prime", "3prime")))
  stop("unexpected `kind` values: ",
       paste(setdiff(unique(sv$kind), c("5prime", "3prime")), collapse = ", "),
       "\n  a side table built before the 5prime/3prime rename must be ",
       "regenerated.", call. = FALSE)
sv[, endlab := kind]
sv[, state := fcase(side_verdict %like% "^positive", "bounded",
                    side_verdict %like% "^negative", "unbounded",
                    default = "untested")]
sv[, bounded := paste0(state, ".", endlab)]

## exontest writes one row per side as `comparison`; side 1 is diff1.
res[, side := fifelse(comparison == "diff1_vs_ref", 1L, 2L)]
n0 <- nrow(res)
res <- merge(res, sv[, .(gene, event, side, contrast, bounded)],
             by = c("gene", "event", "side", "contrast"), all.x = TRUE)
if (nrow(res) != n0)
  stop("the join changed the row count (", n0, " -> ", nrow(res),
       "); the side verdict is not unique per (gene, event, side, contrast)")

## Every SIGNIFICANT external row must have a verdict: the walk was run on
## exactly that set, so a gap means the two inputs are from different runs.
miss <- res[significant == TRUE & is.na(bounded)]
if (nrow(miss))
  stop(nrow(miss), " significant rows have no verdict, e.g. ",
       paste(head(miss[, paste(gene, event, comparison, contrast)], 3),
             collapse = "; "),
       "\n  --results and --sides look like different runs.", call. = FALSE)

res[, side := NULL]
cat(sprintf("rows %d; with a verdict %d (%.1f%%)\n", nrow(res),
            sum(!is.na(res$bounded)), 100 * mean(!is.na(res$bounded))))
print(res[!is.na(bounded), .(rows = .N), by = .(contrast, bounded)][order(contrast, -rows)])

## ---- one file per class ---------------------------------------------------
bg <- if (isTRUE(a$keep_universe)) res[is.na(bounded)] else NULL
if (!is.null(bg))
  cat(sprintf("--keep_universe: adding %d non-significant rows to every class file\n",
              nrow(bg)))
classes <- sort(unique(res$bounded[!is.na(res$bounded)]))
for (cl in classes) {
  d <- res[bounded == cl]
  if (!is.null(bg)) d <- rbind(d, bg)
  fn <- sprintf("%s.%s.txt", a$out_prefix, cl)
  fwrite(d, fn, sep = "\t", quote = FALSE)
  cat(sprintf("wrote %-78s %7d rows\n", fn, nrow(d)))
}
nd <- sum(is.na(res$bounded))
if (!isTRUE(a$keep_universe) && nd)
  cat(sprintf("NOT written: %d rows in no class (not a significant external call).\n",
              nd),
      "  The class files therefore do not sum to the input. Pass --keep_universe\n",
      "  if a downstream step needs the full background.\n", sep = "")
