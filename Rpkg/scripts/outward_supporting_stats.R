#!/usr/bin/env Rscript
## Supporting statistics for the outward-walk write-up (outward_walk_test.md):
## the 5' / 3' power asymmetry, between-contrast replication of the side
## verdict, and the three candidate fractions against each contrast's base
## rates. Recomputed from the verdict tables, so the document is regenerable
## rather than transcribed from a log.
##
## Replication is restricted to sides judgeable in BOTH contrasts of a pair, so
## a non-replication is a disagreement and not an absence. Quoting it over all
## sides instead counts "could not be asked" as "did not replicate" and gives a
## far lower figure that does not mean what it appears to.
##
## Lineage panel only: it has 10 contrasts, and replication needs more than the
## activation panel's two.
##
## Usage: Rscript outward_supporting_stats.R
suppressPackageStartupMessages(library(data.table))
setwd(Sys.getenv("DICE", "/mnt/data1/home/mirahan/DICE"))
## IG=1 as the first argument reads the gene-edge run under lineage/ig instead.
IG  <- commandArgs(TRUE)[1]; if (is.na(IG)) IG <- "0"
DIR <- if (IG == "1") "intronparts/lineage/ig" else "intronparts/lineage"
cat("reading", DIR, "\n")

w  <- fread(file.path(DIR, "outward_walk_verdict.tsv"))
v  <- fread(file.path(DIR, "outward_walk_verdict.sides.tsv"))
st <- unique(fread(file.path(DIR, "outward_side_table.tsv"),
                   select = c("gene","event","side","kind")))

## ---- 1. the arm is a 3' instrument ---------------------------------------
x <- merge(w, st, by = c("gene","event","side"), all.x = TRUE)
cat("\n### 5' vs 3' termini, by contrast\n")
a <- x[!is.na(kind), .(termini   = .N,
                       below     = 100*mean(term_verdict == "terminus below floor"),
                       zero      = 100*mean(term_verdict == "terminus zero count"),
                       judgeable = 100*mean(t_ok & !is.na(u_dead))),
       by = .(contrast, kind)]
b <- dcast(a, contrast ~ kind, value.var = c("termini","below","zero","judgeable"))
b[, gap_below := below_5prime - below_3prime]
print(b[, .(contrast, n5 = termini_5prime, n3 = termini_3prime,
            below5 = round(below_5prime,1), below3 = round(below_3prime,1),
            gap = round(gap_below,1),
            zero5 = round(zero_5prime,1), zero3 = round(zero_3prime,1),
            judg5 = round(judgeable_5prime,1), judg3 = round(judgeable_3prime,1))])
cat(sprintf("\npooled: below floor 5' %.1f%% vs 3' %.1f%%; zero 5' %.1f%% vs 3' %.1f%%\n",
            100*mean(x[kind=="5prime"]$term_verdict == "terminus below floor"),
            100*mean(x[kind=="3prime"]$term_verdict == "terminus below floor"),
            100*mean(x[kind=="5prime"]$term_verdict == "terminus zero count"),
            100*mean(x[kind=="3prime"]$term_verdict == "terminus zero count")))

## ---- 2. between-contrast replication of the side verdict ------------------
## Restricted to sides judgeable in BOTH contrasts of a pair, so a
## non-replication is a disagreement and not an absence.
cat("\n### between-contrast replication (sides judgeable in both)\n")
j <- v[n_judgeable > 0, .(gene, event, side, contrast, positive)]
ctr <- sort(unique(j$contrast))
rep <- rbindlist(lapply(combn(ctr, 2, simplify = FALSE), function(p) {
  m <- merge(j[contrast == p[1]], j[contrast == p[2]],
             by = c("gene","event","side"), suffixes = c(".a",".b"))
  if (!nrow(m)) return(NULL)
  data.table(a = p[1], b = p[2], n = nrow(m),
             both_pos = sum(m$positive.a & m$positive.b),
             pos_a    = sum(m$positive.a),
             pos_b    = sum(m$positive.b),
             both_neg = sum(!m$positive.a & !m$positive.b),
             neg_a    = sum(!m$positive.a))
}))
rep[, repl_pos := 100 * both_pos / pos_a]
rep[, repl_neg := 100 * both_neg / neg_a]
rep[, agree    := 100 * (both_pos + both_neg) / n]
print(rep[order(-n), .(a, b, n, repl_pos = round(repl_pos,1),
                       repl_neg = round(repl_neg,1), agree = round(agree,1))])
cat(sprintf("\nmedian over %d pairs: positive replication %.1f%%, negative %.1f%%, agreement %.1f%%\n",
            nrow(rep), median(rep$repl_pos), median(rep$repl_neg), median(rep$agree)))

## ---- 3. correlation of the three fractions with the base rates ------------
cat("\n### the three fractions against background and t call rate\n")
bgt <- w[, .(t_rate = 100*mean(t_ok)), by = contrast]
sm  <- v[, .(sides = .N, judgeable = sum(n_judgeable > 0), positive = sum(positive)),
         by = contrast]
sm[, pct_of_judgeable := 100*positive/judgeable]
sm[, pct_of_all       := 100*positive/sides]
sm <- merge(sm, bgt, by = "contrast")
print(sm[order(t_rate), .(contrast, sides, t_rate = round(t_rate,1),
                          pct_of_judgeable = round(pct_of_judgeable,1),
                          pct_of_all = round(pct_of_all,1))])
cat(sprintf("\nr(t_rate, pct_of_judgeable) = %+.2f   r(t_rate, pct_of_all) = %+.2f\n",
            cor(sm$t_rate, sm$pct_of_judgeable), cor(sm$t_rate, sm$pct_of_all)))
cat(sprintf("range: of judgeable %.1f-%.1f (%.1f-fold spread in t rate)\n",
            min(sm$pct_of_judgeable), max(sm$pct_of_judgeable),
            max(sm$t_rate)/min(sm$t_rate)))
cat(sprintf("       of all called %.1f-%.1f\n", min(sm$pct_of_all), max(sm$pct_of_all)))
