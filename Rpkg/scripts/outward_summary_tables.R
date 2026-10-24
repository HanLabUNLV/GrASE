#!/usr/bin/env Rscript
## Emit the summary tables quoted in the outward-walk write-up, so the document
## is regenerable rather than transcribed from a log. One row per contrast.
##
## Usage: Rscript outward_summary_tables.R [PANEL] [IG]   (lineage | activation)
## IG=1 reads the gene-edge run under <panel dir>/ig instead. That run has no
## `gene edge` best_evidence column at all -- every gene-edge terminus now has a
## flank -- so the measured share is taken over whatever evidence columns exist
## rather than over a fixed four.
suppressPackageStartupMessages(library(data.table))
setwd(Sys.getenv("DICE", "/mnt/data1/home/mirahan/DICE"))
PANEL <- commandArgs(TRUE)[1]; if (is.na(PANEL)) PANEL <- "lineage"
IG    <- commandArgs(TRUE)[2]; if (is.na(IG)) IG <- "0"
dir <- if (PANEL == "lineage") "intronparts/lineage" else "intronparts"
if (IG == "1") dir <- file.path(dir, "ig")

w <- fread(file.path(dir, "outward_walk_verdict.tsv"))
v <- fread(file.path(dir, "outward_walk_verdict.sides.tsv"))

tr <- w[, .(termini = .N, t_called = sum(t_ok), t_rate = 100*mean(t_ok)), by = contrast]
sd <- v[, .(sides = .N, judgeable = sum(n_judgeable > 0), positive = sum(positive),
            pct_of_judgeable = 100*sum(positive)/sum(n_judgeable > 0),
            pct_of_all       = 100*sum(positive)/.N), by = contrast]
be <- dcast(v[positive == TRUE, .N, by = .(contrast, best_evidence)],
            contrast ~ best_evidence, value.var = "N", fill = 0)
ev <- setdiff(names(be), "contrast")
be[, measured_pct := 100 * measured / rowSums(.SD), .SDcols = ev]

out <- Reduce(function(a,b) merge(a,b,by="contrast"), list(tr, sd, be[, .(contrast, measured_pct)]))
setorder(out, t_rate)
for (j in c("t_rate","pct_of_judgeable","pct_of_all","measured_pct"))
  set(out, j = j, value = round(out[[j]], 1))
cat("\n### panel:", PANEL, if (IG == "1") "(with gene-edge flanks)" else "", "\n")
cat("    evidence columns:", paste(ev, collapse = ", "), "\n")
print(out[, .(contrast, sides, termini, t_rate, judgeable, positive,
              pct_of_judgeable, pct_of_all, measured_pct)])

f <- file.path(dir, "outward_walk_summary.tsv")
fwrite(out, f, sep = "\t")
cat("\nwrote", f, "\n")
