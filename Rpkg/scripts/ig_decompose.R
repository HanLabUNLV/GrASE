## Three-way decomposition of what changed in the outward walk between the
## pre-flank run and the gene-edge run.
##
## A  old scorer, old run, no flanks   intronparts/*.oldscorer
## B  CURRENT scorer, old run, no flanks   intronparts/outward_walk_verdict.tsv
## C  CURRENT scorer, new run, WITH flanks intronparts/ig/outward_walk_verdict.tsv
##
## B - A isolates the scorer change (the u_dead reorder: a measured-thin flank is
## a positive). C - B isolates the flanks plus the re-run's effect on padj.
## Reporting C - A alone would credit the flanks with all three.
##
## Usage: Rscript ig_decompose.R
suppressPackageStartupMessages(library(data.table))
setwd(Sys.getenv("DICE", "/mnt/data1/home/mirahan/DICE"))
K <- c("gene","event","side","contrast","rank")
rd <- function(f) fread(f, sep="\t", showProgress=FALSE,
                        colClasses=list(character=c("gene","event")))

A <- rd("intronparts/outward_walk_verdict.tsv.oldscorer")
B <- rd("intronparts/outward_walk_verdict.tsv")
C <- rd("intronparts/ig/outward_walk_verdict.tsv")

cat("=== terminus counts\n")
cat(sprintf("  A %d   B %d   C %d\n\n", nrow(A), nrow(B), nrow(C)))

tab <- function(d, lab) {
  cat("===", lab, "-- terminus verdicts by contrast\n")
  p <- dcast(d[, .N, by=.(term_verdict, contrast)], term_verdict ~ contrast, value.var="N", fill=0L)
  print(p)
  cat("\n")
  for (cn in unique(d$contrast)) {
    x <- d[contrast==cn]
    nj <- sum(x$t_ok & !is.na(x$u_dead)); np <- sum(x$term_positive)
    cat(sprintf("  %-20s judgeable %5d of %5d (%.1f%%)  positive %5d (%.1f%% of judgeable)\n",
                cn, nj, nrow(x), 100*nj/nrow(x), np, 100*np/nj))
  }
  cat("\n")
}
tab(A, "A old scorer / old run / no flanks")
tab(B, "B current scorer / old run / no flanks")
tab(C, "C current scorer / new run / WITH flanks")

mv <- function(X, Y, labx, laby) {
  m <- merge(X[, c(K,"term_verdict","term_positive"), with=FALSE],
             Y[, c(K,"term_verdict","term_positive"), with=FALSE],
             by=K, suffixes=c("_x","_y"), all=TRUE)
  cat(sprintf("=== %s -> %s : verdict transitions (only rows that moved)\n", labx, laby))
  d <- m[is.na(term_verdict_x) | is.na(term_verdict_y) | term_verdict_x != term_verdict_y]
  print(d[, .(termini=.N), by=.(term_verdict_x, term_verdict_y)][order(-termini)])
  cat(sprintf("  moved %d of %d (%.1f%%);  positives %d -> %d\n\n",
              nrow(d), nrow(m), 100*nrow(d)/nrow(m),
              sum(m$term_positive_x, na.rm=TRUE), sum(m$term_positive_y, na.rm=TRUE)))
  invisible(m)
}
mv(A, B, "A", "B")
m_bc <- mv(B, C, "B", "C")

cat("=== the gene-edge termini specifically (B's 'drop: gene edge' -> C)\n")
ge <- m_bc[term_verdict_x == "drop: gene edge"]
print(ge[, .(termini=.N), by=.(contrast, term_verdict_y)][order(contrast, -termini)])
cat(sprintf("\n  %d gene-edge termini; still positive %d (%.1f%%); now NEGATIVE %d (%.1f%%); NA %d\n\n",
            nrow(ge), sum(ge$term_positive_y, na.rm=TRUE),
            100*sum(ge$term_positive_y, na.rm=TRUE)/nrow(ge),
            sum(ge$term_verdict_y == "signal continues outward"),
            100*sum(ge$term_verdict_y == "signal continues outward")/nrow(ge),
            sum(ge$term_verdict_y == "outward flank unmeasurable")))

cat("=== side level, A / B / C\n")
sd <- function(f) fread(f, sep="\t", showProgress=FALSE,
                        colClasses=list(character=c("gene","event")))
S <- list(A=sd("intronparts/outward_walk_verdict.sides.tsv.oldscorer"),
          B=sd("intronparts/outward_walk_verdict.sides.tsv"),
          C=sd("intronparts/ig/outward_walk_verdict.sides.tsv"))
for (cn in unique(S$C$contrast)) {
  cat(sprintf("  %s\n", cn))
  for (nm in names(S)) {
    x <- S[[nm]][contrast==cn]
    cat(sprintf("    %s  sides %5d  judgeable %.1f%%  positive/judgeable %.1f%%  positive/all %.1f%%\n",
                nm, nrow(x), 100*mean(x$n_judgeable>0),
                100*sum(x$positive)/sum(x$n_judgeable>0), 100*sum(x$positive)/nrow(x)))
  }
}

cat("\n=== LINEAGE, same three-way decomposition\n")
LA <- rd("intronparts/lineage/outward_walk_verdict.tsv.oldscorer")
L0 <- rd("intronparts/lineage/outward_walk_verdict.tsv")       # B: current scorer, no flanks
L1 <- rd("intronparts/lineage/ig/outward_walk_verdict.tsv")    # C: current scorer, flanks
mAB <- merge(LA[, c(K,"term_verdict","term_positive"), with=FALSE],
             L0[, c(K,"term_verdict","term_positive"), with=FALSE],
             by=K, suffixes=c("_x","_y"), all=TRUE)
dAB <- mAB[term_verdict_x != term_verdict_y]
cat(sprintf("  A->B (scorer only): moved %d of %d (%.1f%%); positives %d -> %d\n",
            nrow(dAB), nrow(mAB), 100*nrow(dAB)/nrow(mAB),
            sum(mAB$term_positive_x, na.rm=TRUE), sum(mAB$term_positive_y, na.rm=TRUE)))
print(dAB[, .(termini=.N), by=.(term_verdict_x, term_verdict_y)][order(-termini)])
cat(sprintf("  termini: %d (all three)\n", nrow(L1)))
mL <- merge(L0[, c(K,"term_verdict","term_positive"), with=FALSE],
            L1[, c(K,"term_verdict","term_positive"), with=FALSE],
            by=K, suffixes=c("_x","_y"), all=TRUE)
geL <- mL[term_verdict_x == "drop: gene edge"]
cat(sprintf("  gene-edge termini %d; still positive %d (%.1f%%); now NEGATIVE %d (%.1f%%); NA %d\n",
            nrow(geL), sum(geL$term_positive_y, na.rm=TRUE),
            100*sum(geL$term_positive_y, na.rm=TRUE)/nrow(geL),
            sum(geL$term_verdict_y == "signal continues outward", na.rm=TRUE),
            100*sum(geL$term_verdict_y == "signal continues outward", na.rm=TRUE)/nrow(geL),
            sum(geL$term_verdict_y == "outward flank unmeasurable", na.rm=TRUE)))
print(geL[, .(termini=.N), by=.(term_verdict_y)][order(-termini)])
SL0 <- sd("intronparts/lineage/outward_walk_verdict.sides.tsv")
SL1 <- sd("intronparts/lineage/ig/outward_walk_verdict.sides.tsv")
cat("\n  side level by contrast (B current scorer no flanks -> C with flanks)\n")
for (cn in sort(unique(SL1$contrast))) {
  a <- SL0[contrast==cn]; b <- SL1[contrast==cn]
  cat(sprintf("    %-24s judgeable %.1f%%->%.1f%%  pos/judg %.1f%%->%.1f%%  positive/all %.1f%%->%.1f%%\n",
    cn, 100*mean(a$n_judgeable>0), 100*mean(b$n_judgeable>0),
    100*sum(a$positive)/sum(a$n_judgeable>0), 100*sum(b$positive)/sum(b$n_judgeable>0),
    100*sum(a$positive)/nrow(a), 100*sum(b$positive)/nrow(b)))
}
cat("\n=== unmeasurable flanks, both panels\n")
for (nm in c("activation","lineage")) {
  d <- if (nm=="activation") C else L1
  cat(sprintf("  %-11s term_verdict == 'outward flank unmeasurable': %d of %d\n",
              nm, sum(d$term_verdict=="outward flank unmeasurable"), nrow(d)))
}
