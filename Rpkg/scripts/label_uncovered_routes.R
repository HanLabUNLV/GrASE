#!/usr/bin/env Rscript
## Compute intron_uncovered1/2 for junction-substituted bipartition sides.
##
## A side whose exonic distinct set is empty is measured on its junctions
## instead. Under the cut rule each ROUTE through the bubble contributes one
## junction; a route with no distinct intron contributes nothing, so D
## undercounts those molecules. `uncovered` counts such routes -- a side with
## uncovered > 0 is only partly represented by its junction measure.
##
## Only substituted sides are relabelled: `uncovered` is undefined elsewhere.
## No statistics are recomputed; this emits a per-side flag for filtering.
##
## Usage:
##   Rscript scripts/label_uncovered_routes.R <results.annotated.txt> \
##           <graphml_dir> <out.tsv> [rpkg_R_dir]

suppressMessages({library(igraph); library(dplyr)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: <results> <graphml_dir> <out.tsv> [R_dir]")
res_path <- args[1]; graphml_dir <- args[2]; out_path <- args[3]
rdir <- if (length(args) >= 4) args[4] else
        file.path(dirname(dirname(normalizePath(sys.frames()[[1]]$ofile %||% "."))), "R")
if (!dir.exists(rdir)) rdir <- "/mnt/data1/home/mirahan/GrASE/Rpkg/R"
for (f in list.files(rdir, pattern = "\\.R$", full.names = TRUE))
  try(source(f), silent = TRUE)
stopifnot(exists("precompute_gene_graph"), exists("route_node_pairs_by_route"))

d <- read.delim(res_path, stringsAsFactors = FALSE)
if (!"comparison" %in% names(d)) d$comparison <- "diff1_vs_ref"
d$sd <- ifelse(grepl("^diff1", d$comparison), d$setdiff1, d$setdiff2)
sub <- d[is.na(d$sd) | d$sd %in% c("", "NA"), ]
sub <- sub[!duplicated(sub[, c("gene", "event", "comparison")]), ]
cat(sprintf("substituted sides to relabel: %d  (%d genes)\n",
            nrow(sub), length(unique(sub$gene))))

sub$intron_uncovered <- NA_integer_
sub$n_routes <- NA_integer_
sub$n_junctions_cut <- NA_integer_
sub$n_junctions_union <- NA_integer_

ge_cache <- new.env(parent = emptyenv())
get_ge <- function(gene) {
  if (!is.null(ge_cache[[gene]])) return(ge_cache[[gene]])
  gp <- file.path(graphml_dir, paste0(gene, ".graphml"))
  if (!file.exists(gp)) { ge_cache[[gene]] <- NA; return(NA) }
  g <- tryCatch(read_graph(gp, format = "graphml"), error = function(e) NULL)
  ge <- if (is.null(g)) NA else tryCatch(precompute_gene_graph(g),
                                         error = function(e) NA)
  ge_cache[[gene]] <- ge
  ge
}

n_ok <- 0L
for (i in seq_len(nrow(sub))) {
  ge <- get_ge(sub$gene[i])
  if (length(ge) == 1L && is.na(ge)) next
  p1 <- sub$path1[i]; p2 <- sub$path2[i]
  if (is.na(p1) || is.na(p2)) next
  bv <- unique(suppressWarnings(as.integer(unlist(
         strsplit(paste(p1, p2, sep = ","), "[-,]")))))
  bv <- bv[!is.na(bv)]
  t1 <- trimws(unlist(strsplit(sub$transcripts1[i], ",")))
  t2 <- trimws(unlist(strsplit(sub$transcripts2[i], ",")))
  base <- list(ge, t1, t2, "chr1", bv,
               route_node_pairs(p1), route_node_pairs(p2))
  k <- tryCatch(do.call(label_bipartition_introns,
        c(base, list(routes1 = route_node_pairs_by_route(p1),
                     routes2 = route_node_pairs_by_route(p2), rule = "cut"))),
        error = function(e) NULL)
  u <- tryCatch(do.call(label_bipartition_introns, c(base, list(rule = "union"))),
        error = function(e) NULL)
  if (is.null(k) || is.null(u)) next
  side <- if (grepl("^diff1", sub$comparison[i])) 1L else 2L
  kk <- if (side == 1L) k$distinct1 else k$distinct2
  uu <- if (side == 1L) u$distinct1 else u$distinct2
  cu <- if (side == 1L) k$uncovered1 else k$uncovered2
  rr <- length(route_node_pairs_by_route(if (side == 1L) p1 else p2))
  sub$intron_uncovered[i]   <- cu
  sub$n_routes[i]           <- rr
  sub$n_junctions_cut[i]    <- if (is.na(kk)) 0L else length(strsplit(kk, ",")[[1]])
  sub$n_junctions_union[i]  <- if (is.na(uu)) 0L else length(strsplit(uu, ",")[[1]])
  n_ok <- n_ok + 1L
  if (n_ok %% 500L == 0L) cat("  ", n_ok, "relabelled\n")
}

keep <- c("gene", "event", "comparison", "intron_uncovered", "n_routes",
          "n_junctions_cut", "n_junctions_union")
write.table(sub[, keep], out_path, sep = "\t", quote = FALSE, row.names = FALSE)

ok <- !is.na(sub$intron_uncovered)
cat(sprintf("\nrelabelled          : %d of %d\n", sum(ok), nrow(sub)))
if (any(ok)) {
  cat(sprintf("uncovered == 0      : %d (%.1f%%)  -> keep\n",
              sum(sub$intron_uncovered[ok] == 0), 100 * mean(sub$intron_uncovered[ok] == 0)))
  cat(sprintf("uncovered  > 0      : %d (%.1f%%)  -> drop\n",
              sum(sub$intron_uncovered[ok] > 0), 100 * mean(sub$intron_uncovered[ok] > 0)))
  cat(sprintf("junctions summed    : union median %.0f -> cut median %.0f\n",
              median(sub$n_junctions_union[ok]), median(sub$n_junctions_cut[ok])))
}
cat("wrote", out_path, "\n")
