#!/usr/bin/env Rscript
## Four-panel figure for one GrASE bipartition: gene model, every test at the
## locus, the tested bipartition as graph routes, and the path proportion
## across conditions.
##
## R port of scripts/plot_test_results.py, kept deliberately close to it so the
## two can be diffed: same panel order, same colours, same ordering rules, same
## call rule (the `significant` column exontest.R wrote, never re-derived here).
##
## Usage: Rscript scripts/plot_test_results.R --gene ENSG... --event 11 \
##          --kind TSSTTS --analysis <dir> --conditions A,B,C --out fig.pdf
suppressPackageStartupMessages({ library(igraph); library(optparse) })

CB <- c(S = "#4477AA", D1 = "#CCBB44", D2 = "#EE6677", other = "#DDDDDD")

opt_list <- list(
  make_option("--gene", type = "character"),
  make_option("--event", type = "character"),
  make_option("--analysis", type = "character"),
  make_option("--kind", type = "character", default = "TSSTTS"),
  make_option("--conditions", type = "character"),
  make_option("--comparison", type = "character", default = NULL),
  make_option("--gffdir", type = "character", default = "~/DICE/dexseq.gff"),
  make_option("--bipdir", type = "character", default = "~/DICE/bipartition.filtered"),
  make_option("--graphmldir", type = "character", default = "~/DICE/graphml.v34"),
  make_option("--sjcounts", type = "character", default = ""),
  make_option("--stack_kinds", type = "character", default = "internal,TSSTTS"),
  make_option("--max_stack", type = "integer", default = 0L),
  make_option("--feature_contrast", type = "character", default = NULL),
  make_option("--mark_parts", type = "character", default = ""),
  make_option("--symbol", type = "character", default = NULL),
  make_option("--note", type = "character", default = ""),
  make_option("--desc", type = "character", default = NULL),
  make_option("--ylab", type = "character", default = NULL),
  make_option("--raw", action = "store_true", default = FALSE),
  make_option("--out", type = "character")
)
a <- parse_args(OptionParser(option_list = opt_list))
for (req in c("gene", "event", "analysis", "conditions", "out"))
  if (is.null(a[[req]])) stop("--", req, " is required")

## ---------------------------------------------------------------- helpers
parts <- function(cell) {
  if (is.null(cell) || is.na(cell) || cell %in% c("", "NA", "NaN")) return(character(0))
  x <- trimws(strsplit(cell, ",")[[1]]); x[nzchar(x) & x != "NA"]
}
set_length <- function(cell, coords) {
  p <- parts(cell); p <- p[p %in% names(coords)]
  if (!length(p)) return(NA_real_)
  sum(vapply(p, function(e) coords[[e]][2] - coords[[e]][1] + 1, numeric(1)))
}
pi_perbase <- function(pi, len_d, len_s) {
  den <- pi * len_s + (1 - pi) * len_d
  ifelse(den > 0, pi * len_s / den, NA_real_)
}
split_junctions <- function(cell) {
  if (is.null(cell) || is.na(cell) || !nzchar(trimws(cell)) || trimws(cell) == "NA")
    return(list())
  out <- list()
  for (j in strsplit(cell, ",")[[1]]) {
    f <- strsplit(trimws(j), ":")[[1]]
    if (length(f) < 3) next
    v <- suppressWarnings(as.integer(f[(length(f) - 1):length(f)]))
    if (!any(is.na(v))) out[[length(out) + 1]] <- sort(v)
  }
  out
}
read_tsv <- function(p) read.table(p, header = TRUE, sep = "\t", quote = "",
                                   stringsAsFactors = FALSE, comment.char = "",
                                   colClasses = "character")

## graph: node name -> genomic position, and the COARSE edge type per node pair.
## An ex_part edge must never overwrite a coarse ex/in/R/L edge, because vertex
## path hops are always coarse.
read_graph_attrs <- function(path) {
  g <- igraph::read_graph(path, format = "graphml")
  nm <- igraph::V(g)$name; ps <- igraph::V(g)$position
  pos <- setNames(as.list(ps), nm)
  el <- igraph::as_edgelist(g, names = TRUE)
  et <- igraph::E(g)$ex_or_in
  keep <- et != "ex_part"
  edge_type <- setNames(as.list(et[keep]), paste(el[keep, 1], el[keep, 2], sep = "\r"))
  list(pos = pos, edge_type = edge_type)
}
etype <- function(ET, u, v) { k <- paste(u, v, sep = "\r"); if (!is.null(ET[[k]])) ET[[k]] else NA_character_ }

gffdir <- path.expand(a$gffdir); bipdir <- path.expand(a$bipdir)
graphmldir <- path.expand(a$graphmldir); andir <- path.expand(a$analysis)
base_id <- sub("\\..*$", "", a$gene)

hits <- sort(Sys.glob(file.path(gffdir, paste0(base_id, "*.dexseq.gff"))))
if (!length(hits)) stop("no dexseq gff for ", base_id, " in ", gffdir)
gff <- hits[1]; gv <- sub("\\.dexseq.*$", "", basename(gff))

## exonic part coordinates and strand
gl <- readLines(gff); gl <- gl[grepl("\texonic_part\t", gl)]
f <- strsplit(gl, "\t")
pn <- vapply(f, function(x) sub('".*$', "", sub('^.*exonic_part_number "', "", x[9])), "")
coords <- setNames(lapply(f, function(x) as.integer(c(x[4], x[5]))), paste0("E", pn))
strand <- f[[1]][7]
names_e <- names(coords)[order(as.integer(sub("^E", "", names(coords))))]
idx <- setNames(seq_along(names_e) - 1L, names_e)
conds <- trimws(strsplit(a$conditions, ",")[[1]])

G <- list(pos = list(), edge_type = list())
gpath <- file.path(graphmldir, paste0(gv, ".graphml"))
if (file.exists(gpath)) G <- read_graph_attrs(gpath)

boundary_x <- new.env(parent = emptyenv())
for (e in names_e) {
  assign(as.character(coords[[e]][1]), idx[[e]], envir = boundary_x)
  assign(as.character(coords[[e]][2]), idx[[e]] + 1, envir = boundary_x)
}
bx_pos <- sort(as.numeric(ls(boundary_x)))
bx_val <- vapply(as.character(bx_pos), function(k) get(k, envir = boundary_x), numeric(1))

node_x <- function(nm) {
  if (identical(nm, "R")) return(if (strand == "-") length(names_e) + 0.4 else -0.4)
  if (identical(nm, "L")) return(if (strand == "-") -0.4 else length(names_e) + 0.4)
  v <- suppressWarnings(as.integer(G$pos[[nm]]))
  if (is.na(v)) return(NA_real_)
  k <- as.character(v); if (exists(k, envir = boundary_x)) get(k, envir = boundary_x) else NA_real_
}
gx <- function(p) {                      # interpolate onto the exonic-part scale
  if (!length(bx_pos)) return(NA_real_)
  if (p <= bx_pos[1]) return(bx_val[1])
  if (p >= bx_pos[length(bx_pos)]) return(bx_val[length(bx_val)])
  i <- findInterval(p, bx_pos)
  g1 <- bx_pos[i]; g2 <- bx_pos[i + 1]; x1 <- bx_val[i]; x2 <- bx_val[i + 1]
  if (g2 == g1) x1 else x1 + (x2 - x1) * (p - g1) / (g2 - g1)
}
route_junctions <- function(nodes) {
  out <- list()
  if (length(nodes) < 2) return(out)
  for (i in seq_len(length(nodes) - 1)) {
    if (!identical(etype(G$edge_type, nodes[i], nodes[i + 1]), "in")) next
    pf <- suppressWarnings(as.integer(G$pos[[nodes[i]]]))
    pt <- suppressWarnings(as.integer(G$pos[[nodes[i + 1]]]))
    if (is.na(pf) || is.na(pt)) next
    out[[length(out) + 1]] <- sort(c(pf, pt))
  }
  out
}
hop_parts <- function(u, v) {
  if (!identical(etype(G$edge_type, u, v), "ex")) return(character(0))
  pf <- suppressWarnings(as.integer(G$pos[[u]])); pt <- suppressWarnings(as.integer(G$pos[[v]]))
  if (is.na(pf) || is.na(pt)) return(character(0))
  lo <- min(pf, pt); hi <- max(pf, pt)
  names_e[vapply(names_e, function(e) coords[[e]][1] >= lo && coords[[e]][2] <= hi, logical(1))]
}

## merged exon+SJ metadata: intron_distinct1/2 and intron_shared per event
sj_meta <- list()
if (nzchar(a$sjcounts)) {
  d <- path.expand(a$sjcounts)
  cand <- c(file.path(d, paste0(gv, ".bipartition.exoncnt.txt")),
            sort(Sys.glob(file.path(d, paste0(base_id, ".*.bipartition.exoncnt.txt")))))
  for (fp in cand) {
    if (!file.exists(fp)) next
    tb <- read_tsv(fp)
    if (!"intron_distinct1" %in% names(tb)) break
    keep <- !duplicated(tb$event)
    sj_meta <- setNames(lapply(which(keep), function(i)
      as.list(tb[i, c("intron_distinct1", "intron_distinct2", "intron_shared")])),
      tb$event[keep])
    break
  }
}

read_gene_rows <- function(kind) {
  for (p in c(file.path(andir, sprintf("test_bipartition.%s_betabinom_EBmap.annotated.txt", kind)),
              file.path(andir, "exontest_results",
                        sprintf("test_bipartition.%s_betabinom_EBmap.annotated.txt", kind)))) {
    if (!file.exists(p)) next
    tb <- read_tsv(p)
    return(tb[sub("\\..*$", "", tb$gene) == base_id, , drop = FALSE])
  }
  NULL
}
is_sig <- function(v) {
  if (any(is.na(v)) || !length(v)) stop("no `significant` column: re-run exontest.R")
  v == "TRUE"
}

stack_kinds <- trimws(strsplit(a$stack_kinds, ",")[[1]])
all_rows <- setNames(lapply(stack_kinds, read_gene_rows), stack_kinds)
if (!a$kind %in% names(all_rows)) all_rows[[a$kind]] <- read_gene_rows(a$kind)
rows <- all_rows[[a$kind]]
rows <- rows[rows$event == as.character(a$event), , drop = FALSE]
if (!nrow(rows)) stop("no tests for ", base_id, " event ", a$event, " (", a$kind, ")")

num <- function(x) suppressWarnings(as.numeric(x))
if (!is.null(a$comparison)) {
  comp <- a$comparison
} else {
  eff <- ifelse(!is.na(num(rows$padj)) & num(rows$padj) < 0.05, abs(num(rows$delta_pi)), -1)
  eff[is.na(eff)] <- -1
  comp <- if (max(eff) > 0) rows$comparison[which.max(eff)] else {
    pj <- num(rows$padj); pj[is.na(pj)] <- 1
    rows$comparison[which.min(pj)]
  }
}
rows <- rows[rows$comparison == comp, , drop = FALSE]
if (!nrow(rows)) stop("no rows for comparison ", comp)

focal <- rows[1, ]
S <- focal$ref_ex_part; D1 <- focal$setdiff1; D2 <- focal$setdiff2
which_D <- if (grepl("diff1", comp)) "D1" else "D2"
tested_D <- if (which_D == "D1") D1 else D2

pivals <- list(); best <- list()
for (i in seq_len(nrow(rows))) {
  cc <- rows$contrast[i]
  if (!grepl("_vs_", cc)) next
  tt <- sub("_vs_.*$", "", cc); rf <- sub("^.*?_vs_", "", cc)
  pit <- num(rows$pi_trt[i]); pir <- num(rows$pi_ref[i]); pa <- num(rows$padj[i])
  if (is.na(pit) || is.na(pir) || is.na(pa)) next
  pivals[[tt]] <- pit; pivals[[rf]] <- pir; best[[cc]] <- pa
}
miss <- setdiff(conds, names(pivals))
if (length(miss)) stop("no pi for conditions: ", paste(miss, collapse = ", "))

perbase <- !isTRUE(a$raw)
len_D <- set_length(tested_D, coords); len_S <- set_length(S, coords)
if (perbase && (is.na(len_D) || is.na(len_S))) {
  message("  per-base pi undefined: ",
          if (is.na(len_D)) "distinct set is a junction (no length)"
          else "reference set did not resolve against the GFF",
          "; plotting the raw ratio")
  perbase <- FALSE
}
if (perbase) {
  conv <- lapply(pivals, pi_perbase, len_d = len_D, len_s = len_S)
  if (any(vapply(conv, function(v) is.na(v), logical(1)))) {
    message("  per-base pi undefined for some conditions; plotting raw")
  } else {
    pivals <- conv
    message(sprintf("  per-base pi: len_D %d bp, len_S %d bp", len_D, len_S))
  }
}

sig_contrasts <- rows$contrast[rows$contrast %in% names(best) & is_sig(rows$significant)]
nsig <- length(sig_contrasts)

## bipartition stack
read_bips <- function(kinds) {
  out <- NULL
  for (k in kinds) {
    p <- file.path(bipdir, sprintf("%s.bipartition.%s.txt", gv, k))
    if (!file.exists(p)) next
    tb <- read_tsv(p); tb$`_kind` <- k
    out <- if (is.null(out)) tb else rbind(out, tb)
  }
  out
}
bips <- read_bips(stack_kinds)
bkey <- function(s, d1, d2)
  paste(paste(sort(parts(s)), collapse = "|"), paste(sort(parts(d1)), collapse = "|"),
        paste(sort(parts(d2)), collapse = "|"), sep = "~")
focal_key <- bkey(S, D1, D2)

stack <- list()
if (!is.null(bips)) for (i in seq_len(nrow(bips))) {
  pp <- c(parts(bips$ref_ex_part[i]), parts(bips$setdiff1[i]), parts(bips$setdiff2[i]))
  ii <- idx[pp[pp %in% names(idx)]]
  if (length(ii)) stack[[length(stack) + 1]] <-
    list(lo = min(ii), hi = max(ii), b = bips[i, ])
}
matched <- NULL
if (!is.null(bips)) for (i in seq_len(nrow(bips)))
  if (bkey(bips$ref_ex_part[i], bips$setdiff1[i], bips$setdiff2[i]) == focal_key) {
    matched <- bips[i, ]; break
  }

sig_D <- list()
for (k in names(all_rows)) {
  rs <- all_rows[[k]]
  if (is.null(rs) || !nrow(rs)) next
  sel <- which(rs$significant == "TRUE")
  for (i in sel) {
    kk <- bkey(rs$ref_ex_part[i], rs$setdiff1[i], rs$setdiff2[i])
    wd <- if (grepl("diff1", rs$comparison[i])) "D1" else "D2"
    sig_D[[kk]] <- union(sig_D[[kk]], wd)
  }
}

ord <- order(vapply(stack, `[[`, numeric(1), "lo"),
             -(vapply(stack, `[[`, numeric(1), "hi") - vapply(stack, `[[`, numeric(1), "lo")))
stack <- stack[ord]
if (a$max_stack > 0 && length(stack) > a$max_stack) {
  fi <- which(vapply(stack, function(s) bkey(s$b$ref_ex_part, s$b$setdiff1, s$b$setdiff2) ==
                       focal_key, logical(1)))
  cidx <- if (length(fi)) fi[1] else length(stack) %/% 2
  half <- a$max_stack %/% 2
  lo_i <- max(1, cidx - half); hi_i <- min(length(stack), lo_i + a$max_stack - 1)
  n_drop <- length(stack) - (hi_i - lo_i + 1)
  stack <- stack[lo_i:hi_i]
  message(sprintf("note: panel B truncated to %d of %d bipartitions",
                  length(stack), length(stack) + n_drop))
}

## ------------------------------------------------------------------ draw
## Panel E's row pitch is made to match panel B's: both use a 1.0-unit row
## pitch with 0.9 units of margin, so equal inches-per-unit needs
## ratio_E / range_E == ratio_B / range_B.
route_list <- function(cell) lapply(Filter(nzchar, trimws(strsplit(cell, ",")[[1]])),
                                    function(p) strsplit(p, "-")[[1]])
n_tx <- function(cell) length(Filter(nzchar, trimws(strsplit(cell, ",")[[1]])))
e_groups <- list()
if (!is.null(matched)) {
  for (nd in route_list(matched$path1)) e_groups[[length(e_groups) + 1]] <- list("D1", nd)
  for (nd in route_list(matched$path2)) e_groups[[length(e_groups) + 1]] <- list("D2", nd)
}
e_n_rows <- length(e_groups)

adjacent_S <- character(0)
if (!is.null(matched)) {
  src <- matched$source; snk <- matched$sink
  for (k in names(G$edge_type)) {
    if (!identical(G$edge_type[[k]], "ex")) next
    uv <- strsplit(k, "\r")[[1]]
    if (uv[2] == src || uv[1] == snk)
      adjacent_S <- union(adjacent_S, intersect(hop_parts(uv[1], uv[2]), parts(S)))
  }
}

ratio_A <- 0.42; ratio_B <- 2.05; ratio_CD <- 1.15
A_MAI <- c(0.34, 0.55, 0.24, 0.25)   # panel A margins, inches

## panel C needs to know whether a junction caption will be printed before the
## layout is fixed, so resolve the junction-sourced sides here
sj_side <- list()
m_ev <- sj_meta[[as.character(a$event)]]
if (!is.null(m_ev)) for (g3 in list(c("D1", D1, "intron_distinct1"),
                                    c("D2", D2, "intron_distinct2"))) {
  if (length(parts(g3[2]))) next
  jj <- split_junctions(m_ev[[g3[3]]])
  if (length(jj)) sj_side[[g3[1]]] <- jj
}
range_B <- length(stack) + 0.9; range_E <- max(e_n_rows, 1) + 0.9
ratio_E <- ratio_B * range_E / range_B
per_unit_h <- 8.6 / (ratio_A + ratio_B + 1.3 + ratio_CD)
fig_h <- per_unit_h * (ratio_A + ratio_B + ratio_E + ratio_CD)
## Panel A's share must cover its margins AND leave a drawing region the size
## matplotlib gives it, so express it in inches and convert back to a ratio.
## Same for panel C: ratio_E was derived to match panel B's row pitch, which is
## a DRAWING height, so its margins have to be added on top or the routes are
## squeezed to nothing.
C_MAI_B <- if (length(sj_side)) 0.72 else 0.52
ratio_A <- (per_unit_h * ratio_A + A_MAI[1] + A_MAI[3]) / per_unit_h
ratio_E <- (per_unit_h * ratio_E + C_MAI_B + 0.24) / per_unit_h
fig_h <- per_unit_h * (ratio_A + ratio_B + ratio_E + ratio_CD)
sym <- if (!is.null(a$symbol)) a$symbol else gv

out <- path.expand(a$out)
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
ext <- tolower(tools::file_ext(out))
if (ext == "pdf") {
  pdf(out, width = 11.5, height = fig_h)
} else if (ext == "png") {
  png(out, width = 11.5, height = fig_h, units = "in", res = 200)
} else if (ext %in% c("eps", "ps")) {
  setEPS(); postscript(out, width = 11.5, height = fig_h)
} else {
  stop("unsupported output extension: ", ext)
}

layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 5), nrow = 4, byrow = TRUE),
       heights = c(ratio_A, ratio_B, ratio_E, ratio_CD), widths = c(3.1, 1))
par(xaxs = "i", yaxs = "i")
## Margins in inches. matplotlib puts the title and the x-label in the
## gridspec gap; here they come out of the panel, so keep them tight or the
## gene model in panel A is squeezed into a sliver of its allocation.
TOP <- 0.26        # room for one title line
mai_panel <- function(bottom, top = TOP) par(mai = c(bottom, 0.55, top, 0.25))
NX <- length(names_e)
rectp <- function(x, y, w, h, col, border = NA, lwd = 0.4)
  rect(x, y, x + w, y + h, col = col, border = border, lwd = lwd)

## ---- A: gene model -------------------------------------------------------
EH <- 1 / 6
top <- EH + 0.52
mark <- list()
if (nzchar(a$mark_parts)) for (tok in strsplit(a$mark_parts, ",")[[1]]) {
  if (!grepl(":", tok)) next
  kv <- strsplit(tok, ":", fixed = TRUE)[[1]]
  pn <- trimws(kv[1]); lb <- trimws(paste(kv[-1], collapse = ":"))
  if (pn %in% names(idx)) mark[[lb]] <- c(mark[[lb]], pn)
}
if (length(mark)) top <- max(top, EH + 0.66 + 0.30)
par(mai = A_MAI)
plot(NA, xlim = c(-0.6, NX), ylim = c(-0.34, top), axes = FALSE, xlab = "", ylab = "")
segments(0, EH / 2, NX, EH / 2, col = "#888888", lwd = 0.7)
for (e in names_e) {
  col <- CB[["other"]]
  if (e %in% parts(S)) col <- CB[["S"]]
  else if (e %in% parts(D1)) col <- CB[["D1"]]
  else if (e %in% parts(D2)) col <- CB[["D2"]]
  rectp(idx[[e]], 0, 0.86, EH, col, border = "#555555")
}
lbl <- unique(c(head(parts(S), 1), head(parts(D2), 1), head(parts(D1), 1)))
for (j in seq_along(lbl)) {
  e <- lbl[j]; if (!e %in% names(idx)) next
  dy <- if (j %% 2 == 1) EH + 0.30 else EH + 0.13
  text(idx[[e]] + 0.43, dy, e, cex = 0.5)
  segments(idx[[e]] + 0.43, EH + 0.02, idx[[e]] + 0.43, dy - 0.03, col = "#777777", lwd = 0.5)
}
if (length(mark)) {
  y0 <- EH + 0.66
  for (lb in names(mark)) {
    pp <- mark[[lb]]
    lo <- min(idx[pp]); hi <- max(idx[pp]) + 0.86
    segments(lo, y0, hi, y0, col = "#333333", lwd = 1.0, lend = "butt")
    segments(c(lo, hi), y0 - 0.05, c(lo, hi), y0, col = "#333333", lwd = 1.0)
    text((lo + hi) / 2, y0 + 0.03, lb, cex = 0.62, font = 2, col = "#222222", adj = c(0.5, 0))
  }
  mtext("brackets: named annotation regions", side = 1, line = 0.1, adj = 0.012,
        cex = 0.53, col = "#555555")
}
mtext(sprintf("A   %s exonic parts; the tested bipartition (%s strand)", sym, strand),
      side = 3, line = 0.2, adj = 0, cex = 0.78)
setlab <- function(x, cell = NULL) {
  pp <- parts(x)
  if (length(pp)) return(paste(pp, collapse = ","))
  m <- sj_meta[[as.character(a$event)]]
  if (!is.null(cell) && !is.null(m)) {
    nj <- length(split_junctions(m[[cell]]))
    if (nj) return(sprintf("%d junction%s, no exonic part", nj, if (nj > 1) "s" else ""))
  }
  "none"
}
legend(x = -0.6, y = -0.06, horiz = TRUE, bty = "n", cex = 0.58, xpd = NA, x.intersp = 0.5,
       fill = CB[c("S", "D2", "D1", "other")], border = NA,
       legend = c(sprintf("shared S (%s)", setlab(S, "intron_shared")),
                  sprintf("distinct D2 (%s)", setlab(D2, "intron_distinct2")),
                  sprintf("distinct D1 (%s)", setlab(D1, "intron_distinct1")),
                  "not in this test"))

## ---- B: every bipartition at the locus -----------------------------------
mai_panel(0.60)
plot(NA, xlim = c(-0.6, NX), ylim = c(0.2, length(stack) + 1.1), axes = FALSE,
     xlab = "", ylab = "")
n_marked <- 0
for (y in seq_along(stack)) {
  s_ <- stack[[y]]; yy <- length(stack) - y + 1
  kk <- bkey(s_$b$ref_ex_part, s_$b$setdiff1, s_$b$setdiff2)
  if (kk == focal_key) rect(-0.6, yy - 0.46, NX, yy + 0.46, col = "#F2F2F2", border = NA)
  segments(s_$lo, yy, s_$hi + 0.86, yy, col = "#BBBBBB", lwd = 0.8)
  if (!is.null(sig_D[[kk]])) {
    text(-0.38, yy, paste0("*", paste(sort(sig_D[[kk]]), collapse = "/")),
         adj = c(1, 0.5), cex = 0.6, col = "#222222", xpd = NA)
    n_marked <- n_marked + 1
  }
  for (ci in list(c("ref_ex_part", "S"), c("setdiff1", "D1"), c("setdiff2", "D2")))
    for (p in parts(s_$b[[ci[1]]]))
      if (p %in% names(idx)) rectp(idx[[p]], yy - 0.32, 0.86, 0.64, CB[[ci[2]]])
  if (kk == focal_key) {
    arrows(s_$hi + 2.6, yy, s_$hi + 1.0, yy, length = 0.05, lwd = 0.9, col = "#333333", xpd = NA)
    text(s_$hi + 2.8, yy, "tested here", adj = c(0, 0.5), cex = 0.62, xpd = NA)
  }
}
mtext(sprintf("B   All %d bipartitions tested at this locus, ordered by span containment",
              length(stack)), side = 3, line = 0.2, adj = 0, cex = 0.78)
mtext("exonic parts, genomic order  ->", side = 1, line = 1.1, cex = 0.62)
mtext(sprintf(paste0("each row = one tested bipartition; rows nest where spans are contained | ",
                     "*D1 / *D2 / *D1/D2 = significant (exontest.R call rule) in at least one ",
                     "contrast [%d of %d shown]"), n_marked, length(stack)),
      side = 1, line = 2.1, adj = 0, cex = 0.55, col = "#555555")

## ---- C (panel 3): the tested bipartition as graph routes -----------------
mai_panel(C_MAI_B)
if (is.null(matched)) {
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  text(0.02, 0.5, "no matching bipartition row found for panel C", adj = 0,
       cex = 0.65, col = "#888888")
} else {
  plot(NA, xlim = c(-0.6, NX), ylim = c(0.2, e_n_rows + 0.9), axes = FALSE,
       xlab = "", ylab = "")
  box_h <- 0.64
  for (y in seq_along(e_groups)) {
    grp <- e_groups[[y]][[1]]; nodes <- e_groups[[y]][[2]]
    yy <- e_n_rows - y + 1
    real_n <- setdiff(nodes, c("R", "L"))
    xs_ok <- vapply(real_n, node_x, numeric(1)); xs_ok <- xs_ok[!is.na(xs_ok)]
    if (!length(xs_ok)) next
    all_xs <- c(xs_ok, if (length(adjacent_S)) idx[adjacent_S] else numeric(0))
    outer_lo <- min(all_xs); outer_hi <- max(all_xs) + 0.86
    segments(outer_lo, yy, outer_hi, yy, col = "#BBBBBB", lwd = 0.8)
    for (e in adjacent_S) rectp(idx[[e]], yy - box_h / 2, 0.86, box_h, CB[["S"]])
    hollow <- grp %in% names(sj_side)
    if (length(nodes) > 1) for (i in seq_len(length(nodes) - 1)) {
      for (e in hop_parts(nodes[i], nodes[i + 1])) {
        if (e %in% parts(S)) rectp(idx[[e]], yy - box_h / 2, 0.86, box_h, CB[["S"]])
        else if (hollow) rect(idx[[e]], yy - box_h / 2, idx[[e]] + 0.86, yy + box_h / 2,
                              col = NA, border = CB[[grp]], lwd = 0.6)
        else rectp(idx[[e]], yy - box_h / 2, 0.86, box_h, CB[[grp]])
      }
    }
    if (hollow) {
      rjs <- route_junctions(nodes)
      for (jp in sj_side[[grp]]) {
        on_route <- any(vapply(rjs, function(r) identical(as.integer(r), as.integer(jp)),
                               logical(1)))
        if (!on_route) next
        x0 <- gx(jp[1]); x1 <- gx(jp[2])
        if (is.na(x0) || is.na(x1)) next
        if (x1 < x0) { tmp <- x0; x0 <- x1; x1 <- tmp }
        tt <- seq(0, base::pi, length.out = 40)
        lines(x0 + (x1 - x0) * (1 - cos(tt)) / 2, yy + 0.34 * sin(tt),
              col = CB[[grp]], lwd = 1.4)
      }
    }
    start_edge <- if (strand == "-") outer_hi else outer_lo
    end_edge   <- if (strand == "-") outer_lo else outer_hi
    if (nodes[1] == "R") {
      r_x <- node_x("R")
      segments(min(r_x, start_edge), yy, max(r_x, start_edge), yy, col = CB[[grp]],
               lwd = 0.9, lty = 3)
      if (r_x <= start_edge) text(r_x, yy, "R ->", adj = c(1.1, 0.5), cex = 0.45, col = "#555555")
      else text(r_x, yy, "<- R", adj = c(-0.1, 0.5), cex = 0.45, col = "#555555")
    }
    if (nodes[length(nodes)] == "L") {
      l_x <- node_x("L")
      segments(min(end_edge, l_x), yy, max(end_edge, l_x), yy, col = CB[[grp]],
               lwd = 0.9, lty = 3)
      if (l_x >= end_edge) text(l_x, yy, "-> L", adj = c(-0.1, 0.5), cex = 0.45, col = "#555555")
      else text(l_x, yy, "L <-", adj = c(1.1, 0.5), cex = 0.45, col = "#555555")
    }
    text(-0.9, yy, grp, adj = c(1, 0.5), cex = 0.45, col = CB[[grp]], xpd = NA)
  }
  n_d1 <- sum(vapply(e_groups, function(g) g[[1]] == "D1", logical(1)))
  n_d2 <- sum(vapply(e_groups, function(g) g[[1]] == "D2", logical(1)))
  mtext(sprintf(paste0("C   Tested bipartition %s-%s as graph routes along the exons: ",
                       "source (blue) diverges into D1 vs D2 routes"),
                matched$source, matched$sink), side = 3, line = 0.2, adj = 0, cex = 0.78)
  xlab <- sprintf(paste0("exonic parts, genomic order  ->   (D1: %d route(s), %d transcript(s) ",
                         "total; D2: %d route(s), %d transcript(s) total)"),
                  n_d1, n_tx(matched$transcripts1), n_d2, n_tx(matched$transcripts2))
  mtext(xlab, side = 1, line = 1.1, cex = 0.62)
  if (length(sj_side))
    mtext(paste0(paste(sort(names(sj_side)), collapse = " / "),
                 " has no exonic distinct set: hollow boxes are parts the route TRAVERSES, ",
                 "arcs are the counted junctions, drawn on the routes that carry them"),
          side = 1, line = 2.6, cex = 0.52, col = "#555555")
}

## ---- D: path proportion across conditions --------------------------------
par(mai = c(0.45, 0.78, TOP, 0.15))
vals <- vapply(conds, function(c_) pivals[[c_]], numeric(1))
span <- max(vals) - min(vals); pad <- max(span, 0.05)
plot(seq_along(conds) - 1, vals, type = "o", pch = 19, cex = 1.1, lwd = 2,
     col = CB[["D2"]], axes = FALSE, xlab = "", ylab = "",
     xlim = c(-0.4, length(conds) - 0.6),
     ylim = c(min(vals) - pad * 0.30, max(vals) + pad * 0.30))
points(seq_along(conds) - 1, vals, pch = 21, bg = CB[["D2"]], col = "white", cex = 1.1, lwd = 1.2)
text(seq_along(conds) - 1, vals + pad * 0.11, sprintf("%.3f", vals), cex = 0.65)
axis(1, at = seq_along(conds) - 1, labels = sub("\\.CLASSIC", "", conds),
     cex.axis = 0.72, tick = FALSE, line = -0.4)
axis(2, las = 1, cex.axis = 0.7)
grid(nx = NA, ny = NULL, col = "#CCCCCC", lty = 1)
yl <- if (!is.null(a$ylab)) a$ylab else
  sprintf("%s path proportion pi\nof distinct %s",
          if (perbase) "per-base" else "raw", which_D)
mtext(yl, side = 2, line = 2.6, cex = 0.68)
mtext("D   Path proportion across conditions", side = 3, line = 0.2, adj = 0, cex = 0.78)

## ---- D text block --------------------------------------------------------
par(mai = c(0.45, 0.10, TOP, 0.10))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
text(0, 0.95, paste0(sym, if (nzchar(a$note)) paste0("\n(", a$note, ")") else ""),
     adj = c(0, 1), cex = 0.72)
desc <- if (!is.null(a$desc)) a$desc else
  sprintf("distinct %s %s\nagainst shared S %s", which_D, tested_D, S)
text(0, 0.64, desc, adj = c(0, 1), cex = 0.65)
if (nsig) {
  shown <- NULL
  if (!is.null(a$feature_contrast) && a$feature_contrast %in% sig_contrasts) {
    shown <- a$feature_contrast
  } else {
    if (!is.null(a$feature_contrast))
      message("note: --feature_contrast ", a$feature_contrast,
              " is not among the called contrasts; using the largest effect")
    dsub <- rows[rows$contrast %in% sig_contrasts, , drop = FALSE]
    shown <- dsub$contrast[which.max(abs(num(dsub$delta_pi)))]
  }
  t_c <- sub("_vs_.*$", "", shown); r_c <- sub("^.*?_vs_", "", shown)
  head_txt <- sprintf("pi %.3f (%s) -> %.3f (%s)   [%s]",
                      pivals[[r_c]], r_c, pivals[[t_c]], t_c, gsub("_vs_", " vs ", shown))
  if (nsig > 1) head_txt <- paste0(head_txt, sprintf("\nand %d other called contrast(s)", nsig - 1))
} else {
  head_txt <- sprintf("pi %.3f (%s) -> %.3f (%s)   [no contrast called; range shown]",
                      pivals[[conds[1]]], conds[1], pivals[[conds[length(conds)]]], conds[length(conds)])
}
text(0, 0.34, sprintf("%s\nsignificant (exontest.R call) in %d of %d contrasts",
                      head_txt, nsig, length(best)), adj = c(0, 1), cex = 0.65)

invisible(dev.off())
message(sprintf("wrote %s   (%d bipartitions, comparison %s, %d sig contrasts)",
                a$out, length(stack), comp, nsig))
