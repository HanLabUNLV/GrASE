#!/usr/bin/env Rscript
## GrASE, rMATS and DEXSeq compared on exonic parts, the one unit all three reach through the
## splicing graph, with an UpSet plot of which tools call each part.
##
## Usage:
##   Rscript compare_parts_upset.R \
##     --grase  <annotated results, comma separated, e.g. internal and terminal> \
##     --map_dirs <rMATS mapping directories, comma separated> \
##     --rmats_dir <dir with {SE,A3SS,A5SS,RI}.MATS.JCEC.txt> \
##     --dexseq <DEXSeq results table> --outdir <dir> [--contrast NAME] [options]
##
## One adjusted-p cutoff (--pval) for all three tools, and an effect size for each. A part is
##   DEXSeq  called  padj < pval and |log2 fold change of exon usage| >= --min_l2fc
##   GrASE   called  in the distinct sets (setdiff1 and setdiff2) of a bipartition with either side
##                   called by grase::is_significant(), the rule behind the `significant` column, with
##                   pval as its p-value threshold (--min_dpi, --min_dpi_sj, --min_reads). A
##                   junction-sourced side has no exonic part.
##   rMATS   called  in the distinct set of an event matched to a bipartition, with
##                   FDR < pval and |IncLevelDifference| >= --min_dpsi
## A part is tested when the tool gave it a p-value: DEXSeq a padj, GrASE a bipartition with a p-value
## whose distinct sets hold it, rMATS a matched event whose distinct set holds it.
##
## --grase takes one GrASE results file per class (internal, terminal), named by --class. --map_dirs
## gives the rMATS mappings; every one given is used, so the graph mapping may be given once per class
## and the bipartition mapping as map.internal,map.terminal. Classes are a property of GrASE's
## bipartitions only, so rMATS and DEXSeq are never split by class; the comparisons that use the
## classes put both GrASE classes in the same plot against the same rMATS and DEXSeq calls.
##
## Writes exonic_parts_three_tools.tsv: one row per part in the DEXSeq table or tested by GrASE or
## rMATS, with each tool's tested and called flags (GrASE per class and combined), the UpSet category
## of each comparison, how many of the three tools tested the part, whether all three did, and the
## evidence behind each call: DEXSeq padj and fold change; per GrASE class, the events whose distinct
## sets hold the part, those called, and their smallest padj; the rMATS events likewise, with their
## smallest FDR. GrASE and rMATS parts that find no DEXSeq row (e.g. removed by a count filter before
## DEXSeq or satuRn) are kept with dex_tested FALSE; their number is printed.
##
## Two comparisons, two UpSet plots each:
##   upset_exonic_parts             GrASE (either class), rMATS, DEXSeq
##   upset_exonic_parts.by_class    GrASE internal, GrASE terminal, rMATS, DEXSeq (with two classes)
## <stem>.{pdf,png} covers every part called by at least one set, each bar split by how many of the
## three tools tested the part (3, 2 or 1): in a bar, parts tested by all three were examined by the
## tools outside it and not called, parts tested only by the tools in it could not be examined by the
## others. Its counts are in <stem>.tsv. <stem>.tested_by_all.{pdf,png} restricts to parts all three
## tools tested; its bars are the "3" segments of the first plot, on their own scale. The rows behind
## each bar are in categories/{all,by_class}.<category>.tsv.
## Axes: unless --ymax, --setmax and their _tested_by_all forms are given, each axis runs to the largest
## value over every plot drawn in the run, plus 5%, so all of them share one scale.
## --count writes the tables only. --draw <per-part tables, comma separated> draws from tables
## written earlier, all on one scale, each into its table's directory, without reading any tool's
## results: count each contrast with --count, then draw them together.
suppressPackageStartupMessages({library(optparse); library(data.table); library(UpSetR); library(ggplot2); library(patchwork); library(grase)})
o <- parse_args(OptionParser(option_list = list(
  make_option("--grase",     type = "character", help = "annotated exontest.R results, comma separated"),
  make_option("--map_dirs",  type = "character", help = "rMATS mapping directories (map_rmats_graph.R or map_rmats_bipartitions.R), comma separated"),
  make_option("--rmats_dir", type = "character"),
  make_option("--dexseq",    type = "character", help = "DEXSeq results, or a satuRn table on DEXSeq exonic parts"),
  make_option("--outdir",    type = "character"),
  make_option("--contrast",  type = "character", default = NULL, help = "keep only this contrast of the GrASE results"),
  make_option("--lfc_col",   type = "character", default = NULL, help = "DEXSeq fold change column [first log2fold_*]"),
  make_option("--etypes",    type = "character", default = "SE,A3SS,A5SS,RI", help = "rMATS event types [%default]"),
  make_option("--pval",      type = "double", default = 0.05, help = "adjusted p cutoff for all three tools [%default]"),
  make_option("--min_dpsi",  type = "double", default = 0.1,  help = "rMATS |IncLevelDifference| [%default]"),
  make_option("--min_l2fc",  type = "double", default = 0.5,  help = "DEXSeq |log2 fold change| [%default]"),
  make_option("--min_dpi",   type = "double", default = 0.1,  help = "GrASE |delta_pi| [%default]"),
  make_option("--min_dpi_sj", type = "double", default = 0.05, help = "GrASE |delta_pi| for a junction-sourced side [%default]"),
  make_option("--class",    type = "character", default = NULL, help = "GrASE class of each --grase file, comma separated [internal,terminal]"),
  make_option("--font",     type = "character", default = "Arial", help = "font of every plot, embedded in the PDFs [%default]"),
  make_option("--ymax",     type = "double", default = NULL, help = "bar axis limit of the stacked plots, to share a scale across runs [per plot]"),
  make_option("--setmax",   type = "double", default = NULL, help = "parts-called axis limit of the stacked plots [per plot]"),
  make_option("--ymax_tested_by_all",   type = "double", default = NULL, help = "bar axis limit of the tested_by_all plots [per plot]"),
  make_option("--setmax_tested_by_all", type = "double", default = NULL, help = "set size axis limit of the tested_by_all plots [per plot]"),
  make_option("--count",    action = "store_true", default = FALSE, help = "count and write the tables only, drawing nothing"),
  make_option("--draw",     type = "character", default = NULL, help = "per-part tables written earlier, comma separated: draw them on one scale, each into its directory, reading no tool results"),
  make_option("--min_reads", type = "double", default = 10,   help = "GrASE reads on the distinct set in a group [%default]"))))
if (is.null(o$draw)) for (r in c("grase", "map_dirs", "rmats_dir", "dexseq", "outdir"))
  if (is.null(o[[r]])) stop("--", r, " is required", call. = FALSE)
split <- function(x) trimws(strsplit(x, ",")[[1]])
if (is.null(o$draw)) dir.create(o$outdir, showWarnings = FALSE, recursive = TRUE)
## both devices fall back silently to the system sans font when the requested one is not installed
font_used <- systemfonts::match_fonts(o$font)$path
if (!grepl(o$font, basename(font_used), ignore.case = TRUE))
  warning(sprintf("font %s is not installed; plots use %s", o$font, basename(font_used)), call. = FALSE)
if (is.null(o$draw)) cat(sprintf("pval %g | rMATS |dPSI| >= %g | DEXSeq |log2FC| >= %g | GrASE reads >= %g\n", o$pval, o$min_dpsi, o$min_l2fc, o$min_reads))

cat_of <- function(D) apply(as.matrix(D), 1, function(r) { x <- names(D)[r]; if (length(x)) paste(x, collapse = "+") else "none" })
## the two comparisons: GrASE as one tool, and GrASE split by class, each against the same rMATS and
## DEXSeq; the classes are read off the per-part table's grase_<class>_called columns
make_cmp <- function(P) {
  cls <- sub("^grase_(.*)_called$", "\\1", grep("^grase_.+_called$", names(P), value = TRUE))
  cmp <- list(all = setNames(P[, .(grase_called, rmats_called, dex_called)], c("GrASE", "rMATS", "DEXSeq")))
  if (length(cls) > 1)
    cmp$by_class <- setNames(P[, c(paste0("grase_", cls, "_called"), "rmats_called", "dex_called"), with = FALSE],
                             c(paste0("GrASE_", cls), "rMATS", "DEXSeq"))
  cmp
}

if (is.null(o$draw)) {
  ## DEXSeq: exonic parts with a DEXSeq row. A table without a pvalue column is read as satuRn
  ## (ids = gene:part, empirical_FDR, estimates), as in compare_dexseq_grase.R
  hd <- names(fread(o$dexseq, nrows = 0))
  if ("pvalue" %in% hd) {
    lfc <- if (!is.null(o$lfc_col)) o$lfc_col else grep("^log2fold_", hd, value = TRUE)[1]
    if (is.na(lfc) || !lfc %in% hd) stop("no DEXSeq fold change column (--lfc_col); columns: ", paste(hd, collapse = " "))
    d <- fread(o$dexseq, select = c("groupID", "featureID", "padj", lfc), colClasses = list(character = c("groupID", "featureID")))
    setnames(d, c("groupID", "featureID", lfc), c("gene", "part", "l2fc"))
  } else {
    if (!all(c("ids", "empirical_FDR", "estimates") %in% hd)) stop("--dexseq is neither DEXSeq nor satuRn; columns: ", paste(hd, collapse = " "))
    cat("--dexseq read as satuRn (empirical FDR, estimates)\n")
    d <- fread(o$dexseq, select = c("ids", "empirical_FDR", "estimates"))
    d[, c("gene", "part") := tstrsplit(ids, ":", fixed = TRUE)]
    d <- d[, .(gene, part, padj = empirical_FDR, l2fc = estimates)]
  }
  d[, gene := gsub('"', "", gene)][, part := gsub('"', "", part)]
  cat(sprintf("DEXSeq groups merging genes: %d\n", uniqueN(grep("+", d$gene, fixed = TRUE, value = TRUE))))
  d[, `:=`(dex_tested = !is.na(padj), dex_called = !is.na(padj) & padj < o$pval & !is.na(l2fc) & abs(l2fc) >= o$min_l2fc)]
  P <- d[, .(gene, part, dex_tested, dex_called)]
  min_na <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  ids <- function(x) if (length(x)) paste(sort(unique(x)), collapse = ",") else NA_character_


  gfiles <- split(o$grase); mdirs <- split(o$map_dirs)
  if (length(gfiles) != length(mdirs)) stop("--grase and --map_dirs must pair one to one")
  labels <- if (is.null(o$class)) c("internal", "terminal")[seq_along(gfiles)] else split(o$class)
  et <- split(o$etypes)
  P <- copy(d[, .(gene, part, dex_tested, dex_called, dex_padj = padj, dex_l2fc = l2fc)])

  ## GrASE, per class: parts of the distinct sets of both sides of a bipartition
  for (k in seq_along(gfiles)) {
    cl <- labels[k]
    g <- fread(gfiles[k], colClasses = list(character = c("gene", "setdiff1", "setdiff2")))
    if (!is.null(o$contrast)) g <- g[contrast == o$contrast]
    if (!nrow(g)) stop("no rows for --contrast ", o$contrast, "; the file has: ", paste(unique(fread(gfiles[k], select = "contrast")$contrast), collapse = " "))
    g[, call := is_significant(g, o$pval, 0, o$min_dpi, o$min_dpi_sj, o$min_reads)]
    ## tested = has a p-value; padj is NA for tests in genes not selected at the gene level
    g <- g[!is.na(p.value)]
    ## each row (one side of a bipartition) spreads over both distinct sets of its bipartition
    gp <- rbind(g[, .(gene, event, padj, call, part = setdiff1)], g[, .(gene, event, padj, call, part = setdiff2)])
    gp <- gp[!is.na(part) & part != "NA" & nzchar(part)]
    gp <- gp[, .(part = trimws(unlist(strsplit(part, ",")))), by = .(gene, event, padj, call)]
    gp <- gp[, .(called = any(call), events = ids(event), called_events = ids(event[call]), padj = min_na(padj)), by = .(gene, part)]
    cat(sprintf("GrASE %s: %d bipartition-contrasts tested, %d called; %d of %d parts have no DEXSeq row\n", cl,
                nrow(unique(g[, .(gene, event, contrast)])), nrow(unique(g[call == TRUE, .(gene, event, contrast)])),
                nrow(fsetdiff(gp[, .(gene, part)], d[, .(gene, part)])), nrow(gp)))
    setnames(gp, c("called", "events", "called_events", "padj"), paste0("grase_", cl, c("_called", "_events", "_called_events", "_padj")))
    P <- merge(P, gp, by = c("gene", "part"), all = TRUE)
  }

  ## rMATS: parts of the events' distinct sets, over every mapping given (the same map may be given per class)
  rt <- unique(rbindlist(lapply(unique(mdirs), function(md) rbindlist(lapply(et, function(t) {
    m <- fread(file.path(o$rmats_dir, paste0(t, ".MATS.JCEC.txt")), select = c("ID", "FDR", "IncLevelDifference"))
    m[, tag := paste0(t, "_", ID)]
    x <- fread(file.path(md, sprintf("combined.dexseq.%s.mapped.txt", t)), colClasses = "character")
    setnames(x, 1:3, c("GeneID", "DexseqFragment", "ids"))
    mp <- x[!is.na(ids) & ids != "NA", .(tag = unlist(strsplit(ids, ","))), by = .(gene = gsub('"', "", GeneID), part = DexseqFragment)]
    merge(mp, m[, .(tag, FDR, IncLevelDifference)], by = "tag")
  })))))
  rt <- rt[!is.na(FDR)]
  rt[, call := FDR < o$pval & abs(IncLevelDifference) >= o$min_dpsi]
  rp <- rt[, .(rmats_called = any(call), rmats_events = ids(tag), rmats_called_events = ids(tag[call]),
               rmats_fdr = min_na(FDR)), by = .(gene, part)]
  cat(sprintf("rMATS: %d parts; %d have no DEXSeq row\n", nrow(rp), nrow(fsetdiff(rp[, .(gene, part)], d[, .(gene, part)]))))
  P <- merge(P, rp, by = c("gene", "part"), all = TRUE)

  ## rows are the union of the parts each tool tested: a part missing from a tool's table is untested by it.
  ## GrASE as one tool tests or calls a part when either class does.
  gcl <- paste0("grase_", labels, "_called"); gte <- paste0("grase_", labels, "_tested")
  P[, (gte) := lapply(.SD, function(x) !is.na(x)), .SDcols = gcl]
  P[, rmats_tested := !is.na(rmats_called)]
  for (v in c("dex_tested", "dex_called", "rmats_called", gcl)) set(P, which(is.na(P[[v]])), v, FALSE)
  P[, grase_tested := Reduce(`|`, .SD), .SDcols = gte][, grase_called := Reduce(`|`, .SD), .SDcols = gcl]
  P[, tested_by_all := dex_tested & grase_tested & rmats_tested]
  ## tools that tested the part, of the three
  P[, n_tested := grase_tested + rmats_tested + dex_tested]
  stopifnot(P[, all((!grase_called | grase_tested) & (!rmats_called | rmats_tested) & (!dex_called | dex_tested))])

  cmp <- make_cmp(P)
  P[, category := cat_of(cmp$all)]
  if (!is.null(cmp$by_class)) P[, category_by_class := cat_of(cmp$by_class)]

  cat(sprintf("exonic parts (union) %d | tested: DEXSeq %d, GrASE %d (%s), rMATS %d, all three %d\n", nrow(P),
              sum(P$dex_tested), sum(P$grase_tested), paste(sprintf("%s %d", labels, colSums(P[, ..gte])), collapse = ", "),
              sum(P$rmats_tested), sum(P$tested_by_all)))
  cat(sprintf("called: DEXSeq %d, GrASE %d (%s), rMATS %d\n", sum(P$dex_called), sum(P$grase_called),
              paste(sprintf("%s %d", labels, colSums(P[, ..gcl])), collapse = ", "), sum(P$rmats_called)))
  setcolorder(P, intersect(c("gene", "part", "category", "category_by_class", "tested_by_all", "n_tested"), names(P)))
  fwrite(P, file.path(o$outdir, "exonic_parts_three_tools.tsv"), sep = "\t", na = "NA", quote = FALSE)
  dir.create(file.path(o$outdir, "categories"), showWarnings = FALSE)
  ## a category empty in this run must not keep a file from an earlier one
  unlink(Sys.glob(file.path(o$outdir, "categories", "*.tsv")))
  for (nm in names(cmp)) {
    col <- if (nm == "all") "category" else "category_by_class"
    for (cg in setdiff(unique(P[[col]]), "none"))
      fwrite(P[get(col) == cg], file.path(o$outdir, "categories", sprintf("%s.%s.tsv", nm, cg)), sep = "\t", na = "NA", quote = FALSE)
  }
  if (o$count) quit(save = "no")
  Ps <- list(P); outdirs <- o$outdir
} else {
  tables <- split(o$draw)
  Ps <- lapply(tables, fread); outdirs <- dirname(tables)
}

th <- theme_minimal(base_size = 11, base_family = o$font) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
ylab <- "exonic parts"

## UpSetR plot of the parts all three tools tested
up_tba <- function(D, title, file) {
  B <- as.data.frame(lapply(D, as.integer)); B <- B[rowSums(B) > 0, ]
  cat(title, ": ", nrow(B), " parts called by at least one tool\n", sep = "")
  print(as.data.table(B)[, .N, by = names(B)][order(-N)])
  draw <- function() print(upset(B, sets = names(B), keep.order = TRUE, order.by = "freq", nintersects = NA,
                                 mainbar.y.label = ylab, sets.x.label = "parts called", text.scale = 1.3,
                                 mainbar.y.max = o$ymax_tested_by_all, set_size.scale_max = o$setmax_tested_by_all))
  cairo_pdf(file, width = 8, height = 5, family = o$font); draw(); dev.off()
  png(sub("\\.pdf$", ".png", file), width = 2400, height = 1500, res = 300, type = "cairo", family = o$font); draw(); dev.off()
}

## every part called by any tool, each bar split by how many of the three tools tested the part
up_split <- function(D, n_tested, title, file) {
  cg <- cat_of(D); keep <- cg != "none"
  lev <- c("3", "2", "1")
  S <- data.table(category = cg[keep], others = as.character(n_tested[keep]))[, .N, by = .(category, others)]
  tot <- S[, .(N = sum(N)), by = category][order(-N)]
  fwrite(setnames(dcast(S, category ~ factor(others, levels = lev), value.var = "N", fill = 0L)[match(tot$category, category)],
                  lev, paste0("tested_by_", lev), skip_absent = TRUE),
         sub("\\.pdf$", ".tsv", file), sep = "\t", quote = FALSE)
  cat(title, ": ", sum(tot$N), " parts called by at least one tool\n", sep = "")
  print(S[order(match(category, tot$category), match(others, lev))])
  S[, category := factor(category, levels = tot$category)][, others := factor(others, levels = lev)]
  ## UpSet layout: stacked intersection bars on top, the set matrix below, each set's calls at its left
  sets <- names(D)                                           # matrix rows bottom to top
  M <- CJ(category = factor(tot$category, levels = tot$category), set = factor(sets, levels = sets))
  M[, on := mapply(function(c, s) s %in% strsplit(c, "+", fixed = TRUE)[[1]], as.character(category), as.character(set))]
  L <- M[(on), .(lo = min(as.integer(set)), hi = max(as.integer(set))), by = category]
  ss <- data.table(set = factor(sets, levels = sets), N = colSums(D))
  ## columns fixed to the bar order: a layer missing a category would otherwise move it to the end
  xs <- scale_x_discrete(limits = tot$category)
  top <- ggplot(S, aes(category, N, fill = others)) +
    geom_col(width = 0.55, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
    geom_text(data = tot[, .(category = factor(category, levels = tot$category), N)], aes(category, N, label = N),
              inherit.aes = FALSE, vjust = -0.4, size = 3, colour = "grey20", family = o$font) + xs +
    scale_fill_manual(values = setNames(c("#104281", "#2a78d6", "#86b6ef"), lev), drop = FALSE,
                      name = "Tools that tested the part") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)), limits = if (is.null(o$ymax)) NULL else c(0, o$ymax)) +
    guides(fill = guide_legend(reverse = TRUE, nrow = 1)) +
    labs(x = NULL, y = ylab) + th +
    theme(axis.text.x = element_blank(), legend.position = "top", legend.justification = "left")
  mat <- ggplot(M, aes(category, set)) + xs +
    geom_point(data = M[(!on)], size = 3, colour = "grey85") +
    geom_segment(data = L, aes(x = category, xend = category, y = lo, yend = hi), inherit.aes = FALSE,
                 linewidth = 0.8, colour = "grey20") +
    geom_point(data = M[(on)], size = 3, colour = "grey20") +
    labs(x = NULL, y = NULL) + th +
    theme(axis.text.x = element_blank(), panel.grid.major.y = element_blank(), axis.text.y = element_text(hjust = 0.5))
  left <- ggplot(ss, aes(N, set)) + geom_col(width = 0.6, fill = "grey35") +
    scale_x_reverse(expand = expansion(mult = c(0.02, 0)), labels = function(x) ifelse(x >= 1000, paste0(x / 1000, "K"), x),
                    limits = if (is.null(o$setmax)) NULL else c(o$setmax, 0)) +
    labs(x = "parts called", y = NULL) + th +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
  g <- patchwork::wrap_plots(patchwork::plot_spacer(), top, left, mat, ncol = 2, widths = c(1, 4),
                             heights = c(3, 0.35 * length(sets)))
  w <- max(8, 2.5 + 0.45 * nrow(tot))
  ggsave(file, g, width = w, height = 5, device = cairo_pdf)
  ggsave(sub("\\.pdf$", ".png", file), g, width = w, height = 5, dpi = 300)
}

## one scale for every plot drawn: the largest bar and set of each plot type, plus 5%, unless given
ax <- c("ymax", "setmax", "ymax_tested_by_all", "setmax_tested_by_all")
if (any(vapply(ax, function(a) is.null(o[[a]]), logical(1)))) {
  lim <- function(D) { cg <- cat_of(D); cg <- cg[cg != "none"]; c(if (length(cg)) max(table(cg)) else 0, max(colSums(D))) }
  L <- do.call(rbind, lapply(Ps, function(P) do.call(rbind, lapply(make_cmp(P), function(D) c(lim(D), lim(D[P$tested_by_all]))))))
  auto <- ceiling(apply(L, 2, max) * 1.05)
  for (k in seq_along(ax)) if (is.null(o[[ax[k]]])) o[[ax[k]]] <- auto[k]
}
cat(sprintf("axes: bars %g, sets %g; tested_by_all bars %g, sets %g\n", o$ymax, o$setmax, o$ymax_tested_by_all, o$setmax_tested_by_all))

for (i in seq_along(Ps)) {
  P <- Ps[[i]]; cmp <- make_cmp(P)
  for (nm in names(cmp)) {
    stem <- if (nm == "all") "upset_exonic_parts" else "upset_exonic_parts.by_class"
    cat("\n=====", outdirs[i], nm, "=====\n")
    up_split(cmp[[nm]], P$n_tested, paste(nm, "called, any tool"), file.path(outdirs[i], paste0(stem, ".pdf")))
    up_tba(cmp[[nm]][P$tested_by_all], paste(nm, "called, parts tested by all three"),
           file.path(outdirs[i], paste0(stem, ".tested_by_all.pdf")))
  }
}
