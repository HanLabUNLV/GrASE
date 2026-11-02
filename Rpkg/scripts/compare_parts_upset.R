#!/usr/bin/env Rscript
## GrASE, rMATS and DEXSeq compared on exonic parts, the one unit all three reach through the
## splicing graph, with an UpSet plot of which tools call each part.
##
## Usage:
##   Rscript compare_parts_upset.R \
##     --grase  <annotated results, comma separated, e.g. internal and terminal> \
##     --map_dirs <map_rmats_bipartitions.R --outdir of each class, comma separated> \
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
## --grase and --map_dirs are paired by position, one pair per class (internal, terminal), and
## each class is compared on its own: rMATS has no TSS/TTS events, so pooling the terminal class
## would make the GrASE-only bar large for a structural reason. Per class, it writes
## exonic_parts_three_tools.<class>.tsv (one row per part in the DEXSeq table or tested by GrASE or rMATS, with each tool's tested/called
## flags, its UpSet category, whether all three tested it, and the evidence behind each call:
## DEXSeq padj and fold change; the GrASE events whose distinct sets hold the part, those called,
## and their smallest padj; the rMATS events likewise, with their smallest FDR) and
## upset_exonic_parts.<class>{,.tested_by_all}.{pdf,png}. The second restricts to parts
## all three tested, so a part one tool never examined is not counted against it. The rows of each
## category called by at least one tool are also written to categories/<class>.<category>.tsv.
## GrASE and rMATS parts that find no DEXSeq row (e.g. removed by a count filter before DEXSeq or satuRn)
## are kept with dex_tested FALSE; their number is printed.
suppressPackageStartupMessages({library(optparse); library(data.table); library(UpSetR); library(grase)})
o <- parse_args(OptionParser(option_list = list(
  make_option("--grase",     type = "character", help = "annotated exontest.R results, comma separated"),
  make_option("--map_dirs",  type = "character", help = "map_rmats_bipartitions.R output directories, comma separated"),
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
  make_option("--class",    type = "character", default = NULL, help = "class names, comma separated [internal,terminal]"),
  make_option("--min_reads", type = "double", default = 10,   help = "GrASE reads on the distinct set in a group [%default]"))))
for (r in c("grase", "map_dirs", "rmats_dir", "dexseq", "outdir"))
  if (is.null(o[[r]])) stop("--", r, " is required", call. = FALSE)
split <- function(x) trimws(strsplit(x, ",")[[1]])
dir.create(o$outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("pval %g | rMATS |dPSI| >= %g | DEXSeq |log2FC| >= %g | GrASE reads >= %g\n", o$pval, o$min_dpsi, o$min_l2fc, o$min_reads))

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
sets <- c("GrASE", "rMATS", "DEXSeq")
et <- split(o$etypes)

for (k in seq_along(gfiles)) {
  cl <- labels[k]; cat("\n=====", cl, "=====\n")
  P <- copy(d[, .(gene, part, dex_tested, dex_called, dex_padj = padj, dex_l2fc = l2fc)])

  ## GrASE: parts of the distinct sets of both sides of a bipartition
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
  gp <- gp[, .(grase_called = any(call), grase_events = ids(event), grase_called_events = ids(event[call]),
               grase_padj = min_na(padj)), by = .(gene, part)]
  cat(sprintf("GrASE: %d bipartition-contrasts tested, %d called; %d of %d parts have no DEXSeq row\n",
              nrow(unique(g[, .(gene, event, contrast)])), nrow(unique(g[call == TRUE, .(gene, event, contrast)])),
              nrow(fsetdiff(gp[, .(gene, part)], d[, .(gene, part)])), nrow(gp)))
  P <- merge(P, gp, by = c("gene", "part"), all = TRUE)

  ## rMATS: parts of the matched events' distinct sets
  rt <- rbindlist(lapply(et, function(t) {
    m <- fread(file.path(o$rmats_dir, paste0(t, ".MATS.JCEC.txt")), select = c("ID", "FDR", "IncLevelDifference"))
    m[, tag := paste0(t, "_", ID)]
    x <- fread(file.path(mdirs[k], sprintf("combined.dexseq.%s.mapped.txt", t)), colClasses = "character")
    setnames(x, 1:3, c("GeneID", "DexseqFragment", "ids"))
    mp <- x[!is.na(ids) & ids != "NA", .(tag = unlist(strsplit(ids, ","))), by = .(gene = gsub('"', "", GeneID), part = DexseqFragment)]
    merge(mp, m[, .(tag, FDR, IncLevelDifference)], by = "tag")
  }))
  rt <- rt[!is.na(FDR)]
  rt[, call := FDR < o$pval & abs(IncLevelDifference) >= o$min_dpsi]
  rp <- rt[, .(rmats_called = any(call), rmats_events = ids(tag), rmats_called_events = ids(tag[call]),
               rmats_fdr = min_na(FDR)), by = .(gene, part)]
  cat(sprintf("rMATS: %d parts; %d have no DEXSeq row\n", nrow(rp), nrow(fsetdiff(rp[, .(gene, part)], d[, .(gene, part)]))))
  P <- merge(P, rp, by = c("gene", "part"), all = TRUE)
  ## rows are the union of the parts each tool tested: a part missing from a tool's table is untested by it
  P[, `:=`(grase_tested = !is.na(grase_called), rmats_tested = !is.na(rmats_called))]
  for (v in c("dex_tested", "dex_called", "grase_called", "rmats_called")) set(P, which(is.na(P[[v]])), v, FALSE)

  cat(sprintf("exonic parts (union) %d | tested: DEXSeq %d, GrASE %d, rMATS %d, all three %d\n", nrow(P),
              sum(P$dex_tested), sum(P$grase_tested), sum(P$rmats_tested), sum(P$dex_tested & P$grase_tested & P$rmats_tested)))
  cat(sprintf("called: DEXSeq %d, GrASE %d, rMATS %d\n", sum(P$dex_called), sum(P$grase_called), sum(P$rmats_called)))
  P[, category := apply(cbind(ifelse(grase_called, "GrASE", NA), ifelse(rmats_called, "rMATS", NA),
                              ifelse(dex_called, "DEXSeq", NA)), 1, function(x) {
    x <- x[!is.na(x)]; if (length(x)) paste(x, collapse = "+") else "none" })]
  P[, tested_by_all := dex_tested & grase_tested & rmats_tested]
  setcolorder(P, c("gene", "part", "category", "tested_by_all"))
  fwrite(P, file.path(o$outdir, sprintf("exonic_parts_three_tools.%s.tsv", cl)), sep = "\t", na = "NA", quote = FALSE)
  dir.create(file.path(o$outdir, "categories"), showWarnings = FALSE)
  for (cg in setdiff(unique(P$category), "none"))
    fwrite(P[category == cg], file.path(o$outdir, "categories", sprintf("%s.%s.tsv", cl, cg)), sep = "\t", na = "NA", quote = FALSE)

  up <- function(D, title, file) {
    B <- data.frame(GrASE = as.integer(D$grase_called), rMATS = as.integer(D$rmats_called), DEXSeq = as.integer(D$dex_called))
    B <- B[rowSums(B) > 0, ]
    cat(title, ": ", nrow(B), " parts called by at least one tool\n", sep = "")
    print(as.data.table(B)[, .N, by = .(GrASE, rMATS, DEXSeq)][order(-N)])
    draw <- function() print(upset(B, sets = sets, keep.order = TRUE, order.by = "freq",
                                   mainbar.y.label = paste("exonic parts,", cl), sets.x.label = "parts called", text.scale = 1.3))
    pdf(file, width = 7, height = 5); draw(); dev.off()
    png(sub("\\.pdf$", ".png", file), width = 2100, height = 1500, res = 300); draw(); dev.off()
  }
  up(P, paste(cl, "called, any tool"), file.path(o$outdir, sprintf("upset_exonic_parts.%s.pdf", cl)))
  up(P[dex_tested & grase_tested & rmats_tested], paste(cl, "called, parts tested by all three"),
     file.path(o$outdir, sprintf("upset_exonic_parts.%s.tested_by_all.pdf", cl)))
}
