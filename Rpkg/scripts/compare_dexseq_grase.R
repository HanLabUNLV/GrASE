#!/usr/bin/env Rscript
## GrASE against an exon-level tool (DEXSeq, or satuRn run on DEXSeq exonic parts) on exonic parts.
##
## Each GrASE test is spread over the parts of its distinct set (setdiff1 for diff1_vs_ref,
## setdiff2 for diff2_vs_ref); a part in several tests keeps the one with the smallest padj.
## The exon-level table is read as DEXSeq (groupID, featureID, pvalue, padj, log2fold_*) or,
## when it has no pvalue column, as satuRn (ids = gene:part, empirical_pval, empirical_FDR,
## estimates), the same detection compare_rmats_dexseq.R uses.
##
## Calls use adjusted p-values on both sides: GrASE padj and DEXSeq padj / satuRn empirical FDR,
## each below --pval. The p-value plots use the unadjusted values.
##
## Writes <output> (one row per part both tools tested), <output>.discrepancies.txt (parts called
## by one tool only) and <output>.pdf.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-g", "--grase"), type = "character", default = NULL,
              help = "GrASE annotated exontest results, comma separated (e.g. internal and terminal)"),
  make_option(c("-e", "--exon_results"), type = "character", default = NULL,
              help = "DEXSeq results table, or a satuRn table on DEXSeq exonic parts"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "output table; the PDF and discrepancy list are named after it"),
  make_option("--contrast", type = "character", default = NULL,
              help = "GrASE contrast to compare; required when the GrASE results hold more than one"),
  make_option("--pval", type = "double", default = 0.05,
              help = "cutoff on GrASE padj and on DEXSeq padj / satuRn empirical FDR [default= %default]")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
for (nm in c("grase", "exon_results", "output")) {
  if (is.null(opt[[nm]])) {
    print_help(opt_parser)
    stop("--", nm, " is required and has no default", call. = FALSE)
  }
}
grase_files <- trimws(strsplit(opt$grase, ",")[[1]])
for (f in c(grase_files, opt$exon_results)) if (!file.exists(f)) stop("not found: ", f)
output_file <- opt$output
plot_file <- sub("\\.[^.]+$", ".pdf", output_file)
if (plot_file == output_file) plot_file <- paste0(output_file, ".pdf")
discrep_file <- sub("\\.[^.]+$", ".discrepancies.txt", output_file)
if (discrep_file == output_file) discrep_file <- paste0(output_file, ".discrepancies.txt")

## --- Exon-level tool ---
read_exon_results <- function(path) {
  x <- fread(path)
  if ("pvalue" %in% names(x)) {
    lfc <- grep("^log2fold_", names(x), value = TRUE)[1]
    if (is.na(lfc)) stop("DEXSeq table has no log2fold_* column")
    message("Exon-level table read as DEXSeq (effect: ", lfc, ")")
    out <- x[, .(gene = groupID, part = featureID, exon_p = pvalue, exon_padj = padj, exon_effect = get(lfc))]
  } else {
    need <- c("ids", "empirical_pval", "empirical_FDR", "estimates")
    if (!all(need %in% names(x))) stop("exon-level table is neither DEXSeq nor satuRn; columns: ", paste(names(x), collapse = " "))
    message("Exon-level table read as satuRn (empirical p-values)")
    ids <- tstrsplit(x$ids, ":", fixed = TRUE)
    out <- data.table(gene = ids[[1]], part = ids[[2]], exon_p = x$empirical_pval,
                      exon_padj = x$empirical_FDR, exon_effect = x$estimates)
  }
  out[, gene := gsub('"', "", gene)][, part := gsub('"', "", part)]
  out[]
}
ex <- read_exon_results(opt$exon_results)
message(" - exonic parts in the exon-level table: ", nrow(ex),
        "; groups merging genes: ", uniqueN(grep("+", ex$gene, fixed = TRUE, value = TRUE)))

## --- GrASE ---
g <- rbindlist(lapply(grase_files, fread, colClasses = list(character = c("gene", "event", "setdiff1", "setdiff2"))),
               use.names = TRUE, fill = TRUE)
cts <- unique(g$contrast)
if (is.null(opt$contrast)) {
  if (length(cts) > 1) stop("GrASE results hold ", length(cts), " contrasts; pick one with --contrast: ", paste(cts, collapse = " "))
} else {
  g <- g[contrast == opt$contrast]
  if (!nrow(g)) stop("no rows for --contrast ", opt$contrast, "; the files have: ", paste(cts, collapse = " "))
}
g <- g[!is.na(p.value)]
g[, part := ifelse(grepl("diff1", comparison), setdiff1, setdiff2)]
g <- g[!is.na(part) & part != "NA" & nzchar(part)]
gp <- g[, .(part = trimws(unlist(strsplit(part, ",")))),
        by = .(gene, event, comparison, grase_p = p.value, grase_padj = padj, grase_effect = effect_size)]
setorder(gp, gene, part, grase_padj, grase_p, na.last = TRUE)
gp <- gp[, .SD[1], by = .(gene, part)]
message(" - GrASE parts: ", nrow(gp), "; without an exon-level row: ", nrow(fsetdiff(gp[, .(gene, part)], ex[, .(gene, part)])))

## --- Match ---
m <- merge(gp, ex, by = c("gene", "part"))
m[, `:=`(grase_sig = !is.na(grase_padj) & grase_padj < opt$pval,
         exon_sig  = !is.na(exon_padj)  & exon_padj  < opt$pval)]
m[, category := fifelse(grase_sig & exon_sig, "Both_Sig",
                fifelse(grase_sig, "GrASE_Only", fifelse(exon_sig, "Exon_Only", "Neither_Sig")))]
message("Parts tested by both: ", nrow(m))

cor_log <- function(a, b) cor(-log10(a + 1e-300), -log10(b + 1e-300), use = "complete.obs")
stats_summary <- list(
  total_common_parts   = nrow(m),
  pearson_log_pval     = cor_log(m$grase_p, m$exon_p),
  spearman_pval        = cor(m$grase_p, m$exon_p, method = "spearman", use = "complete.obs"),
  pearson_effect_size  = cor(m$grase_effect, m$exon_effect, use = "complete.obs"),
  spearman_effect_size = cor(m$grase_effect, m$exon_effect, method = "spearman", use = "complete.obs"),
  n_sig_grase = sum(m$grase_sig),
  n_sig_exon  = sum(m$exon_sig),
  n_sig_both  = sum(m$grase_sig & m$exon_sig),
  jaccard     = sum(m$grase_sig & m$exon_sig) / max(1, sum(m$grase_sig | m$exon_sig))
)
message("--- Summary (adjusted p < ", opt$pval, ") ---")
print(stats_summary)
print(table(m$category))

## --- Plots ---
pdf(plot_file, width = 10, height = 5)
pd <- copy(m)[, `:=`(log_grase = -log10(grase_p + 1e-300), log_exon = -log10(exon_p + 1e-300))]
print(ggplot(pd, aes(log_grase, log_exon, color = category)) + geom_point(alpha = 0.3) +
        geom_abline(intercept = 0, slope = 1, col = "red", linetype = "dashed") +
        labs(title = "P-value correlation (-log10, unadjusted)",
             subtitle = paste("Colour: adjusted p <", opt$pval),
             x = "-log10(GrASE p-value)", y = "-log10(exon-level p-value)") + theme_minimal())
print(ggplot(pd, aes(grase_effect, exon_effect, color = category)) + geom_point(alpha = 0.3) +
        labs(title = "Effect size", x = "GrASE effect size", y = "exon-level effect") + theme_minimal())
invisible(dev.off())
message("Plots saved to: ", plot_file)

## --- Save ---
fwrite(m, output_file, sep = "\t")
fwrite(m[category %in% c("GrASE_Only", "Exon_Only")][order(grase_padj)], discrep_file, sep = "\t")
message("Results saved to: ", output_file, "; discrepancies: ", discrep_file)
