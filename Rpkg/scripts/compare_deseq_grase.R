library(tidyverse)
library(optparse)

option_list <- list(
  make_option("--deseq_dir", type="character",
              help="directory containing deseq2_<contrast>.txt files"),
  make_option("--grase_dir", type="character",
              help="directory containing GrASE mincomb.annotated.txt files"),
  make_option("--outdir", type="character", default=NULL,
              help="output directory for result tables [default: deseq_dir]"),
  make_option("--padj", type="double", default=0.05,
              help="FDR threshold [default: 0.05]"),
  make_option("--contrasts", type="character", default=NULL,
              help="comma-separated list of contrasts; default: auto-detect from deseq_dir"),
  make_option("--dge_prefix", type="character", default="deseq2_",
              help="filename prefix for DGE classification files [default: deseq2_]")
)

opt <- parse_args(OptionParser(option_list=option_list))

DESEQ_DIR  <- path.expand(opt$deseq_dir)
GRASE_DIR  <- path.expand(opt$grase_dir)
OUTDIR     <- if (is.null(opt$outdir)) DESEQ_DIR else path.expand(opt$outdir)
PADJ_THR   <- opt$padj
DGE_PREFIX <- opt$dge_prefix
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)

strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

if (!is.null(opt$contrasts)) {
  contrasts <- trimws(strsplit(opt$contrasts, ",")[[1]])
} else {
  pat <- paste0("^", DGE_PREFIX, ".+\\.txt$")
  contrasts <- sub(paste0("^", DGE_PREFIX, "(.+)\\.txt$"), "\\1",
                   list.files(DESEQ_DIR, pattern=pat))
}
message("Contrasts: ", paste(contrasts, collapse=", "))

grase_int <- read.table(
  file.path(GRASE_DIR, "test_bipartition.internal_betabinom_EBmap.annotated.txt"),
  header=TRUE, sep="\t", stringsAsFactors=FALSE)
grase_tss <- read.table(
  file.path(GRASE_DIR, "test_bipartition.TSSTTS_betabinom_EBmap.annotated.txt"),
  header=TRUE, sep="\t", stringsAsFactors=FALSE)

summary_rows <- list()

for (ct in contrasts) {
  deseq  <- read.table(file.path(DESEQ_DIR, paste0(DGE_PREFIX, ct, ".txt")),
                       sep="\t", header=TRUE, row.names=1)
  d_sig  <- strip_ver(rownames(deseq)[!is.na(deseq$padj) & deseq$padj < PADJ_THR])
  gi_sig <- strip_ver(unique(grase_int$gene[grase_int$contrast == ct & grase_int$significant]))
  gt_sig <- strip_ver(unique(grase_tss$gene[ grase_tss$contrast == ct & grase_tss$significant]))

  # primary biological split: DGE vs no-DGE
  g_any_sig    <- union(gi_sig, gt_sig)
  DGE_only     <- setdiff(d_sig,     g_any_sig)
  splicing_only <- setdiff(g_any_sig, d_sig)
  both         <- intersect(d_sig,   g_any_sig)

  # within splicing_only: which assay(s) detected it
  sp_i_only <- setdiff(gi_sig, union(d_sig, gt_sig))
  sp_t_only <- setdiff(gt_sig, union(d_sig, gi_sig))
  sp_i_and_t <- setdiff(intersect(gi_sig, gt_sig), d_sig)

  # within both: which assay(s) detected splicing
  both_i_only <- setdiff(intersect(d_sig, gi_sig), gt_sig)
  both_t_only <- setdiff(intersect(d_sig, gt_sig), gi_sig)
  both_i_and_t <- intersect(d_sig, intersect(gi_sig, gt_sig))

  categories <- list(
    DGE_only      = DGE_only,
    splicing_only = splicing_only,
    both          = both,
    # secondary breakdown -- splicing_only by assay
    sp_i_only     = sp_i_only,
    sp_t_only     = sp_t_only,
    sp_i_and_t    = sp_i_and_t,
    # secondary breakdown -- both by assay
    both_i_only   = both_i_only,
    both_t_only   = both_t_only,
    both_i_and_t  = both_i_and_t
  )

  for (cat_name in names(categories)) {
    genes <- categories[[cat_name]]
    outfile <- file.path(OUTDIR, paste0(ct, ".", cat_name, ".genes.txt"))
    writeLines(genes, outfile)
  }

  summary_rows[[ct]] <- data.frame(
    contrast      = ct,
    DESeq2        = length(d_sig),
    GrASE_i       = length(gi_sig),
    GrASE_t       = length(gt_sig),
    DGE_only      = length(DGE_only),
    splicing_only = length(splicing_only),
    both          = length(both),
    sp_i_only     = length(sp_i_only),
    sp_t_only     = length(sp_t_only),
    sp_i_and_t    = length(sp_i_and_t),
    both_i_only   = length(both_i_only),
    both_t_only   = length(both_t_only),
    both_i_and_t  = length(both_i_and_t)
  )
  message(ct, ": done")
}

summary_df <- bind_rows(summary_rows)
print(summary_df)

outfile <- file.path(OUTDIR, "deseq_grase_overlap_summary.txt")
write.table(summary_df, outfile, sep="\t", quote=FALSE, row.names=FALSE)
message("Summary written to ", outfile)
message("Per-contrast gene lists written to ", OUTDIR)
