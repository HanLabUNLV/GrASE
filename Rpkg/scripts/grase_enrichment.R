#!/usr/bin/env Rscript
# grase_enrichment.R
# GO enrichment for GrASE significant AS genes, controlling for:
#   (1) number of exon-event tests per gene (n_events)
#   (2) gene expression level (total_expr, from baseMean in test file or
#       an external count file supplied via --expr_file)
#
# Both covariates are fitted jointly via logistic regression with natural
# splines: glm(is_sig ~ ns(log(n_events)) + ns(log(total_expr)), binomial).
# Fitted probabilities serve as per-gene weights fed to goseq's Wallenius
# non-central hypergeometric test (same as the nullp() PWF, but multi-variate).
# If --expr_file is omitted, total_expr = sum of baseMean across all events
# for that gene (a reasonable within-file proxy).
#
# Required packages: goseq, org.Hs.eg.db, GO.db, AnnotationDbi,
#                    splines, dplyr, ggplot2, optparse
#
# Usage:
#   Rscript grase_enrichment.R -f test_results.txt -o outdir/
#   Rscript grase_enrichment.R -f test_results.txt -o outdir/ \
#       --contrast CD4_N_STIM_vs_TH1 --ont BP

suppressPackageStartupMessages({
  library(goseq)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(splines)
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-f", "--file"), type="character",
              help="GrASE test result file (output of exontest.R)"),
  make_option(c("-o", "--outdir"), type="character", default=".",
              help="output directory [default: .]"),
  make_option(c("--contrast"), type="character", default=NULL,
              help=paste("restrict to one contrast (e.g. 'CD4_N_STIM_vs_TH1').",
                         "Default: pool all contrasts -- a gene is significant",
                         "if it is significant in ANY contrast.")),
  make_option(c("--expr_file"), type="character", default=NULL,
              help=paste("optional tab-separated file with per-gene expression.",
                         "Two columns: gene (Ensembl ID, with or without version)",
                         "and count (raw or normalized counts). If omitted,",
                         "total_expr = sum(baseMean) across events from --file.")),
  make_option(c("--ont"), type="character", default="BP",
              help="GO ontology: BP, MF, or CC [default: BP]"),
  make_option(c("--alpha"), type="double", default=0.05,
              help="FDR threshold for enrichment results [default: 0.05]"),
  make_option(c("--min_geneset"), type="integer", default=10,
              help="minimum tested genes per GO term [default: 10]"),
  make_option(c("--max_geneset"), type="integer", default=500,
              help="maximum tested genes per GO term [default: 500]"),
  make_option(c("--top_n"), type="integer", default=30,
              help="top N terms to include in the dotplot [default: 30]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$file)) stop("--file is required")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ---------- 1. Load test results and build per-gene summary ----------

cat("Reading", opt$file, "\n")
res <- read.table(opt$file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                  quote="")
required_cols <- c("gene", "event", "significant")
missing <- setdiff(required_cols, names(res))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse=", "))

if (!is.null(opt$contrast)) {
  res <- res[res$contrast == opt$contrast, ]
  if (nrow(res) == 0) stop("No rows match contrast: ", opt$contrast)
  cat("Restricted to contrast:", opt$contrast, "\n")
}

res$gene_base <- sub("\\.[0-9]+$", "", res$gene)

# Per gene: significant if ANY event is significant; n_events = tests done.
# total_expr = mean(baseMean) across events.  Using mean (not sum) keeps
# total_expr orthogonal to n_events; sum grows with n_events and acts as
# a near-duplicate of that covariate (cor ~0.41 vs ~0.13 for mean).
# NOTE: mean(baseMean) reflects read depth / statistical power.  Adding it
# as a covariate corrects for the power bias (highly expressed genes detect
# AS more easily) but also removes real housekeeping GO signal.  Use
# --expr_file with external normalized gene counts if you want to separate
# expression level from within-event coverage depth.
gene_summary <- res %>%
  group_by(gene_base, event) %>%
  summarise(event_mean = mean(baseMean, na.rm=TRUE), .groups="drop") %>%
  group_by(gene_base) %>%
  summarise(
    n_events   = n_distinct(event),
    total_expr = mean(event_mean, na.rm=TRUE),
    .groups    = "drop"
  ) %>%
  left_join(
    res %>%
      group_by(gene_base) %>%
      summarise(is_sig = any(significant == TRUE, na.rm=TRUE), .groups="drop"),
    by="gene_base"
  )

n_tested <- nrow(gene_summary)
n_sig    <- sum(gene_summary$is_sig)
cat(sprintf("Genes tested: %d  |  significant: %d  |  fraction: %.3f\n",
            n_tested, n_sig, n_sig / n_tested))

# ---------- 2. Attach external expression if supplied ----------

if (!is.null(opt$expr_file)) {
  cat("Reading expression file:", opt$expr_file, "\n")
  expr_df <- read.table(opt$expr_file, header=TRUE, sep="\t",
                        stringsAsFactors=FALSE, quote="")
  if (ncol(expr_df) < 2)
    stop("--expr_file must have at least two columns: gene and count")
  colnames(expr_df)[1:2] <- c("gene_base", "total_expr_ext")
  expr_df$gene_base <- sub("\\.[0-9]+$", "", expr_df$gene_base)
  expr_df <- expr_df[!duplicated(expr_df$gene_base), c("gene_base","total_expr_ext")]
  gene_summary <- left_join(gene_summary, expr_df, by="gene_base") %>%
    mutate(total_expr = ifelse(!is.na(total_expr_ext), total_expr_ext, total_expr)) %>%
    select(-total_expr_ext)
  cat("External expression matched to",
      sum(!is.na(gene_summary$total_expr)), "genes\n")
}

# ---------- 3. Map Ensembl -> Entrez ----------

cat("Mapping Ensembl -> Entrez IDs...\n")
id_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys=gene_summary$gene_base,
                        columns="ENTREZID",
                        keytype="ENSEMBL")
)
id_map <- id_map[!is.na(id_map$ENTREZID) & !duplicated(id_map$ENSEMBL), ]
colnames(id_map) <- c("gene_base", "ENTREZID")

gene_summary <- inner_join(gene_summary, id_map, by="gene_base") %>%
  arrange(ENTREZID, desc(n_events)) %>%
  filter(!duplicated(ENTREZID))
cat(sprintf("Genes with Entrez mapping: %d / %d (%.1f%%)\n",
            nrow(gene_summary), n_tested,
            100 * nrow(gene_summary) / n_tested))

# ---------- 4. Fit PWF via logistic regression ----------
# Default: P(sig) ~ ns(log(n_events))          [n_events bias only]
# With --expr_file: P(sig) ~ ns(log(n_events)) + ns(log(total_expr))
#
# We use the GLM route (instead of nullp()) so that the two-covariate
# model is handled consistently.  For the single-covariate case the
# GLM spline fit is nearly identical to nullp()'s monotone spline.
#
# Why not always include expression?
#   mean(baseMean) strongly predicts significance (AIC improvement ~2300)
#   because depth drives power.  Correcting for it removes housekeeping
#   GO signal (translation, splicing) that is genuine biology, not artifact.
#   Use --expr_file with external normalized gene counts to include
#   expression as a covariate only when you have a proper measure of it
#   that is independent of within-event read depth.

gene_summary$log_n <- log(pmax(gene_summary$n_events, 1))

if (!is.null(opt$expr_file)) {
  cat("Fitting P(sig) ~ ns(log(n_events)) + ns(log(total_expr))...\n")
  gene_summary$log_exp <- log(pmax(gene_summary$total_expr, 1))
  glm_fit <- glm(is_sig ~ ns(log_n, df=3) + ns(log_exp, df=3),
                 data=gene_summary, family=binomial())
} else {
  cat("Fitting P(sig) ~ ns(log(n_events))...\n")
  glm_fit <- glm(is_sig ~ ns(log_n, df=3),
                 data=gene_summary, family=binomial())
}
cat("GLM converged:", glm_fit$converged, "\n")
gene_summary$pwf_prob <- fitted(glm_fit)

# Build the data frame that goseq() expects (same structure as nullp() output)
gene_vec <- setNames(as.integer(gene_summary$is_sig), gene_summary$ENTREZID)
pwf <- data.frame(
  DEgenes   = gene_vec,
  bias.data = setNames(gene_summary$n_events,  gene_summary$ENTREZID),
  pwf       = setNames(gene_summary$pwf_prob,  gene_summary$ENTREZID),
  row.names = gene_summary$ENTREZID
)

# ---------- 5. PWF diagnostic plot ----------

pwf_plot_file <- file.path(opt$outdir, "pwf_fit.pdf")
use_expr <- !is.null(opt$expr_file)
pdf(pwf_plot_file, width=if (use_expr) 9 else 5, height=4)
par(mfrow=c(1, if (use_expr) 2 else 1))

# Panel 1: P(sig) vs n_events
n_seq  <- exp(seq(min(gene_summary$log_n), max(gene_summary$log_n), length=200))
nd1    <- data.frame(log_n=log(n_seq))
if (use_expr) nd1$log_exp <- median(gene_summary$log_exp)
pred1  <- predict(glm_fit, newdata=nd1, type="response")
plot(n_seq, pred1, type="l", log="x",
     xlab="n events (log scale)", ylab="P(significant)",
     main=if (use_expr) "Effect of n_events\n(expression at median)" else "P(sig) ~ n_events")
rug(jitter(gene_summary$n_events[gene_summary$is_sig],  amount=0.5), col="red")
rug(jitter(gene_summary$n_events[!gene_summary$is_sig], amount=0.5), col="grey")

# Panel 2 (only when expr covariate used): P(sig) vs total_expr
if (use_expr) {
  exp_seq <- exp(seq(min(gene_summary$log_exp), max(gene_summary$log_exp), length=200))
  nd2     <- data.frame(log_n=median(gene_summary$log_n), log_exp=log(exp_seq))
  pred2   <- predict(glm_fit, newdata=nd2, type="response")
  plot(exp_seq, pred2, type="l", log="x",
       xlab="expression (log scale)", ylab="P(significant)",
       main="Effect of total_expr\n(n_events at median)")
  rug(jitter(gene_summary$total_expr[gene_summary$is_sig],  amount=0.5), col="red")
  rug(jitter(gene_summary$total_expr[!gene_summary$is_sig], amount=0.5), col="grey")
}

dev.off()
cat("PWF plot saved to", pwf_plot_file, "\n")

# ---------- 6. Build GO gene sets from org.Hs.eg.db ----------

cat("Fetching GO annotations (ontology:", opt$ont, ")...\n")
go_all <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys=names(gene_vec),
                        columns=c("GO", "ONTOLOGY"),
                        keytype="ENTREZID")
)
go_all <- go_all[!is.na(go_all$GO) & go_all$ONTOLOGY == opt$ont,
                 c("ENTREZID", "GO")]
go_all <- unique(go_all)
colnames(go_all) <- c("gene_id", "go_id")
cat(sprintf("GO terms with at least one tested gene: %d\n",
            length(unique(go_all$go_id))))

# ---------- 7. Run GOseq (Wallenius) ----------

cat("Running GOseq Wallenius test...\n")
goseq_res <- goseq(pwf, gene2cat=go_all)

# ---------- 8. Annotate, filter, BH-correct ----------

goseq_res$Description <- suppressMessages(
  mapIds(GO.db, keys=goseq_res$category, column="TERM",
         keytype="GOID", multiVals="first")
)

goseq_res <- goseq_res[
  !is.na(goseq_res$numInCat) &
  goseq_res$numInCat >= opt$min_geneset &
  goseq_res$numInCat <= opt$max_geneset, ]

goseq_res$p.adjust  <- p.adjust(goseq_res$over_represented_pvalue, method="BH")
goseq_res$GeneRatio <- goseq_res$numDEInCat / n_sig

goseq_res <- goseq_res %>%
  rename(ID=category, pvalue=over_represented_pvalue,
         Count=numDEInCat, BgSize=numInCat) %>%
  arrange(p.adjust) %>%
  select(ID, Description, Count, BgSize, GeneRatio, pvalue, p.adjust)

# Annotate significant genes per term
sig_entrez <- names(gene_vec)[gene_vec == 1]
sym_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys=sig_entrez, columns="SYMBOL", keytype="ENTREZID")
)
sym_map    <- sym_map[!duplicated(sym_map$ENTREZID), ]
sym_lookup <- setNames(sym_map$SYMBOL, sym_map$ENTREZID)

term_gene_map    <- split(go_all$gene_id, go_all$go_id)
goseq_res$geneID <- sapply(goseq_res$ID, function(term) {
  sig_in <- intersect(sig_entrez, term_gene_map[[term]])
  syms   <- sym_lookup[sig_in]
  paste(syms[!is.na(syms)], collapse="/")
})

# ---------- 9. Write results ----------

out_file <- file.path(opt$outdir, sprintf("enrichment_GO_%s.txt", opt$ont))
write.table(goseq_res, out_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Results written to", out_file, "\n")
cat(sprintf("Significant terms (FDR < %.2f): %d\n",
            opt$alpha, sum(goseq_res$p.adjust < opt$alpha, na.rm=TRUE)))

# ---------- 10. Dotplot ----------

top_terms <- goseq_res %>%
  filter(!is.na(p.adjust), p.adjust < opt$alpha) %>%
  slice_min(p.adjust, n=opt$top_n, with_ties=FALSE)

if (nrow(top_terms) == 0) {
  cat(sprintf("No terms significant at FDR < %.2f -- skipping dotplot.\n",
              opt$alpha))
} else {
  expr_label <- if (!is.null(opt$expr_file)) "external counts" else "mean(baseMean)"
  p <- ggplot(top_terms,
              aes(x=GeneRatio,
                  y=reorder(Description, -log10(p.adjust)),
                  color=p.adjust, size=Count)) +
    geom_point() +
    scale_color_gradient(low="red", high="blue", name="FDR") +
    scale_size_continuous(name="Gene count", range=c(2, 8)) +
    labs(
      x        = sprintf("Gene ratio (sig in term / %d total sig)", n_sig),
      y        = NULL,
      title    = sprintf("GO %s enrichment", opt$ont),
      subtitle = sprintf("GOseq Wallenius; covariates: n_events + total_expr (%s)",
                         expr_label)
    ) +
    theme_bw(base_size=11) +
    theme(axis.text.y=element_text(size=9))

  plot_file <- file.path(opt$outdir,
                         sprintf("enrichment_GO_%s_dotplot.pdf", opt$ont))
  ggsave(plot_file, p, width=10, height=max(4, 0.3 * nrow(top_terms) + 3))
  cat("Dotplot saved to", plot_file, "\n")
}

cat("Done.\n")
