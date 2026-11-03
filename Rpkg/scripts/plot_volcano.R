#!/usr/bin/env Rscript
## Volcano plots of GrASE bipartition results on three effect scales, one panel each:
##   delta_pi          raw-count path proportion difference (pi_trt - pi_ctrl)
##   delta_pi_perbase  the same after length normalization (grase::pi_perbase)
##   effect_size       log odds ratio of the fitted model
## y is -log10 of the unadjusted p-value (--y p.value, default) or of padj (--y padj), uncapped unless
## --ymax; a value of 0 is drawn at the largest finite one. padj is NA for tests in genes not selected at
## the gene level, so --y padj drops those rows. Rows are the input files (e.g. internal, terminal); only
## rows with the chosen value are drawn. Colour is the call of grase::is_significant() at --pval with --min_dpi/--min_dpi_sj/--min_effect/--min_reads.
## Dashed lines: the padj cutoff (--y padj only), +/- --min_dpi on the two delta_pi scales and +/- --min_effect on effect_size.
## Per-base delta_pi already in the input (e.g. from add_perbase_calls.R) is used as is; otherwise it is computed
## when --gff_dir is given.
## Usage: Rscript plot_volcano.R --grase <internal.annotated.txt>,<terminal.annotated.txt> --labels internal,terminal
##          --contrast <name> --gff_dir <per-gene DEXSeq GFF dir> --out <prefix>
## Output: <prefix>.png and <prefix>.pdf (points rasterized)
suppressPackageStartupMessages({library(optparse); library(data.table); library(ggplot2); library(grase)})
o <- parse_args(OptionParser(option_list = list(
  make_option("--grase",      type = "character", help = "annotated exontest.R results, comma separated"),
  make_option("--labels",     type = "character", default = NULL, help = "row labels, comma separated [file names]"),
  make_option("--contrast",   type = "character", default = NULL, help = "keep only this contrast"),
  make_option("--gff_dir",    type = "character", default = NULL, help = "per-gene DEXSeq GFFs, for delta_pi_perbase"),
  make_option("--read_length", type = "integer", default = NULL, help = "read length, for the effective lengths with --gff_dir"),
  make_option("--out",        type = "character", help = "output prefix"),
  make_option("--pval",       type = "double", default = 0.01, help = "adjusted p cutoff [%default]"),
  make_option("--min_dpi",    type = "double", default = 0.05, help = "|delta_pi| of is_significant(), per-base where defined [%default]"),
  make_option("--min_dpi_sj", type = "double", default = 0.05, help = "|delta_pi| for a junction-sourced side [%default]"),
  make_option("--min_reads",  type = "double", default = 10,   help = "reads on the distinct set in a group [%default]"),
  make_option("--min_effect", type = "double", default = 0.5,  help = "|effect_size| of is_significant() [%default]"),
  make_option("--xmax_effect", type = "double", default = 10,  help = "|effect_size| beyond this is drawn at the edge [%default]"),
  make_option("--y",          type = "character", default = "p.value", help = "p.value or padj [%default]"),
  make_option("--ymax",       type = "double", default = Inf,  help = "-log10(padj) beyond this is drawn at the edge [no cap]"),
  make_option("--font",       type = "character", default = "Arial"))))
for (r in c("grase", "out")) if (is.null(o[[r]])) stop("--", r, " is required")
files <- strsplit(o$grase, ",")[[1]]
labs  <- if (is.null(o$labels)) basename(files) else strsplit(o$labels, ",")[[1]]
stopifnot(length(labs) == length(files))

stopifnot(o$y %in% c("p.value", "padj"))
cols <- c("contrast", "gene", "comparison", "p.value", "padj", "delta_pi", "effect_size", "pi_ctrl", "pi_trt",
          "lfc_diff_net", "d_support", "ref_ex_part", "setdiff1", "setdiff2")
D <- rbindlist(lapply(seq_along(files), function(i) {
  hd <- names(fread(files[i], nrows = 0))
  d <- fread(files[i], select = intersect(c(cols, "pi_ref", "delta_pi_perbase"), hd))
  if ("pi_ref" %in% names(d)) setnames(d, "pi_ref", "pi_ctrl")   # results written before the rename
  if (!is.null(o$contrast)) d <- d[contrast == o$contrast]
  if (!"delta_pi_perbase" %in% names(d) && !is.null(o$gff_dir))
    d <- as.data.table(annotate_pi_perbase(as.data.frame(d), path.expand(o$gff_dir), o$read_length))
  if (!"delta_pi_perbase" %in% names(d)) d[, delta_pi_perbase := NA_real_]
  ## the call reads delta_pi_perbase, so it is made after the column exists
  d[, called := is_significant(d, o$pval, 0, o$min_dpi, o$min_dpi_sj, o$min_reads, min_effect = o$min_effect)]
  d <- d[!is.na(get(o$y))]
  d[, class := labs[i]]
  d
}))
cat(sprintf("%s: %d rows with %s, %d called\n", labs, D[, .N, by = class]$N, o$y, D[, sum(called), by = class]$V1), sep = "")

meas <- c(delta_pi = "delta_pi", delta_pi_perbase = "per-base delta_pi", effect_size = "log odds ratio (effect_size)")
L <- rbindlist(lapply(names(meas), function(m)
  D[, .(class, called, x = get(m), y = -log10(get(o$y)), measure = meas[[m]])]))
## a p of 0 underflowed: draw it at the largest finite value
ytop <- max(L$y[is.finite(L$y)])
n0 <- sum(!is.finite(L$y)) / length(meas)
L[!is.finite(y), y := ytop]
L[, y := pmin(y, o$ymax)]
L <- L[!is.na(x)]
L[measure == meas[["effect_size"]], x := pmax(pmin(x, o$xmax_effect), -o$xmax_effect)]
L[, measure := factor(measure, meas)]
L[, class := factor(class, labs)]
L[, Call := factor(called, c(FALSE, TRUE), c("not called", "called"))]
n_pb <- D[, sum(!is.na(delta_pi_perbase)), by = class]
cat("rows with per-base delta_pi:", paste(n_pb$class, n_pb$V1, collapse = ", "), "\n")

V <- data.table(measure = factor(meas, meas), v = c(o$min_dpi, o$min_dpi, o$min_effect))
p <- ggplot(L[order(called)], aes(x, y, colour = Call)) +
  ggrastr::rasterise(geom_point(size = 0.35, alpha = 0.5), dpi = 300) +
  (if (o$y == "padj") geom_hline(yintercept = -log10(o$pval), linetype = 2, linewidth = 0.3)) +
  geom_vline(data = V, aes(xintercept = v), linetype = 2, linewidth = 0.3) +
  geom_vline(data = V, aes(xintercept = -v), linetype = 2, linewidth = 0.3) +
  facet_grid(class ~ measure, scales = "free_x") +
  scale_colour_manual(values = c("not called" = "#b9c3cf", "called" = "#2a78d6"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = NULL, y = paste0("-log10(", o$y, ")", if (is.finite(o$ymax)) sprintf(", capped at %g", o$ymax) else "",
                            if (n0 > 0) sprintf("; 0 drawn at %.0f", min(ytop, o$ymax)) else ""),
       title = if (is.null(o$contrast)) NULL else o$contrast) +
  theme_bw(base_size = 11, base_family = o$font) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
h <- 1.2 + 2.8 * length(labs)
ggsave(paste0(o$out, ".png"), p, width = 10, height = h, dpi = 200, device = ragg::agg_png)
ggsave(paste0(o$out, ".pdf"), p, width = 10, height = h, device = cairo_pdf)
cat("wrote", paste0(o$out, c(".png", ".pdf")), "\n")
