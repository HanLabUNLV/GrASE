#!/usr/bin/env Rscript
## Turn enumerated outward units into a count table exontest.R can run.
##
## INPUT is enumerate_outward_units.R's anchor table: one row per
## (side, anchor), where an anchor is either a TERMINUS of that side (one per
## route, mean 3.03 per side) or D's outer edge, and `unit_id` is the tiling unit
## immediately outward of it.
##
## ONE TEST PER (gene, unit, ref_ex_part), NOT PER ANCHOR. 72,470 anchors
## collapse to 39,911 test units (44.9%), because different sides and events
## reach the same unit against the same reference -- testing those separately
## would be the same measurement repeated under different labels.
##
## THE COLLAPSE IS EXACT, not approximate. `ref` is the counts of ref_ex_part,
## so two sides sharing a reference part have identical reference counts:
## verified on real data, 0 of 85,340 (gene, ref_ex_part, sample) groups spanning
## more than one event had an inconsistent `ref`. That is also why the reference
## pass keys on (gene, ref_ex_part, sample) rather than (gene, event, sample),
## which cuts it from 403,217 gene/event pairs to 125,407.
##
## UNIT COUNTS COME FROM TWO PLACES and are never recounted: intronic units from
## the intronic pass, exonic units from the per-part DEXSeq counts the exonic arm
## already used. Quantification is universal -- each unit is counted once, in
## whichever annotation owns it.
##
## SEPARATE exontest RUN. Sharing one run with the exonic arm would change
## nested_BH's within-gene burden and move every exonic padj. Synthetic event
## ids are therefore safe: nothing collides, and the side table maps them back.
##
## ref IS UNAFFECTED BY THE JUNCTION RULE -- it is always an exonic part -- so
## the pre-2026-10-04 combined table is valid as a reference source. What the
## junction rerun will change is which sides are CALLED, i.e. the enumeration.
##
## Usage:
##   Rscript scripts/build_outward_counts.R --units intronparts/outward_units.tsv \
##     --intron_counts intronparts/counts --exon_counts count_files \
##     --exon_combined <...exoncnt.combined.txt> \
##     --out intronparts/outward.exoncnt.combined.txt \
##     --out_sides intronparts/outward_side_table.tsv [--ref_cache <file>]
suppressPackageStartupMessages({library(data.table); library(optparse)})
a <- parse_args(OptionParser(option_list = list(
  make_option("--units",         type = "character"),
  make_option("--intron_counts", type = "character"),
  make_option("--exon_counts",   type = "character"),
  make_option("--exon_combined", type = "character"),
  make_option("--out",           type = "character"),
  make_option("--out_sides",     type = "character"),
  make_option("--ref_cache",     type = "character", default = NULL,
              help = "reuse/populate a (gene, ref_ex_part, sample) reference selection"))))
for (r in c("units","intron_counts","exon_counts","exon_combined","out","out_sides"))
  if (is.null(a[[r]])) stop("--", r, " is required")

u <- fread(a$units, sep = "\t", showProgress = FALSE,
           colClasses = list(character = c("gene","event")))
cat(sprintf("anchor rows: %d\n", nrow(u)))

## ---- test units ------------------------------------------------------------
tu <- unique(u[, .(gene, unit_type, unit_id, ref_ex_part)])
setorder(tu, gene, unit_id, ref_ex_part)
tu[, synth_event := as.character(seq_len(.N)), by = gene]
## ONE TEST UNIT PER ROW. An earlier version packed two units sharing a
## (gene, ref_ex_part) into the diff1/diff2 slots, because exontest.R called
## compute_lfc_summary on both arms unconditionally and died on an all-NA diff2.
## That is fixed in the package (the zero-row guard in compute_lfc_summary), so
## the slots no longer have to be filled. Pairing was also a hazard in its own
## right: it put unrelated units in the two slots of one `event`, which would
## make any per-event combined p-value (--combine_p) meaningless.
cat(sprintf("test units: %d over %d genes  (collapse %.1f%% from anchors)\n",
            nrow(tu), uniqueN(tu$gene), 100*(1 - nrow(tu)/nrow(u))))
u <- merge(u, tu[, .(gene, unit_type, unit_id, ref_ex_part, synth_event)],
           by = c("gene","unit_type","unit_id","ref_ex_part"), all.x = TRUE)
if (anyNA(u$synth_event)) stop("anchor row failed to map to a test unit")

## ---- unit counts, two sources ---------------------------------------------
## FILTER ON THE RAW KEY, then transform only the survivors. The count files
## write `"GENE":"UNIT"`; the earlier version stripped the quotes and swapped
## the colon on EVERY row of every file and only then filtered, which on the
## lineage panel meant 757M rows transformed (699 samples x 1.08M rows) to keep
## 67M (8.8%) -- 91% of the string work discarded. Building the wanted keys in
## raw form once costs 95k paste0 calls instead.
##
## The old `!startsWith(k, "_")` guard is gone with it: HTSeq's `_ambiguous`,
## `_no_feature` etc. carry no quotes, so they can never match a raw key.
raw_key <- function(gene, unit) paste0('"', gene, '":"', unit, '"')
read_counts <- function(fn, keep_raw) {
  d <- fread(fn, header = FALSE, sep = "\t", quote = "", showProgress = FALSE,
             col.names = c("k","n"))
  d[k %chin% keep_raw, .(k, n = as.integer(n))]
}
want_i <- unique(tu[unit_type == "intron", raw_key(gene, unit_id)])
want_e <- unique(tu[unit_type == "exon",   raw_key(gene, unit_id)])
cat(sprintf("units to fetch: %d intronic, %d exonic\n", length(want_i), length(want_e)))

ifiles  <- sort(list.files(a$intron_counts, pattern = "_counts\\.txt$", full.names = TRUE))
if (!length(ifiles)) stop("no *_counts.txt under ", a$intron_counts)
samples <- sub("_counts\\.txt$", "", basename(ifiles))
cat(sprintf("samples: %d\n", length(samples)))
## Parallel over FILES: each is independent, and the serial lapply read 1,398
## of them one at a time on a many-core box.
ncore <- as.integer(Sys.getenv("MC_CORES", unset = "16"))
cat(sprintf("reading count files on %d cores\n", ncore))
ic <- rbindlist(parallel::mclapply(seq_along(ifiles), function(i)
  read_counts(ifiles[i], want_i)[, sample := samples[i]], mc.cores = ncore))

ex_path <- function(s) file.path(a$exon_counts, sub("_SRR[0-9]+$", "", s),
                                 paste0(sub("^.*_(SRR[0-9]+)$", "\\1", s), "_counts.txt"))
ep <- vapply(samples, ex_path, "")
if (any(!file.exists(ep)))
  stop(sprintf("no exonic count file for %d samples, e.g. %s", sum(!file.exists(ep)),
               paste(head(samples[!file.exists(ep)], 3), collapse = ", ")))
ec <- rbindlist(parallel::mclapply(seq_along(ep), function(i)
  read_counts(ep[i], want_e)[, sample := samples[i]], mc.cores = ncore))
cnt <- rbind(ic, ec)
## now, and only now, put the key in `GENE|UNIT` form -- once, on the kept rows
cnt[, k := gsub('"', "", k, fixed = TRUE)]
cnt[, k := sub(":", "|", k, fixed = TRUE)]
cat(sprintf("count rows: %d over %d units\n", nrow(cnt), uniqueN(cnt$k)))

## An unresolved unit becomes 0 downstream, not an error, and a zero unit is
## exactly the boundary verdict -- a silently missing unit would manufacture the
## very result this arm measures.
need <- unique(tu[, .(k = paste(gene, unit_id, sep = "|"))])
missing <- need[!k %chin% unique(cnt$k)]
if (nrow(missing))
  stop(sprintf(paste0("FATAL: %d units absent from the count matrices, e.g. %s\n",
                      "A missing unit counts as zero, indistinguishable from a real\n",
                      "boundary. Fix the annotation/counting before continuing."),
               nrow(missing), paste(head(missing$k, 5), collapse = ", ")))
cat(sprintf("all %d units resolve to counts\n", nrow(need)))

## ---- reference counts, keyed on (gene, ref_ex_part, sample) ----------------
con   <- file(a$exon_combined, "r")
hdr   <- strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]]
probe <- strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]]
close(con)
shift <- length(probe) - length(hdr)
if (!shift %in% c(0L,1L)) stop("unexpected combined-table shape")
idx <- setNames(seq_along(hdr) + shift, hdr); ncol_ <- length(hdr) + shift
for (c_ in c("gene","event","ref_ex_part","sample","ref","diff1","diff2","groups"))
  if (!c_ %in% names(idx)) stop("combined table lacks column '", c_, "'")

## THE CACHE IS CHECKED FOR COVERAGE, NOT TRUSTED. Reusing it blindly dropped
## 708 units silently when the unit set grew: the assemble merge is an inner
## join, so a (gene, ref_ex_part) absent from the cache removes that unit with
## no error. It happened because sides whose every anchor sat at a gene edge
## were dropped before the cache was built, and adding terminus_self units
## brought them back. Missing pairs are fetched and appended here.
select_refs <- function(pairs) {
  keyf <- tempfile(); sel <- tempfile()
  fwrite(pairs[, .(k = paste(gene, ref_ex_part, sep = "\t"))], keyf,
         col.names = FALSE, quote = FALSE)
  ## LC_ALL=C: the keys are ASCII, and locale-aware comparison roughly halves
  ## awk's throughput on a 56 GiB single pass.
  cmd <- sprintf(
    "LC_ALL=C awk -F'\\t' -v OFS='\\t' 'NR==FNR{keep[$1\"\\t\"$2]=1; next} FNR==1{next} NF==%d && (($%d\"\\t\"$%d) in keep){print $%d,$%d,$%d,$%d,$%d}' %s %s | sort -u > %s",
    ncol_, idx[["gene"]], idx[["ref_ex_part"]], idx[["gene"]], idx[["ref_ex_part"]],
    idx[["sample"]], idx[["ref"]], idx[["groups"]], shQuote(keyf),
    shQuote(a$exon_combined), shQuote(sel))
  if (system(cmd) != 0) stop("reference selection failed")
  d <- fread(sel, header = FALSE, sep = "\t", showProgress = FALSE,
             colClasses = list(character = 1:3),
             col.names = c("gene","ref_ex_part","sample","ref","groups"))
  d
}

need_pairs <- unique(tu[, .(gene, ref_ex_part)])
refsel <- if (!is.null(a$ref_cache) && file.exists(a$ref_cache)) {
  cat("reading reference cache ", a$ref_cache, "\n", sep = "")
  d <- fread(a$ref_cache, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene","ref_ex_part","sample")))
  miss <- need_pairs[!unique(d[, .(gene, ref_ex_part)]), on = .(gene, ref_ex_part)]
  cat(sprintf("  cache covers %d of %d needed (gene, ref) pairs\n",
              nrow(need_pairs) - nrow(miss), nrow(need_pairs)))
  if (nrow(miss)) {
    cat(sprintf("  fetching %d missing pairs, appending to the cache\n", nrow(miss)))
    d <- rbind(d, select_refs(miss))
    fwrite(d, a$ref_cache, sep = "\t", quote = FALSE)
  }
  d
} else {
  d <- select_refs(need_pairs)
  if (!is.null(a$ref_cache)) { fwrite(d, a$ref_cache, sep = "\t", quote = FALSE)
    cat("wrote reference cache ", a$ref_cache, "\n", sep = "") }
  d
}
refsel <- refsel[sample %chin% samples]
dup <- refsel[, .N, by = .(gene, ref_ex_part, sample)][N > 1]
if (nrow(dup))
  stop(sprintf("FATAL: ref not unique per (gene, ref_ex_part, sample) for %d keys -- the collapse would be wrong", nrow(dup)))
cat(sprintf("reference rows: %d over %d (gene, ref) pairs\n",
            nrow(refsel), uniqueN(refsel[, paste(gene, ref_ex_part)])))
## the assemble merge is an inner join: a unit whose (gene, ref) has no
## reference rows vanishes without error, so require full coverage first
gap <- need_pairs[!unique(refsel[, .(gene, ref_ex_part)]), on = .(gene, ref_ex_part)]
if (nrow(gap))
  stop(sprintf(paste0("FATAL: %d (gene, ref_ex_part) pairs have no reference rows, e.g. %s\n",
                      "Those units would be dropped silently by the join."),
               nrow(gap), paste(head(gap[, paste(gene, ref_ex_part)], 3), collapse = ", ")))

## ---- assemble ---------------------------------------------------------------
tu[, k := paste(gene, unit_id, sep = "|")]
rows <- merge(tu[, .(gene, synth_event, unit_type, ref_ex_part, k)], cnt,
              by = "k", allow.cartesian = TRUE)
rows <- merge(rows, refsel, by = c("gene","ref_ex_part","sample"))
cat(sprintf("emitted rows: %d over %d test units, %d samples\n",
            nrow(rows), uniqueN(rows[, paste(gene, synth_event)]), uniqueN(rows$sample)))

o <- as.data.table(matrix("NA", nrow = nrow(rows), ncol = ncol_))
if (shift) o[[1]] <- as.character(seq_len(nrow(rows)))
o[[idx[["gene"]]]]        <- rows$gene
o[[idx[["event"]]]]       <- rows$synth_event
o[[idx[["sample"]]]]      <- rows$sample
o[[idx[["ref"]]]]         <- as.character(rows$ref)
o[[idx[["diff1"]]]]       <- as.character(rows$n)
o[[idx[["diff2"]]]]       <- "NA"   # one-sided arm; compute_lfc_summary guards it
o[[idx[["groups"]]]]      <- rows$groups
o[[idx[["ref_ex_part"]]]] <- rows$ref_ex_part
if ("diff1_source" %in% names(idx)) o[[idx[["diff1_source"]]]] <- rows$unit_type
con <- file(a$out, "w"); writeLines(paste(hdr, collapse = "\t"), con); close(con)
fwrite(o, a$out, sep = "\t", col.names = FALSE, quote = FALSE, append = TRUE)
cat("wrote ", a$out, "\n", sep = "")

## ---- side table: the only record of what each synthetic event tested -------
## exontest.R consumes diff1_source and does not emit it, and its annotation
## columns (source/sink/setdiff/path) come from the EXONIC per-gene split files,
## so in this arm they describe D and not the unit tested. Read the unit here.
um <- rows[, .(unit_mean = mean(n), unit_sum = sum(n), n_samples = .N),
           by = .(gene, synth_event)]
st <- merge(u[, .(gene, event, side, kind, strand, anchor, anchor_kind, gene_edge,
                  unit_type, unit_id, ref_ex_part, synth_event)],
            um, by = c("gene","synth_event"), all.x = TRUE)
st[, comparison := "diff1_vs_ref"]
st[, unit_mean := sprintf("%.6g", as.numeric(unit_mean))]
setorder(st, gene, event, side, anchor_kind, anchor)
fwrite(st, a$out_sides, sep = "\t", quote = FALSE)
cat("wrote ", a$out_sides, "  (", nrow(st), " anchor rows)\n", sep = "")
uu <- unique(st[, .(gene, synth_event, unit_mean)])
cat(sprintf("  units with mean < 0.5 reads: %d of %d\n",
            sum(as.numeric(uu$unit_mean) < 0.5, na.rm = TRUE), nrow(uu)))
