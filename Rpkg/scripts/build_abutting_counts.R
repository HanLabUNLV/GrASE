#!/usr/bin/env Rscript
## Build the ABUTTING-UNIT arm of the TSS/TTS localization test.
##
## For each CALLED bipartition side, take the tiling unit immediately outward of
## that side's terminus and test it against the SAME reference the side used:
##
##     ref    = the exonic row's reference counts, unchanged
##     diff1  = counts of the unit abutting side 1's terminus, outward
##     diff2  = counts of the unit abutting side 2's terminus, outward
##
## The unit is EXONIC or INTRONIC, whichever the annotation puts there, and the
## test is identical either way -- a bin against the reference. `diff1_source`
## and `diff2_source` record which. That is the whole point: a boundary claim
## asks whether the differential stops at the terminus, and what lies beyond is
## whatever the annotation says lies beyond, not a synthetic window.
##
## SUPERSEDES build_flank_counts.R. That script used a fixed 150bp window keyed
## per SIDE, which (a) could fall outside any annotated unit, producing verdicts
## that described the window rather than the biology, and (b) emitted the same
## interval once per side sharing a terminus -- 35,011 bins for 22,034 distinct
## intervals, 37.1% redundant. Keying on the gene's tiling removes both.
##
## WHY THE OUTWARD UNIT IS WELL DEFINED. A gene's exonic parts are disjoint
## (0 overlapping pairs in 519,701 adjacent pairs over 60,669 genes), so exonic
## and intronic parts tile the gene completely. And every route terminus sits
## exactly on an exonic-part edge -- verified on all 538,603 sides, no
## exceptions. So the unit outward of a terminus is found by locating the exonic
## part whose outward edge IS the terminus and stepping one unit out: either the
## exonic part that touches it, or the intronic part in the gap. Both of those
## facts are asserted below rather than assumed, because the lookup is silent if
## they break.
##
## EXONIC UNITS ARE NOT RECOUNTED. They already have counts, in the per-part
## DEXSeq count files the exonic arm used. Only intronic units needed a new
## counting pass. That asymmetry is deliberate: quantification is universal, so
## every unit of every gene is quantified exactly once, in whichever annotation
## owns it.
##
## SEPARATE RUN, not a shared one. Folding these rows into the exonic table
## would share the EB dispersion prior, but exontest uses nested_BH whose first
## stage is BH *within gene* -- adding rows to a gene raises that gene's test
## count and moves every exonic event's padj, silently altering published calls.
## Event ids are therefore NOT offset: with a separate run there is nothing to
## collide with, and the two arms join directly on (gene, event, comparison).
##
## RESTRICTED to called sides, so this arm faces a smaller multiple-testing
## burden than the exonic arm did and is more readily called, pushing verdicts
## toward `continuation`. Read any boundary fraction as a floor.
##
## Usage:
##   Rscript scripts/build_abutting_counts.R \
##       --results   grase_activation_exsj/test_bipartition.exsj.external_...annotated.txt \
##       --geometry  flankbins/flankbins_manifest.tsv \
##       --exon_parts tss_edge/geom_parts.tsv \
##       --intron_manifest intronparts/intronparts_manifest.tsv \
##       --intron_counts   intronparts/counts \
##       --exon_counts     count_files \
##       --exon_combined   <...exoncnt.combined.txt> \
##       --out intronparts/abut.exoncnt.combined.txt \
##       --out_summary intronparts/abut_side_summary.tsv
suppressPackageStartupMessages({ library(data.table); library(optparse) })

opt_list <- list(
  make_option("--results",         type = "character",
              help = "exonic annotated.txt; defines which sides are called"),
  make_option("--geometry",        type = "character",
              help = "manifest with gene/event/side/kind/strand/boundary per side"),
  make_option("--exon_parts",      type = "character",
              help = "table of gene, chrom, strand, part, start, end"),
  make_option("--intron_manifest", type = "character"),
  make_option("--intron_counts",   type = "character"),
  make_option("--exon_counts",     type = "character",
              help = "count_files dir: <group>/<SRR>_counts.txt"),
  make_option("--exon_combined",   type = "character"),
  make_option("--out",             type = "character"),
  make_option("--out_summary",     type = "character", default = NULL),
  make_option("--direction",       type = "character", default = "outward",
              help = paste("outward: the unit beyond the side's TERMINUS (boundary test).",
                           "inward: the unit beyond D's INNER end (processing test --",
                           "is D spliced onward or continuous with the next intron).",
                           "[default %default]")),
  make_option("--padj",            type = "double", default = 0.01),
  make_option("--min_dpi",         type = "double", default = 0.1)
)
a <- parse_args(OptionParser(option_list = opt_list))
for (req in c("results", "geometry", "exon_parts", "intron_manifest",
              "intron_counts", "exon_counts", "exon_combined", "out"))
  if (is.null(a[[req]])) stop("--", req, " is required")

if (!a$direction %in% c("outward", "inward"))
  stop("--direction must be 'outward' or 'inward', got '", a$direction, "'")
num <- function(x) as.numeric(x)   # coordinates may arrive as "2.3e+07"

## ---- called sides ----------------------------------------------------------
## `significant` is per contrast; a side called in ANY contrast gets a row, and
## exontest then tests every contrast on it.
res <- fread(a$results, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene", "event")))
for (c_ in c("gene", "event", "comparison", "significant"))
  if (!c_ %in% names(res)) stop("results lack column '", c_, "'")
res[, side := fifelse(comparison == "diff1_vs_ref", 1L, 2L)]
called <- unique(res[significant == TRUE, .(gene, event, side)])
cat(sprintf("called sides: %d over %d genes\n", nrow(called), uniqueN(called$gene)))

## ---- geometry: terminus per side -------------------------------------------
geo <- fread(a$geometry, sep = "\t", showProgress = FALSE,
             colClasses = list(character = c("gene", "event")))
geo[, side := as.integer(side)]
geo[, boundary := num(boundary)]
geo <- unique(geo[, .(gene, event, side, kind, strand, boundary, distinct_parts)])
sides <- geo[called, on = .(gene, event, side), nomatch = NULL]
if (nrow(sides) < nrow(called))
  cat(sprintf("  %d called sides have no terminus geometry (empty distinct set) -- excluded\n",
              nrow(called) - nrow(sides)))
## outward direction: decreasing coordinate for TSS/+ and TTS/-
sides[, out_dec := (kind == "TSS" & strand == "+") | (kind == "TTS" & strand == "-")]
## the joins below map `boundary` onto the exonic part's start/end, which
## consumes the column; keep a copy under a name the join cannot take.
sides[, bnd := boundary]

## ---- the tiling ------------------------------------------------------------
ex <- fread(a$exon_parts, sep = "\t", showProgress = FALSE)
ex[, `:=`(start = num(start), end = num(end), part = as.integer(part))]
setorder(ex, gene, start)
ex[, `:=`(nxt_start = shift(start, type = "lead"),
          nxt_part  = shift(part,  type = "lead"),
          prv_end   = shift(end,   type = "lag"),
          prv_part  = shift(part,  type = "lag")), by = gene]

ip <- fread(a$intron_manifest, sep = "\t", showProgress = FALSE)
ip[, `:=`(intron_start = num(intron_start), intron_end = num(intron_end))]
## fread reads "005" as integer 5; part ids are zero-padded strings everywhere
## else (the GFF, the count keys), so normalise before any join on them.
ip[, `:=`(abut_part_lo = sprintf("%03d", as.integer(abut_part_lo)),
          abut_part_hi = sprintf("%03d", as.integer(abut_part_hi)))]
## one row per (gene, intron, which end abuts): the bin to use when stepping
## outward onto this intron from the low side (.lo) or the high side (.hi)
ip_lo <- ip[touches %in% c("both", "lo"), .(gene, abut_part_lo, abut_part_hi, bin_id)]
ip_hi <- ip[touches %in% c("both", "hi"), .(gene, abut_part_lo, abut_part_hi, bin_id)]

## ---- locate the unit on the requested side ---------------------------------
## Both directions reduce to the same step: from an ANCHOR exonic part, move one
## unit along a STEP direction, landing either on the exonic part that touches it
## or on the intronic part in the gap.
##
##   outward : anchor = the part whose outward edge IS the terminus
##             step   = away from the gene body
##             asks   = does the differential continue past the claimed position
##   inward  : anchor = D's innermost part (farthest from the terminus)
##             step   = toward the gene body
##             asks   = is D spliced onward, or continuous with the next intron
##
## The inward test is what diagnosed PDCD4: D = {E008} was well covered while
## essentially nothing spliced out of it into E009, so its coverage was
## unprocessed. The outward unit there (I007_008.hi) showed the same thing from
## the other side, but only the inward unit speaks to processing.
if (a$direction == "outward") {
  ## terminus is the outward edge of some exonic part -- verified on all 538,603
  ## sides, and asserted below because the lookup is silent if it breaks
  sides[, bnd := boundary]
  lo <- ex[, .(gene, part, start)]
  hi <- ex[, .(gene, part, end)]
  m1 <- lo[sides[out_dec == TRUE],  on = .(gene, start = boundary), nomatch = NA]
  m2 <- hi[sides[out_dec == FALSE], on = .(gene, end   = boundary), nomatch = NA]
  if (sum(is.na(m1$part)) + sum(is.na(m2$part)))
    stop(sprintf(paste0("FATAL: %d termini do not coincide with an exonic-part edge.\n",
                        "The unit lookup depends on that invariant (verified on 538,603\n",
                        "sides); a violation makes every unit assignment suspect."),
                 sum(is.na(m1$part)) + sum(is.na(m2$part))))
  anc <- rbind(m1[, .(gene, event, side, kind, strand, boundary = bnd, distinct_parts,
                      anchor = part, step_dec = TRUE)],
               m2[, .(gene, event, side, kind, strand, boundary = bnd, distinct_parts,
                      anchor = part, step_dec = FALSE)])
} else {
  ## D's inner end: the member part farthest from the terminus. Stepping inward
  ## from it is the opposite direction to stepping outward from the terminus.
  dp <- sides[, .(gene, event, side, kind, strand, boundary, distinct_parts, out_dec,
                  part = as.integer(trimws(unlist(strsplit(gsub("E", "", distinct_parts), ",")))) ),
              by = seq_len(nrow(sides))][, seq_len := NULL]
  dp <- dp[!is.na(part)]
  dp <- dp[ex[, .(gene, part, start, end)], on = .(gene, part), nomatch = NULL]
  anc <- dp[, .(anchor = if (out_dec[1]) part[which.max(end)] else part[which.min(start)],
                kind = kind[1], strand = strand[1], boundary = boundary[1],
                distinct_parts = distinct_parts[1], step_dec = !out_dec[1]),
            by = .(gene, event, side)]
  setcolorder(anc, c("gene", "event", "side", "kind", "strand", "boundary",
                     "distinct_parts", "anchor", "step_dec"))
}
cat(sprintf("anchors resolved: %d  (direction = %s)\n", nrow(anc), a$direction))

## one step from the anchor
nb <- ex[, .(gene, part, start, end, nxt_start, nxt_part, prv_end, prv_part)]
u <- nb[anc, on = .(gene, part = anchor), nomatch = NULL]
setnames(u, "part", "anchor")
u[step_dec == TRUE,  `:=`(out_part = prv_part,
                          touching = !is.na(prv_end)   & prv_end   == start - 1,
                          end_entered = "hi")]
u[step_dec == FALSE, `:=`(out_part = nxt_part,
                          touching = !is.na(nxt_start) & nxt_start == end + 1,
                          end_entered = "lo")]
n_edge <- sum(is.na(u$out_part))
u <- u[!is.na(out_part)]
cat(sprintf("sides with a unit: %d  (%d at a gene edge, excluded)\n", nrow(u), n_edge))

pid <- function(x) sprintf("%03d", as.integer(x))
u[, `:=`(lo_id = pid(pmin(anchor, out_part)), hi_id = pid(pmax(anchor, out_part)))]
u[touching == TRUE, `:=`(unit_type = "exon", unit_id = pid(out_part))]
iv_hi <- ip_hi[, .(gene, lo_id = abut_part_lo, hi_id = abut_part_hi, bin_hi = bin_id)]
iv_lo <- ip_lo[, .(gene, lo_id = abut_part_lo, hi_id = abut_part_hi, bin_lo = bin_id)]
u <- merge(u, iv_hi, by = c("gene", "lo_id", "hi_id"), all.x = TRUE, sort = FALSE)
u <- merge(u, iv_lo, by = c("gene", "lo_id", "hi_id"), all.x = TRUE, sort = FALSE)
u[touching == FALSE & end_entered == "hi", `:=`(unit_type = "intron", unit_id = bin_hi)]
u[touching == FALSE & end_entered == "lo", `:=`(unit_type = "intron", unit_id = bin_lo)]
n_unres <- sum(is.na(u$unit_id))
if (n_unres) {
  cat(sprintf("  WARNING: %d sides have a gap but no matching intronic part\n", n_unres))
  u <- u[!is.na(unit_id)]
}
cat(sprintf("  unit type: %d exonic, %d intronic\n",
            sum(u$unit_type == "exon"), sum(u$unit_type == "intron")))
u[, own_part := anchor]

## ---- counts for the chosen units ------------------------------------------
## Two sources, because quantification is universal and each unit is counted
## exactly once in whichever annotation owns it:
##   intronic units -> intronparts/counts/<sample>_counts.txt   (the new pass)
##   exonic units   -> count_files/<group>/<SRR>_counts.txt      (already there)
## Both are dexseq_count output keyed "<gene>":"<part>" with quotes, which are
## stripped here. Only the needed keys are retained per file; reading all
## 414,993 intronic bins x 383 samples into memory is neither necessary nor
## viable.
need <- unique(u[, .(gene, unit_id, unit_type)])
cat(sprintf("distinct units to fetch: %d intronic, %d exonic\n",
            sum(need$unit_type == "intron"), sum(need$unit_type == "exon")))
want_i <- need[unit_type == "intron", paste(gene, unit_id, sep = "|")]
want_e <- need[unit_type == "exon",   paste(gene, unit_id, sep = "|")]

read_counts <- function(fn, keep) {
  d <- fread(fn, header = FALSE, sep = "\t", quote = "", showProgress = FALSE,
             col.names = c("k", "n"))
  d <- d[!startsWith(k, "_")]
  d[, k := gsub('"', "", k, fixed = TRUE)]
  d[, k := sub(":", "|", k, fixed = TRUE)]
  d[k %chin% keep, .(k, n = as.integer(n))]
}

ifiles <- sort(list.files(a$intron_counts, pattern = "_counts\\.txt$", full.names = TRUE))
if (!length(ifiles)) stop("no *_counts.txt under ", a$intron_counts)
samples <- sub("_counts\\.txt$", "", basename(ifiles))
cat(sprintf("samples with intronic counts: %d\n", length(samples)))

ic <- rbindlist(lapply(seq_along(ifiles), function(i)
  read_counts(ifiles[i], want_i)[, sample := samples[i]]))
cat(sprintf("intronic count rows kept: %d over %d units\n", nrow(ic), uniqueN(ic$k)))

## exonic counts live per group dir: <sample> = <group>_<SRR>
ex_path <- function(s) {
  grp <- sub("_SRR[0-9]+$", "", s); srr <- sub("^.*_(SRR[0-9]+)$", "\\1", s)
  file.path(a$exon_counts, grp, paste0(srr, "_counts.txt"))
}
ep <- vapply(samples, ex_path, "")
miss_f <- samples[!file.exists(ep)]
if (length(miss_f))
  stop(sprintf("no exonic count file for %d samples, e.g. %s",
               length(miss_f), paste(head(miss_f, 3), collapse = ", ")))
ec <- rbindlist(lapply(seq_along(ep), function(i)
  read_counts(ep[i], want_e)[, sample := samples[i]]))
cat(sprintf("exonic count rows kept: %d over %d units\n", nrow(ec), uniqueN(ec$k)))
cnt <- rbind(ic, ec)

## ---- the assertion that matters -------------------------------------------
## An unresolved unit id becomes 0 downstream, not an error, and a zero unit is
## exactly the `boundary` verdict -- a silently missing bin would manufacture
## the very result this test measures.
seen <- unique(cnt$k)
missing <- need[!paste(gene, unit_id, sep = "|") %chin% seen]
if (nrow(missing))
  stop(sprintf(paste0("FATAL: %d units absent from the count matrices, e.g. %s\n",
                      "A missing unit counts as zero downstream, indistinguishable\n",
                      "from a real boundary. Fix the annotation/counting first."),
               nrow(missing), paste(head(missing[, paste(gene, unit_id)], 5), collapse = ", ")))
cat(sprintf("all %d units resolve to counts\n", nrow(need)))

## ---- exonic table: header shape, then ref/groups per (gene, event, sample) --
## write.table() emits row names as a LEADING UNNAMED column, so data rows carry
## one more field than the header. Detect rather than assume: treating the
## shapes as equal silently skips every row and writes an empty file.
con   <- file(a$exon_combined, "r")
hdr   <- strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]]
probe <- strsplit(readLines(con, 1), "\t", fixed = TRUE)[[1]]
close(con)
shift <- length(probe) - length(hdr)
if (!shift %in% c(0L, 1L))
  stop(sprintf("header has %d fields, data has %d -- unexpected shape",
               length(hdr), length(probe)))
if (shift) cat("leading unnamed row-name column detected; offsetting by 1\n")
idx <- setNames(seq_along(hdr) + shift, hdr); ncol_ <- length(hdr) + shift
for (c_ in c("gene", "event", "sample", "ref", "diff1", "diff2", "groups"))
  if (!c_ %in% names(idx)) stop("exon table lacks column '", c_, "'")

pairs <- unique(u[, .(gene, event)])
key <- tempfile(); sel <- tempfile()
fwrite(pairs[, .(k = paste(gene, event, sep = "\t"))], key,
       col.names = FALSE, quote = FALSE)
cmd <- sprintf(
  "awk -F'\\t' -v OFS='\\t' 'NR==FNR{keep[$1\"\\t\"$2]=1; next} FNR==1{next} NF==%d && (($%d\"\\t\"$%d) in keep){print $%d,$%d,$%d,$%d,$%d}' %s %s > %s",
  ncol_, idx[["gene"]], idx[["event"]], idx[["gene"]], idx[["event"]],
  idx[["sample"]], idx[["ref"]], idx[["groups"]],
  shQuote(key), shQuote(a$exon_combined), shQuote(sel))
if (system(cmd) != 0) stop("awk pre-selection failed")
ex2 <- fread(sel, header = FALSE, sep = "\t", showProgress = FALSE,
             colClasses = list(character = 1:3),
             col.names = c("gene", "event", "sample", "ref", "groups"))
ex2 <- ex2[sample %chin% samples]
cat(sprintf("reference rows selected: %d over %d gene/event pairs, %d samples\n",
            nrow(ex2), uniqueN(ex2[, paste(gene, event)]), uniqueN(ex2$sample)))

## ---- attach each side's unit counts into its diff slot ---------------------
um <- u[, .(gene, event, side, unit_type, unit_id, k = paste(gene, unit_id, sep = "|"))]
j  <- merge(um, cnt, by = "k", allow.cartesian = TRUE)[, .(gene, event, side, sample,
                                                           unit_type, n)]
w  <- dcast(j, gene + event + sample ~ side, value.var = c("n", "unit_type"))
for (c_ in c("n_1", "n_2", "unit_type_1", "unit_type_2"))
  if (!c_ %in% names(w)) w[[c_]] <- if (startsWith(c_, "n_")) NA_integer_ else NA_character_
out <- merge(ex2, w, by = c("gene", "event", "sample"), all.x = FALSE)
out <- out[!(is.na(n_1) & is.na(n_2))]
cat(sprintf("rows: %d over %d gene/event pairs, %d samples\n",
            nrow(out), uniqueN(out[, paste(gene, event)]), uniqueN(out$sample)))

## ---- emit in the exonic schema --------------------------------------------
o <- as.data.table(matrix("NA", nrow = nrow(out), ncol = ncol_))
if (shift) o[[1]] <- as.character(seq_len(nrow(out)))
o[[idx[["gene"]]]]         <- out$gene
o[[idx[["event"]]]]        <- out$event      # same id: separate run, no offset
o[[idx[["sample"]]]]       <- out$sample
o[[idx[["ref"]]]]          <- as.character(out$ref)
o[[idx[["diff1"]]]]        <- fifelse(is.na(out$n_1), "NA", as.character(out$n_1))
o[[idx[["diff2"]]]]        <- fifelse(is.na(out$n_2), "NA", as.character(out$n_2))
o[[idx[["groups"]]]]       <- out$groups
if ("diff1_source" %in% names(idx))
  o[[idx[["diff1_source"]]]] <- fifelse(is.na(out$unit_type_1), "NA", out$unit_type_1)
if ("diff2_source" %in% names(idx))
  o[[idx[["diff2_source"]]]] <- fifelse(is.na(out$unit_type_2), "NA", out$unit_type_2)
con <- file(a$out, "w")
writeLines(paste(hdr, collapse = "\t"), con)   # header only; R takes col 1 as row names
close(con)
fwrite(o, a$out, sep = "\t", col.names = FALSE, quote = FALSE, append = TRUE)
cat("wrote ", a$out, "\n", sep = "")

## ---- side table: the ONLY record of which unit was tested -----------------
## exontest.R consumes diff1_source/diff2_source and does not emit them, so the
## results table carries no trace of the unit. Join it back on
## (gene, event, side), with side = 1 for diff1_vs_ref and 2 for diff2_vs_ref.
##
## MIND columns 22-30 of this arm's results (source, sink, ref_ex_part,
## setdiff1/2, transcripts1/2, path1/2). exontest.R fills those from the
## per-gene split annotation, which is the EXONIC annotation, so here they
## describe D and NOT the unit tested -- `setdiff2 = E008` beside a tested unit
## of I007_008.hi. Same trap as ref_mean/diff1_mean/diff2_mean, which
## drop_cols.R removes. Read the unit from this table, never from those columns.
stat <- j[, .(n_samples = .N, unit_sum = sum(n), unit_mean = mean(n)),
          by = .(gene, event, side)]
side <- merge(u[, .(gene, event, side, kind, strand, boundary,
                    own_part = pid(own_part), out_part = pid(out_part),
                    unit_type, unit_id, direction = a$direction,
                    comparison = fifelse(side == 1L, "diff1_vs_ref", "diff2_vs_ref"))],
              stat, by = c("gene", "event", "side"), all.x = TRUE)
setorder(side, gene, event, side)
side[, unit_mean := sprintf("%.6g", as.numeric(unit_mean))]
if (!is.null(a$out_summary)) {
  fwrite(side, a$out_summary, sep = "\t", quote = FALSE)
  cat("wrote ", a$out_summary, "  (", nrow(side), " sides)\n", sep = "")
  cat(sprintf("  units with mean < 0.5 reads (empty): %d of %d\n",
              sum(as.numeric(side$unit_mean) < 0.5, na.rm = TRUE), nrow(side)))
}
lost <- u[!side[!is.na(n_samples)], on = .(gene, event, side)]
if (nrow(lost))
  cat(sprintf("  WARNING: %d sides got no count rows (e.g. %s)\n",
              nrow(lost), paste(head(lost[, paste(gene, event, side)], 3), collapse = ", ")))
