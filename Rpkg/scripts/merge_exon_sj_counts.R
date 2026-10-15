library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("-e", "--exon_counts"), type="character",
              help="directory of per-gene exon count files (from exoncnt.R)",
              metavar="character"),
  make_option(c("-s", "--sj_counts"), type="character",
              help="directory of per-gene SJ count files",
              metavar="character"),
  make_option(c("-o", "--output"), type="character",
              help="output directory for merged count files",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$exon_counts) || is.null(opt$sj_counts) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("--exon_counts, --sj_counts, and --output are required")
}

exon_dir <- opt$exon_counts
sj_dir   <- opt$sj_counts
out_dir  <- opt$output

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

## Input is the `ex` unit written by exoncnt.R: <gene>.bipartition.ex.<class>.cnt.txt
## Output is the `exsj` unit -- same rows, but a side with no exclusive exonic
## part takes its exclusive junction instead. The UNIT LABEL IN THE NAME is the
## only thing distinguishing the two on disk, and the layouts differ (12 columns
## vs 17), so a mislabelled file silently misaligns every downstream parser.
exon_files <- list.files(exon_dir, pattern="\\.ex\\.[^.]+\\.cnt\\.txt$", full.names=TRUE)
if (length(exon_files) == 0L)
  stop("no .ex.<class>.cnt.txt files found in: ", exon_dir,
       "\n  (pre-rename directories hold .exoncnt.txt and are no longer read; ",
       "see rename_units.py)")

cat(sprintf("Found %d exon count files\n", length(exon_files)))

n_genes_merged <- 0L
n_subst1_total <- 0L
n_subst2_total <- 0L

for (f in exon_files) {
  bname   <- basename(f)
  ## the junction counts are unit-agnostic (they ARE the junction source), so
  ## their name carries no unit: <gene>.bipartition.sjcnt.txt
  sj_file <- file.path(sj_dir, sub("\\.ex\\.[^.]+\\.cnt\\.txt$", ".sjcnt.txt", bname))
  ## ex -> exsj on the way out; everything else about the name is preserved
  oname   <- sub("\\.ex\\.", ".exsj.", bname, fixed = FALSE)
  if (identical(oname, bname))
    stop("input name does not carry the `ex` unit label, refusing to write an ",
         "unlabelled exsj file: ", bname)

  ex <- read.table(f, header=TRUE, sep="\t", row.names=NULL,
                   stringsAsFactors=FALSE, na.strings=c("NA",""))
  # original exon count files have a row-name column; drop it
  if ("row.names" %in% names(ex)) ex <- ex[, names(ex) != "row.names"]

  # initialize SJ annotation columns as NA; will be filled if SJ file exists.
  # transcripts1/2 and path1/2 are intentionally excluded: they contain spaces
  # in the transcript ID lists which break exontest.R's whitespace-delimited
  # read.table, and they are not used by the statistical test.
  ex$intron_distinct1 <- NA_character_
  ex$intron_distinct2 <- NA_character_
  ex$intron_shared    <- NA_character_
  ex$diff1_source     <- "exon"
  ex$diff2_source     <- "exon"

  if (!file.exists(sj_file)) {
    write.table(ex, file.path(out_dir, oname), sep="\t", quote=FALSE, row.names=FALSE)
    next
  }

  sj <- read.table(sj_file, header=TRUE, sep="\t", row.names=NULL,
                   stringsAsFactors=FALSE, na.strings=c("NA",""))

  # if SJ samples have an alignment-pass tag (e.g. _pass2_), keep only the
  # pass2 rows and strip the tag so sample names match the exon count file.
  if (any(grepl("_pass[0-9]+_", sj$sample))) {
    sj <- sj[grepl("_pass2_", sj$sample), ]
    sj$sample <- sub("_pass2_", "_", sj$sample)
  }

  # one row per event: metadata that does not vary by sample.
  # bipartition_sjcnt.R now writes it ONCE to a .sjmeta.txt sidecar instead of
  # repeating it on every sample row (with 1,255 samples that was 99% of the
  # file). Fall back to the in-file columns for sjcnt output made before that
  # change, so old and new trees both work.
  meta_file <- sub("\\.sjcnt\\.txt$", ".sjmeta.txt", sj_file)
  if (file.exists(meta_file)) {
    sj_meta <- read.table(meta_file, header=TRUE, sep="\t", row.names=NULL,
                          stringsAsFactors=FALSE, na.strings=c("NA",""))
    sj_meta <- sj_meta %>%
      distinct(gene, event, intron_distinct1, intron_distinct2, intron_shared)
  } else {
    sj_meta <- sj %>%
      distinct(gene, event, intron_distinct1, intron_distinct2, intron_shared)
  }

  # drop placeholder SJ columns, then join from SJ file
  ex <- ex %>%
    select(-intron_distinct1, -intron_distinct2, -intron_shared) %>%
    left_join(sj_meta, by=c("gene","event"))

  # determine which events need substitution (per event, not per sample)
  event_flags <- ex %>%
    distinct(gene, event, setdiff1, setdiff2, intron_distinct1, intron_distinct2)

  subst1_events <- event_flags %>%
    filter(is.na(setdiff1), !is.na(intron_distinct1)) %>%
    pull(event)

  subst2_events <- event_flags %>%
    filter(is.na(setdiff2), !is.na(intron_distinct2)) %>%
    pull(event)

  if (length(subst1_events) > 0L) {
    sj_d1 <- sj %>%
      filter(event %in% subst1_events) %>%
      select(gene, event, sample, sj_diff1=diff1)
    ex <- ex %>%
      left_join(sj_d1, by=c("gene","event","sample")) %>%
      mutate(
        diff1_source = if_else(event %in% subst1_events, "sj", diff1_source),
        diff1        = if_else(event %in% subst1_events, sj_diff1, diff1)
      ) %>%
      select(-sj_diff1)
  }

  if (length(subst2_events) > 0L) {
    sj_d2 <- sj %>%
      filter(event %in% subst2_events) %>%
      select(gene, event, sample, sj_diff2=diff2)
    ex <- ex %>%
      left_join(sj_d2, by=c("gene","event","sample")) %>%
      mutate(
        diff2_source = if_else(event %in% subst2_events, "sj", diff2_source),
        diff2        = if_else(event %in% subst2_events, sj_diff2, diff2)
      ) %>%
      select(-sj_diff2)
  }

  write.table(ex, file.path(out_dir, oname), sep="\t", quote=FALSE, row.names=FALSE)

  n_genes_merged  <- n_genes_merged  + 1L
  n_subst1_total  <- n_subst1_total  + length(subst1_events)
  n_subst2_total  <- n_subst2_total  + length(subst2_events)
}

cat(sprintf("Genes with SJ data: %d\n",       n_genes_merged))
cat(sprintf("Events with diff1 from SJ: %d\n", n_subst1_total))
cat(sprintf("Events with diff2 from SJ: %d\n", n_subst2_total))
