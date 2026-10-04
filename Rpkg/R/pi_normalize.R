## Length-normalized pi.
##
## pi as computed by the test is a RAW COUNT ratio, D/(D+S). That is the right
## quantity for the beta-binomial likelihood, which is defined on counts and
## cannot take normalized values. But it is NOT a molecular proportion whenever
## the distinct set and the reference differ in length: a 221 bp distinct set
## against a 67 bp reference inflates pi by the length ratio.
##
## LCP2 event 26 is the worked case -- reported pi_ref 0.548 corresponds to
## about 25% on a per-base basis.
##
## These functions add the per-base reading alongside the raw one. They change
## no call: significance, direction and the reported pi are untouched.


#' Convert a raw-count pi to a per-base pi
#'
#' Closed form, so an existing results table can be annotated without
#' recounting anything:
#'
#'   pi_perbase = pi * len_s / (pi * len_s + (1 - pi) * len_d)
#'
#' Derivation: pi_perbase = (D/len_d) / (D/len_d + S/len_s); multiply through by
#' len_d*len_s and divide by (D+S).
#'
#' @param pi Raw-count pi, in [0, 1].
#' @param len_d,len_s Lengths in bases of the distinct set and the reference.
#' @return Per-base pi, or NA where a length is missing or non-positive.
#' @export
pi_perbase <- function(pi, len_d, len_s) {
  bad <- is.na(pi) | is.na(len_d) | is.na(len_s) | len_d <= 0 | len_s <= 0
  out <- rep(NA_real_, length(pi))
  num <- pi * len_s
  den <- pi * len_s + (1 - pi) * len_d
  ok <- !bad & den > 0
  out[ok] <- num[ok] / den[ok]
  out
}


#' Exonic part lengths for one gene
#'
#' @param gff_path Path to a per-gene \code{<gene>.dexseq.gff}.
#' @return Named numeric vector, part number (as character) to length in bases.
#' @export
part_lengths_from_gff <- function(gff_path) {
  if (!file.exists(gff_path)) return(numeric(0))
  ## quote = "" is REQUIRED: GFF attributes are double-quoted, and read.delim's
  ## default quoting mangles the attribute field so the part number never parses.
  x <- utils::read.delim(gff_path, header = FALSE, comment.char = "#",
                         quote = "", stringsAsFactors = FALSE)
  x <- x[x$V3 == "exonic_part", , drop = FALSE]
  if (!nrow(x)) return(numeric(0))
  n <- sub('.*exonic_part_number "([^"]+)".*', "\\1", x$V9)
  stats::setNames(x$V5 - x$V4 + 1, as.character(as.integer(n)))
}


#' Total length of a comma-separated exonic part set
#'
#' @param part_str e.g. "E019,E020,E021", or NA / "" / "NA" for an empty set.
#' @param lens Named vector from \code{part_lengths_from_gff}.
#' @return Summed length, or NA for an empty or unresolvable set.
#' @export
feature_length <- function(part_str, lens) {
  if (is.na(part_str) || part_str %in% c("", "NA")) return(NA_real_)
  p <- as.integer(sub("^E", "", trimws(strsplit(part_str, ",")[[1]])))
  p <- p[!is.na(p)]
  if (!length(p)) return(NA_real_)
  v <- lens[as.character(p)]
  if (all(is.na(v))) return(NA_real_)
  sum(v, na.rm = TRUE)
}


#' Locate a gene's DEXSeq GFF under either supported layout
#'
#' Two layouts are in use and both must work, because the same package annotates
#' both projects:
#'   flat    \code{<dir>/<gene>.dexseq.gff}            (GrASE_simulation/dexseq.gff)
#'   nested  \code{<dir>/<gene>/<gene>.dexseq.gff}     (DICE grase_results/gene_files)
#'
#' @param gff_dir Directory holding per-gene GFFs in either layout.
#' @param gene Gene id.
#' @return Path to the GFF, or NA_character_ if neither layout resolves.
#' @export
gene_gff_path <- function(gff_dir, gene) {
  flat <- file.path(gff_dir, paste0(gene, ".dexseq.gff"))
  if (file.exists(flat)) return(flat)
  nested <- file.path(gff_dir, gene, paste0(gene, ".dexseq.gff"))
  if (file.exists(nested)) return(nested)
  NA_character_
}


#' Annotate a results table with per-base pi
#'
#' Adds \code{len_D}, \code{len_S}, \code{pi_ref_perbase},
#' \code{pi_trt_perbase} and \code{delta_pi_perbase}. Rows whose distinct set
#' is a JUNCTION rather than exonic parts are left NA: a junction is a point
#' feature with no length, so a per-base reading is undefined for it.
#'
#' @param res Results data frame with pi_ref, pi_trt, ref_ex_part, setdiff1,
#'   setdiff2, comparison, gene.
#' @param gff_dir Directory of per-gene \code{<gene>.dexseq.gff} files.
#' @return \code{res} with the new columns appended.
#' @export
annotate_pi_perbase <- function(res, gff_dir) {
  side <- ifelse(grepl("^diff1", res$comparison), 1L, 2L)
  sd   <- ifelse(side == 1L, res$setdiff1, res$setdiff2)
  res$len_D <- NA_real_
  res$len_S <- NA_real_
  cache <- new.env(parent = emptyenv())
  ## A gene with no GFF yields the same NA as a junction side, and the two mean
  ## different things -- one is undefined by nature, the other is missing input.
  ## Track it so the gap is reported rather than absorbed into the NA count.
  unresolved <- character(0)
  for (g in unique(res$gene)) {
    lens <- cache[[g]]
    if (is.null(lens)) {
      p <- gene_gff_path(gff_dir, g)
      if (is.na(p)) unresolved <- c(unresolved, g)
      lens <- if (is.na(p)) numeric(0) else part_lengths_from_gff(p)
      cache[[g]] <- lens
    }
    if (!length(lens)) next
    i <- which(res$gene == g)
    res$len_D[i] <- vapply(sd[i], feature_length, numeric(1), lens = lens)
    res$len_S[i] <- vapply(res$ref_ex_part[i], feature_length, numeric(1),
                           lens = lens)
  }
  if (length(unresolved))
    warning(sprintf(paste0("no DEXSeq GFF found for %d of %d genes under '%s' ",
                           "(%d rows left NA for want of lengths, not because ",
                           "the side is a junction); e.g. %s"),
                    length(unresolved), length(unique(res$gene)), gff_dir,
                    sum(res$gene %in% unresolved),
                    paste(utils::head(unresolved, 3), collapse = ", ")),
            call. = FALSE)
  res$pi_ref_perbase   <- pi_perbase(res$pi_ref, res$len_D, res$len_S)
  res$pi_trt_perbase   <- pi_perbase(res$pi_trt, res$len_D, res$len_S)
  res$delta_pi_perbase <- res$pi_trt_perbase - res$pi_ref_perbase
  res
}
