# Common utility functions

#' Count items in comma-separated string
#' @param x Character string. A comma-separated list of items (e.g., \code{"E001,E002,E003"}).
#'   An empty string returns \code{0}.
#' @export
#' @examples
#' count_items("E001,E002,E003")
#' count_items("")
#' count_items("E001")
count_items <- function(x) {
  if (x == "") return(0)
  length(stringr::str_split(x, ",")[[1]])
}

# Parse comma-separated exonic part string into integer exon numbers.
# Returns integer(0) for empty/NA input.
exon_nums_from_str <- function(x) {
  if (is.na(x) || trimws(x) == "" || x == "NA") return(integer(0))
  parts <- trimws(unlist(strsplit(x, ",")))
  as.integer(sub("^E0*", "", parts))
}

# Minimum distance between ref exon numbers and the union of all setdiff exon numbers.
# setdiff_strs: character vector of one or more comma-separated exonic part strings.
# Returns Inf when either side is empty.
ref_proximity_to_setdiffs <- function(ref_str, ...) {
  r <- exon_nums_from_str(ref_str)
  s <- unique(unlist(lapply(c(...), exon_nums_from_str)))
  if (length(r) == 0 || length(s) == 0) return(Inf)
  min(outer(r, s, function(a, b) abs(a - b)))
}

#' Intersect all elements in a list efficiently
#' @param lst A list of vectors. All elements are intersected
#'   together.
#' @export
intersect_stream <- function(lst) {
  if (length(lst) == 0L) return(vector(mode = "integer", length = 0))
  # sort by length so we intersect the smallest sets first
  ord <- order(vapply(lst, length, integer(1)))
  curr <- unique(lst[[ord[1]]])
  for (i in ord[-1]) {
    # filter in place
    curr <- curr[curr %in% lst[[i]]]
    if (length(curr) == 0L) break
  }
  curr
}
