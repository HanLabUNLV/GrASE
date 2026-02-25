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

#' Intersect all elements in a list efficiently
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
