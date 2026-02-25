.onLoad <- function(libname, pkgname) {
  op <- options()
  op.grase <- list(
    grase.debug = FALSE
  )
  toset <- !(names(op.grase) %in% names(op))
  if (any(toset)) options(op.grase[toset])
}

log_debug <- function(msg) {
  if (getOption("grase.debug", default = FALSE)) {
    message("[DEBUG] ", msg)
  }
}
