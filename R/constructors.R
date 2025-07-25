#' @noRd
new_gdmm <- function(fit) {
  class(fit) <- "gdmm"
  fit
}

#' @noRd
new_bbgdmm <- function(fit) {
  class(fit) <- "bbgdmm"
  fit
}
