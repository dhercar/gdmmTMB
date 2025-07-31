#' @title  Ispline Basis with Degree 2 by Default
#'
#' @description
#' A thin wrapper around `splines2::isp()`. Uses `degree = 2`
#' unless otherwise specified. All other arguments are passed through.
#'
#' @inheritParams splines2::isp
#' @param degree Spline degree; defaults to 2 here (instead of 3).
#' @export

isp <- function(x, degree = 2, ...) {
  splines2::isp(x, degree = degree, ...)
}
