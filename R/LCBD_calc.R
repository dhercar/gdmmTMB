#' @title Compute Sum of Squares or LCBD from a Symmetric Dissimilarity Matrix
#'
#' @description
#' Computes the per-site sum of squares (SS) from a symmetric dissimilarity matrix using
#' double-centering. Optionally returns LCBD values (relative contributions to total SS).
#'
#' @param x A numeric symmetric matrix (e.g., a squared Euclidean dissimilarity matrix).
#' @param LCBD Logical; if `TRUE`, returns normalized values summing to 1 (analogous to LCBD).
#'
#' @details
#' The function double-centers the dissimilarity matrix \code{x} using the Gower method
#' (i.e., \code{A = -0.5 * H %*% x %*% H}, where \code{H = I - 1/n}), and extracts the diagonal
#' of the resulting matrix as per-site sum of squares. If \code{LCBD = TRUE}, the diagonal is
#' normalized to sum to 1, representing each site's relative contribution to total variation (Legendre & De Cáceres, 2013).
#'
#' @return A numeric vector of length \code{nrow(x)}. If \code{LCBD = TRUE}, the values sum to 1.
#'
#' @references
#' Legendre, P., & De Cáceres, M. (2013). Beta diversity as the variance of community data:
#' dissimilarity coefficients and partitioning. \emph{Ecology Letters}, 16(8), 951–963.
#' \doi{10.1111/ele.12141}
#'
#' @examples
#' # Create a simple symmetric dissimilarity matrix
#' d <- as.matrix(dist(matrix(rnorm(20), nrow = 5)))
#' SS_calc(d)         # raw per-site SS
#' SS_calc(d, TRUE)   # LCBD-like values
#'
#' @seealso \code{\link[adespatial]{LCBD.comp}}
#'
#' @export
SS_calc <- function(x, LCBD = FALSE){
  stopifnot(is.matrix(x), isSymmetric(x))

  n <- nrow(x)
  H <- diag(n) - matrix(1, n, n)/n
  A <- -0.5 * H %*% x %*% H
  SS <- diag(A)

  if (LCBD) {
    return(SS/sum(SS))
  } else {
    return(SS)
  }
}
