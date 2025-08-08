#' Nice titles for summary tables
#'
#' @param text Title to display
#' @param width Size of the console. Defaults to current width.
#' @param symb Symbol used above and below the title
#'
#' @return A `character` string (invisibly returned)
#'
#' @examples print_title('hello world', '=')

print_title <- function(text, symb = '-', width = getOption("width")) {
  # Full-width line of dashes
  separator_line <- strrep(symb, width)
  # Centered text
  padding <- floor((width - nchar(text)) / 2)
  if (padding < 0) padding <- 0
  centered_text <- paste0(strrep(" ", padding), text)

  # Print all
  cat(separator_line, "\n")
  cat(centered_text, "\n")
  cat(separator_line, "\n\n")
}


#' Centered title with decorated side
#'
#' @param text Title to display
#' @param width Size of the console. Defaults to current width.
#' @param symb Symbol used left and right of the title
#'
#' @return A `character` string (invisibly returned)
#'
#' @examples print_title2(' hello world ', ':')

print_title2 <- function(text, symb = '-', width = getOption("width")) {
  total_padding <- width - nchar(text)
  if (total_padding < 0) total_padding <- 0

  left_padding <- floor(total_padding / 2)
  right_padding <- ceiling(total_padding / 2)

  decorated_text <- paste0(strrep(symb, left_padding), text, strrep(symb, right_padding))
  cat(decorated_text, "\n\n")
}


#' Bray-Curtis dissimilarity components
#' @description
#' Calculate the numerator and denominator components of the Bray-Curtis dissimilarity index.
#' The Bray-Curtis dissimilarity is calculated as \eqn{\frac{\sum_{i} |x_i - y_i|}{\sum_{i} (x_i + y_i)}}
#'
#' @param x numeric vector of the first sample to be compared
#' @param y numeric vector of the second sample to be compared
#'
#' @return A named numeric vector with two elements:
#' \itemize{
#'   \item num - numerator: sum of absolute differences between x and y
#'   \item den - denominator: sum of x and y values
#' }
#' @examples bray1(c(0,1,1), c(1,1,0))

bray1 <- function(x,y){
  num = sum(abs(x - y))
  den = sum(x,y)
  return(c(num = num, den = den))
}



#' Inverse logit
#'
#' @param x 	A numeric value or vector
#' @description
#' Performs inverse of logit: exp(x) / (1 + exp(x))
#'
#' @return A (transformed) numeric value or vector
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}


#' Scale a numeric vector to [0, 1]
#'
#' Takes a numeric vector and linearly rescales it so that its minimum becomes 0 and its maximum becomes 1.
#'
#' @param x A numeric vector.
#' @return A numeric vector of the same length as \code{x}, with values in the interval \[0, 1\].
#' @examples
#' \dontrun{
#'   scale01(c(10, 20, 30))
#'   #> 0.0, 0.5, 1.0
#' }
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}


#' Scale a [0, 1] vector into an arbitrary subinterval
#'
#' Given a numeric vector \code{x}, this function re-scales \code{x} into the interval \[\code{min}, \code{max}\].
#'
#' @param x A numeric vector with values expected in \[0, 1\].
#' @param new Numeric scalar: the upper nad lower bound of the target interval. Default is \[\code{0.01}, \code{0.99}\].
#' @return A numeric vector of the same length as \code{x}, with values in the interval \[\code{min}, \code{max}\].
#' @examples
#' \dontrun{
#'   # first scale to [0,1], then to [0.1, 0.9]
#'   v <- c(5, 15, 25)
#'   v01 <- scale01(v)
#'   scale_dist(v01, min = 0.1, max = 0.9)
#' }
scale_dist <- function(x, new = c(0.01, 0.99)) {
  min(new) + x * (max(new) - min(new))
}

