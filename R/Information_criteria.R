#' Calculate AICc
#'
#' Calculates the Akaike Information Criterion corrected for small sample sizes (AICc).
#'
#' @param log_likelihood The log-likelihood of the model
#' @param n The number of observations
#' @param k The number of estimated parameters

#' @examples
#' # Example usage
#' log_likelihood <- -150
#' n <- 50
#' k <- 5
#' AICc(log_likelihood, n, k)


AICc <- function(log_likelihood, # model log-likelihood,
                 n,              # n samples
                 k) {            # k parameters estimated

  AIC <- -2 * log_likelihood + 2 * k
  correct <- (2 * k * (k + 1)) / (n - k - 1)

  return(list(AIC = AIC, AICc = AIC + correct))
}


#' Bayesian Information Criterion (BIC)
#'
#' Calculates the Bayesian Information Criterion (BIC).
#'
#' @param log_likelihood The log-likelihood of the model
#' @param n The number of observations
#' @param k The number of estimated parameters
#'
#' @examples
#' # Example usage of BICcalc function
#' log_likelihood <- -150
#' n <- 100
#' k <- 5
#' BICcalc(log_likelihood, n, k)


BIC2 <- function(log_likelihood, # model log-likelihood,
                 n,              # n samples
                 k) {             # k parameters estimated

  BIC <- -2 * log_likelihood + k * log(n)
  return(BIC)
}
