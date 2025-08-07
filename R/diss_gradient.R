#' @title Predicted Effects of Environmental Gradients on Dissimilarity
#'
#' @description
#' Computes the predicted effect of one or more covariates on compositional dissimilarity
#' from a fitted `gdmm` or `bbgdmm` model. Optionally returns uncertainty bands (credible or confidence intervals)
#' using either bootstrapped samples or parametric bootstrapping.
#'
#' @param m A fitted model object of class `gdmm` or `bbgdmm`, representing a generalized dissimilarity (mixed) model.
#' @param var Character vector of variable names to evaluate. Use `"all"` (default) to compute for all predictors in the model.
#' @param n Integer. Number of points to evaluate along each gradient. Defaults to 100.
#' @param CI Logical. Whether to compute confidence intervals. Default is `TRUE`.
#' @param n_sim Number of draws used to obtain CI. If `NULL`, uses `m$n_boot` for `bbgdmm` or 1000 for `gdmm`.
#' @param CI_quant Numeric vector specifying the quantile coverage (e.g., `0.95` for 95% CI). Default is `0.95`.
#'
#' @details
#' The function generates predicted dissimilarities along a gradient of each specified predictor,
#' while holding all other predictors constant at their minimum observed values. If `CI = TRUE`,
#' it uses bootstrap samples (for `bbgdmm`) or simulates from the joint covariance matrix (for `gdmm`)
#' to construct confidence or credible intervals.
#'
#' This approach is useful for visualizing how compositional changes to individual
#' environmental or spatial gradients, as part of a generalized dissimilarity modeling framework (Ferrier et al., 2007).
#'
#' @return A named list of data frames, one for each variable in `var`. Each data frame contains:
#' \itemize{
#'   \item \code{var}: The variable name.
#'   \item \code{f_x}: The predicted dissimilarity response.
#'   \item \code{x}: The values of the focal gradient.
#'   \item \code{CI lower / CI upper}: (Optional) Lower and upper bounds of the confidence interval.
#' }
#'
#' @references
#' Ferrier, S., Manion, G., Elith, J., & Richardson, K. (2007).
#' Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment.
#' \emph{Diversity and Distributions}, 13(3), 252–264. \doi{10.1111/j.1472-4642.2007.00341.x}
#'
#' @importFrom hardhat forge
#' @importFrom TMB sdreport
#'
#' @export
diss_gradient <- function(m,
                          var = 'all',
                          n = 100,
                          CI = TRUE,
                          n_sim = NULL,
                          CI_quant = c(0.95)){

  all_vars <- all.vars(m$diss_formula)

  if (var == 'all') {
    var <- all.vars(m$diss_formula)
  }

  stopifnot(all(var %in% all.vars(m$diss_formula)))

  # Mean beta
  if (class(m) == 'gdmm') {
    sdr <- TMB::sdreport(m$obj)
    mean_beta <- sdr$value[names(sdr$value) == 'e_beta']
  } else if (class(m) == 'bbgdmm') {
    mean_beta <- apply(m$boot_samples[,colnames(m$boot_samples) %in% c('e_beta'), drop = F], 2, mean)
  }

  # SIM
  if (CI) {
    # Param combinations
    if (class(m) == 'bbgdmm') {
      if (is.null(n_sim)) n_sim = m$n_boot
      sims <- m$boot_samples[,colnames(m$boot_samples) %in% c('e_beta'), drop = F]

    } else if (class(m) == 'gdmm')  {
      if (is.null(n_sim)) n_sim = 1000
      sims <- MASS::mvrnorm(n_sim, sdr$value, sdr$cov)
    }
    sims <- sims[,colnames(sims) == 'e_beta', drop = F]
  }

  # prep x data
  x_mean <- apply(m$X[,all_vars, drop = F], 2, min)
  x_min <- matrix(rep(x_mean, each = n), nrow = n, ncol =  length(x_mean), byrow = F)
  colnames(x_min) <- all_vars

  # mean f(x)
  out_list <- list()

  for (v in var) {
    x_new <- x_min
    x_new[,v] <- seq(min(m$X[,v]), max(m$X[,v]), length.out = n)
    x_new <- apply(x_new, 2, as.numeric)
    suppressWarnings({
      f_x <- as.matrix(hardhat::forge(x_new, m$form_X$blueprint)$predictors) %*% mean_beta
    })

    out_i <- data.frame(var = v,
                        f_x = f_x,
                        x = x_new[,v])

    if (CI) {
      quantiles <- sort(c((1 - CI_quant)*0.5, CI_quant + (1 - CI_quant)*0.5))
      samples <- do.call(cbind, lapply(1:n_sim, function(i){
        beta_i <- sims[i,,drop = F]
        suppressWarnings({
          as.matrix(hardhat::forge(x_new, m$form_X$blueprint)$predictors) %*% c(beta_i)
        })
      }))

      CI_i <- t(apply(samples, 1, function(x){quantile(x, probs = quantiles)}))
      colnames(CI_i) <- c(paste0('CI ', colnames(CI_i)))

      out_i <- cbind(out_i, CI_i)
    }
    out_list[[v]] <- out_i
  }

  return(out_list)
}






