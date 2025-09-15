#' Predict method for gdmm and bbgdmm objects
#'
#' @description
#' Generate predictions from fitted generalized dissimilarity mixed models.
#' This method can predict dissimilarity or uniqueness for new data,
#' with optional confidence intervals. Works with both regular gdmm and bootstrap bbgdmm objects.
#'
#' @param object A fitted model object of class 'gdmm' or 'bbgdmm'
#' @param new_X A data frame or matrix of new predictor variables for the dissimilarity component.
#'   If NULL (default), uses the original data from the fitted model.
#' @param new_X_pair A data frame or matrix of new pairwise predictor variables for the dissimilarity component, usually distances.
#' @param new_W A data frame or matrix of new predictor variables for the uniqueness component.
#' @param new_re A data frame of new random effect variables. If NULL, uses random effects in the fitted model
#' @param new_D Matrix of pairwise combination indices. Required when new_X_pair is provided.
#' @param re_sd A character indicating which random effects should be treated as having
#'   standard deviation (used for simulation).
#' @param component Character string specifying what to predict. Either "dissimilarity" (default)
#'   for pairwise dissimilarities or "uniqueness" for site-level uniqueness values.
#' @param scale_uniq Logical indicating whether uniqueness values should be transformed into proportions as in LCBD. Default is TRUE.
#' @param type Character string specifying the type of prediction. Either "response" (default)
#'   for predictions on the response scale, or "link" for predictions on the linear predictor scale.
#' @param CI Logical indicating whether to compute confidence intervals. Default is TRUE.
#' @param CI_quant Numeric vector of confidence levels to compute. Default is \code{c(0.95, 0.5)}
#'   for 95% and 50% confidence intervals.
#' @param n_sim Integer specifying the number of simulations for confidence intervals.
#'   If NULL, uses 1000 for 'gdmm' objects and the number of bootstrap samples for 'bbgdmm' objects.
#' @param n_cores Integer specifying the number of cores to use for parallel computation.
#'   Currently not implemented. Default is 1.
#' @param ... Additional parameters (not used).
#' @method predict gdmm
#' @return
#' If \code{CI = FALSE}, returns a numeric vector of predictions.
#' If \code{CI = TRUE}, returns a matrix with columns for the mean prediction and
#' confidence interval bounds, named according to the confidence levels specified in \code{CI_quant}.
#'
#' @details
#' This function generates predictions from fitted gdmm models by:
#' \itemize{
#'   \item Computing expected values using model coefficients
#'   \item Optionally generating confidence intervals through bootstrapping
#' }
#'
#' For 'gdmm' objects, confidence intervals are computed using parametric bootstrapping.
#'
#' When \code{component = "uniqueness"}, the function computes site-level uniqueness
#' values based on the full dissimilarity matrix. When \code{component = "dissimilarity"},
#' it returns pairwise dissimilarity predictions.
#' If site-level random effects are included via \code{uniq_formula} when fitting the model (e.g.,\code{uniq_formula = (1|site)}),
#' confidence intervals can incorporate site-level uncertainty via \code{re_sd = c('site')}.
#'
#' @importFrom hardhat forge
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm quantile
#' @importFrom utils combn
#' @export
#'
predict.gdmm <- function(object,
                         new_X = NULL,
                         new_W = NULL,
                         new_X_pair = NULL,
                         new_re = NULL,
                         new_D = NULL,
                         re_sd = logical(0),
                         component = 'dissimilarity',
                         scale_uniq = TRUE,
                         type = 'response',
                         CI = TRUE,
                         CI_quant = c(0.95, 0.5),
                         n_sim = NULL,
                         n_cores = 1,
                         ...) {

  # conditional predictions if no data is provided
  if (is.null(new_X) & is.null(new_W) & is.null(new_re) & is.null(new_X_pair)) {
    new_X <- object$X
    new_X_pair <- object$X_pair
    new_D <- object$D
    new_W <- object$X
    new_re <- object$X[,object$re_vars, drop = FALSE]
  }

  if (is.null(new_W)) {
    new_W <- new_X
    message("'new_W' not provided: using new_W = new_X")
  }

  n <- max(nrow(new_X), nrow(new_W), nrow(new_re))

  # sample combinations
  if (is.null(new_D)) {
    D <- t(utils::combn(n, 2))
  } else {
    if (!is.matrix(new_D) && !is.data.frame(new_D)) {
      stop("'new_D' must be a matrix or data frame with two columns")
    }
    if (ncol(new_D) != 2) {
      stop("'new_D' must have exactly two columns")
    }

    D <- as.matrix(new_D)

    if (!all(new_D %in% 1:nrow(new_X))) {
      stop("'new_D' should only contain integers between 1 and nrow(X_new)")
    }
  }


  # ----- EXPECTED VALUE -----
  if (inherits(object, 'gdmm')) {
    beta <- object$obj$report()$e_beta
    beta_p <- object$obj$report()$e_beta_p
    lambda <- object$obj$report()$lambda
    intercept <- object$obj$report()$intercept
    u <- object$obj$report()$u
    sigma <- object$obj$report()$sigma_re[object$re_vars %in% re_sd]
  } else if (inherits(object, 'bbgdmm')) {
    mean_par <- colMeans(object$boot_samples[,colnames(object$boot_samples) %in% c('intercept', 'e_beta', 'lambda', 'u', 'sigma_re')])
    beta <- mean_par[names(mean_par) == 'e_beta']
    beta_p <-mean_par[names(mean_par) == 'e_beta_p']
    lambda <- mean_par[names(mean_par) == 'lambda']
    intercept <- mean_par[names(mean_par) == 'intercept']
    sigma <- mean_par[names(mean_par) == 'sigma_re']
    u <- mean_par[names(mean_par) == 'u']
  }

  # Mean value
  out <- coef_to_pred(
    object = object,
    intercept = intercept,
    beta = beta,
    beta_p = beta_p,
    lambda = lambda,
    sigma = logical(0),
    u = u,
    D = D,
    n = n,
    new_re = new_re,
    new_X = new_X,
    new_X_pair = new_X_pair,
    new_W = new_W,
    type = type,
    scale_uniq = scale_uniq,
    component = component
  )

  if (CI == TRUE) {
    # Param combinations
    if (inherits(obj, 'bbgdmm')) {
      if (is.null(n_sim)) n_sim = object$n_boot
      sims <- object$boot_samples[,colnames(object$boot_samples) %in% c('e_beta','e_beta_p', 'lambda', 'intercept')]

    } else if (inherits(obj, 'gdmm'))  {
      if (is.null(n_sim)) n_sim = 1000
      sdr <- TMB::sdreport(object$obj)
      sims <- MASS::mvrnorm(n_sim, sdr$value,  sdr$cov)
    }

    pred_list <- list()
    for (i in 1:n_sim) {
      beta_i <- sims[i, colnames(sims) == 'e_beta']
      beta_p_i <- sims[i, colnames(sims) == 'e_beta_p']
      lambda_i <- sims[i, colnames(sims) == 'lambda']
      intercept_i <- sims[i, colnames(sims) == 'intercept']
      u_i <- u #sims[i, colnames(sims) == 'u']

      pred_list[[length(pred_list)+1]] <- coef_to_pred(object = object,
                                intercept = intercept_i,
                                beta = beta_i,
                                beta_p = beta_p_i,
                                lambda = lambda_i,
                                sigma = sigma,
                                u = u_i,
                                n = n,
                                new_re = new_re,
                                new_X = new_X,
                                new_X_pair = new_X_pair,
                                D = D,
                                new_W = new_W,
                                type = type,
                                scale_uniq = scale_uniq,
                                component = component)
    }
    quantiles <- sort(c((1-CI_quant)*0.5, CI_quant + (1-CI_quant)*0.5))
    CI <- t(apply(do.call(cbind, pred_list), 1, function(x) stats::quantile(x, quantiles)))
    out <- cbind(mean = out, CI = CI)
    colnames(out) <- c('mean', paste0('CI ', colnames(CI)))
  }
  return(out)
}

#' Predict method for gdmm and bbgdmm objects
#'
#' @describeIn predict.gdmm Method for bootstrap gdmm objects
#' @export
#' @method predict bbgdmm
predict.bbgdmm <- predict.gdmm


#' Convert model coefficients to predictions
#'
#' Internal function that transforms model parameters and new data into
#' dissimilarity or uniqueness predictions.
#'
#' @param object Model object containing fitted model components
#' @param intercept Numeric, model intercept
#' @param beta Numeric vector, coefficients for dissimilarity component
#' @param beta_p Numeric vector, coefficients for dissimilarity component (pairwise)
#' @param lambda Numeric vector, coefficients for uniqueness component
#' @param u Numeric vector, random effect values
#' @param new_W Matrix/data.frame, new data for uniqueness component
#' @param new_X Matrix/data.frame, new data for dissimilarity component
#' @param new_X_pair Matrix/data.frame, new data for dissimilarity component (pairwise)
#' @param new_re Matrix/data.frame, new random effect data
#' @param D Matrix, pairwise combination indices
#' @param n Integer, number of observations
#' @param component Character, "dissimilarity" or "uniqueness"
#' @param scale_uniq Logical, whether to scale uniqueness values
#' @param type Character, "response" or "link" scale
#' @param sigma Numeric, random effect standard deviations
#'
#' @return Numeric vector of predictions (dissimilarities or uniqueness values)
#'
#' @keywords internal
coef_to_pred <- function(object, intercept, beta, beta_p, lambda, u,
                         new_W, new_X, new_X_pair, new_re, D, n, component, scale_uniq, type, sigma) {

  if ((length(lambda) > 0) & !is.null(new_W)) {
    form_W_new <- as.matrix(hardhat::forge(new_W, object$form_W$blueprint)$predictors)
    uniq_comp_pair <- (form_W_new[D[,1],,drop = F] + form_W_new[D[,2],,drop = F]) %*% lambda
    #uniq_comp <- form_W_new %*% lambda
  } else {
    #uniq_comp = 0
    uniq_comp_pair = 0
  }


  if ((length(beta) > 0) & !is.null(new_X)) {
    form_X_new <- as.matrix(hardhat::forge(new_X, object$form_X$blueprint)$predictors)
    diss_comp = abs(form_X_new[D[,1],, drop = F] - form_X_new[D[,2],,drop = F]) %*% c(beta)
  } else {
    diss_comp = 0
  }

  if ((length(beta_p) > 0) & !is.null(new_X_pair)) {
    form_X_p_new <- as.matrix(hardhat::forge(new_X_pair, object$form_X_pair$blueprint)$predictors)
    diss_p_comp = form_X_p_new %*% beta_p
  } else {
    diss_p_comp = 0
  }


  # Random effects matrix
  if (!is.null(new_re) & (length(new_re) > 0)) {
    if (!(all(colnames(new_re) %in% object$re_vars))) stop('random variable(s) "', paste(colnames(new_re)[which(!(colnames(new_re) %in% object$re_vars))], sep = ', '), '" in `new_re` not found in model object')
    new_re_vars <- colnames(new_re)
    new_Z <- sparse.model.matrix(eval(parse(text = paste0('~ 0 + ', paste0(new_re_vars, collapse = ' + ')))),
                                 data = new_re)

    if (!(all(colnames(new_Z) %in% colnames(object$Z_design)))) stop('new levels for random variables are not allowed\n','   new levels found:', paste0(colnames(new_Z)[which(!(colnames(new_Z) %in% colnames(object$Z_design)))], sep = ', '))
    u <- u[match(colnames(new_Z), colnames(object$Z_design))]
    re_comp <- new_Z %*% (u)
    re_comp_pair <- re_comp[D[,1]] + re_comp[D[,2]]

  } else if (length(sigma) > 0 && !is.null(sigma)) {
    re_comp_pair <- rowSums(do.call(cbind, lapply(rep(sigma,2), function(x) stats::rnorm(nrow(D), 0, x))))
  } else {
    re_comp_pair <- 0
  }

  pred_diss <- intercept + diss_comp + diss_p_comp + uniq_comp_pair + re_comp_pair


  if (type == 'response') {
    pred_diss <- switch(object$link,
                        'identity' = pred_diss,
                        'logit' = inv_logit(pred_diss),
                        'neg_exp' = 1 - exp(-pred_diss),
                        'neg_gaus' = 1 - exp(-(pred_diss^2)))
  }

  if (component == 'uniqueness') {
    diss_matrix <- matrix(ncol = n, nrow = n)
    diss_matrix[cbind(D[,1], D[,2])] <- pred_diss
    diss_matrix[cbind(D[,2], D[,1])] <- pred_diss
    diag(diss_matrix) <- 0
    pred_uniq <- SS_calc(diss_matrix, LCBD = scale_uniq)
    names(pred_uniq) <- NULL
    return(pred_uniq)
  } else {
    names(pred_diss) <- NULL
    return(pred_diss)
  }
}



