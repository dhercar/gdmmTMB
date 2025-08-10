#' Fit Generalized Dissimilarity Mixed Model (GDMM)
#'
#' @param Y Response matrix (e.g., species x site matrix). Either Y or Y_diss must be provided.
#' @param Y_diss Pre-calculated dissimilarity values between site pairs. If provided, D must also be supplied.
#' @param Y_den Denominator values for Jaccard, Sorensen or Bray-Curtis dissimilarities when using Y_diss. Required when family = "binomial" and Y_diss is provided.
#' @param D Data frame or matrix with two columns specifying site pairs (indices) corresponding to Y_diss values. Only required when Y_diss is provided.
#' @param X Data frame containing predictor variables as site-level values.
#' @param X_pair Data frame containing predictor variables as pairwise distances.
#' @param diss_formula Formula specifying predictors for dissimilarity gradients. All variables must be numeric. `isp(c)` can be used in combination with `mono = TRUE` to fit monotonic I-splines (e.g., `~ elevation + temperature`).
#' @param uniq_formula Formula specifying predictors and random effects for uniqueness. Uses lme4-style syntax for random effects (e.g., `~ treatment + (1|site)`).
#' @param mono Logical. If TRUE, enforces monotonic (non-decreasing) dissimilarity effects. Default is FALSE.
#' @param family Distribution family for the response. One of "normal", "binomial", or "beta". Default is "normal".
#' @param link Link function. If NULL (default), automatically chosen based on family: "identity" for normal, "logit" for binomial/beta.
#' @param scale_diss Numeric vector of length 2 specifying the range of the re-scaling of dissimilarity values.
#' @param binary Logical. Whether to treat response as binary data before calculating dissimilarity from Y. Default is FALSE.
#' @param method Dissimilarity method applied to Y. One of "bray", "sorensen", "jaccard". Default is "bray".
#' @param control List of control parameters passed to nlminb optimizer (e.g., `control = list(rel.tol = 1e-8, iter.max = 500)`).
#' @param trace Logical. If TRUE, prints optimization  information. Default is FALSE.
#' @param bboot Logical. If TRUE, performs Bayesian bootstrapping. Default is FALSE.
#' @param n_boot Integer. Number of bootstrap samples when `bboot = TRUE`. Default is 1000.
#' @param n_cores Integer. Number of cores for parallel processing during bootstrapping. If NULL, uses `detectCores() - 2`.
#'
#' @return For `bboot = FALSE`: An object of class "gdmm" containing model fit, parameters, and diagnostics.
#'         For `bboot = TRUE`: An object of class "bbgdmm" containing bootstrap samples and model information.
#' @export
#'
#' @importFrom hardhat mold default_formula_blueprint
#' @importFrom lme4 nobars
#' @importFrom Matrix sparseMatrix sparse.model.matrix
#' @importFrom TMB MakeADFun sdreport
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom gtools rdirichlet
#' @importFrom stats nlminb terms
#' @useDynLib gdmmTMB, .registration = TRUE

gdmm <- function(Y = NULL,
                 Y_diss = NULL,
                 Y_den = 0,
                 D = NULL,
                 X = NULL,
                 X_pair = NULL,
                 diss_formula = NULL,
                 uniq_formula = NULL,
                 pair_formula = NULL,
                 mono = FALSE,
                 mono_pair = FALSE,
                 family = 'normal',
                 link = NULL,
                 scale_diss = NULL,
                 binary = FALSE,
                 method = 'bray',
                 control = NULL,
                 trace = FALSE,
                 bboot = FALSE,
                 n_boot = 1000,
                 n_cores = NULL) {

  #### Checks ####
    if (!is.null(pair_formula) && is.null(X_pair)) {
    stop("You supplied 'pair_formula' but no 'X_pair'. Please supply a pairwise data.frame/matrix aligned to D/Y_diss, or remove 'pair_formula'.")
  }

  # Y_den present if needed
  if (!is.null(Y_diss) && family == "binomial" && (is.null(Y_den) || Y_den == 0)) {
    stop("When using 'Y_diss' with binomial family, 'Y_den' must be provided")
  }

  # Check Y format
  if (!is.null(Y)) {
    if (!is.matrix(Y) && !is.data.frame(Y)) {
      stop("'Y' must be a matrix or data frame")
    }
  }

  if (!(is.null(scale_diss) || (length(scale_diss) == 2))) {
    stop("'scale_diss' must be a vector of length 2")
  }

  # Check family
  valid_families <- c("normal", "binomial", "beta")
  if (!family %in% valid_families) {
    stop(paste("'family' must be one of:", paste(valid_families, collapse = ", ")))
  }

  # Check correct method for binomial
  if (family == "binomial" && !method %in% c("bray", "sorensen", "jaccard")) {
    stop("Binomial family only valid with 'bray', 'sorensen', or 'jaccard' methods")
  }

  # Check numeric variables
  if (!is.null(X) && !is.null(diss_formula)) {
    vars_form <- all.vars(diss_formula)
    is_numeric_X <- apply(X[,vars_form, drop = F], 2, is.numeric)
    if (!all(is_numeric_X)) {
      stop(paste0("Variables included in 'diss_formula' are not numeric: ", paste(vars_form[!is_numeric_X], collapse = ', ')))
    }
  }

  if (!is.null(X_pair) && !is.null(pair_formula)) {
    vars_form <- all.vars(pair_formula)
    is_numeric_X_pair <- apply(X_pair[,vars_form, drop = F], 2, is.numeric)
    if (!all(is_numeric_X_pair)) {
      stop(paste0("Variables included in 'pair_formula' are not numeric: ", paste(vars_form[!is_numeric_X_pair], collapse = ', ')))
    }
  }

  #### Model specs ####

  # Match family with link
  if (is.null(link))  {
    link <- switch(family,
                   'normal' = 'identity',
                   'binomial' = 'logit',
                   'beta' = 'logit')
  }

  # Check link
  valid_link <- c("identity", "logit", "neg_exp", 'neg_gaus')
  if (!link %in% valid_link ) {
    stop(paste("'family' must be one of:", paste(valid_link, collapse = ", ")))
  }

  # Numeric codes for fam and link
  family_num <- switch(family,
                       'normal' = 0,
                       'binomial' = 1,
                       'beta' = 2)

  link_num <- switch(link,
                     'identity' = 0,
                     'logit' = 1,
                     'neg_exp' = 2,
                     'neg_gaus' = 3)

  #### VARS PREP ####

  # 1. Diss. model matrix (pair)
  if (!is.null(pair_formula) & !is.null(X_pair)) {
  form_X_pair <- hardhat::mold(pair_formula,
                          X_pair,
                          blueprint = default_formula_blueprint(intercept = FALSE))
  } else {
  form_X_pair <- list(predictors =  matrix(ncol = 0, nrow = 0))
  }

  # 2. Diss. model matrix (site)
  if (!is.null(diss_formula) & !is.null(X)) {
    form_X <- hardhat::mold(diss_formula,
                            X,
                            blueprint = default_formula_blueprint(intercept = FALSE))
  } else {
    form_X <- list(predictors =  matrix(ncol = 0, nrow = 0))
  }

  # 3. Uniq. model matrix
  if (!is.null(uniq_formula) & !is.null(X)) {
    fix_terms_W <- lme4::nobars(uniq_formula)
    re_vars <- all.vars(uniq_formula)[!(all.vars(uniq_formula) %in% all.vars(fix_terms_W))]

    if ( length(attr(terms(fix_terms_W), 'variables')) > 1 ) {
      form_W <- hardhat::mold(fix_terms_W,
                              X,
                              blueprint = default_formula_blueprint(intercept = FALSE, indicators = 'traditional'))
    } else {
      form_W <- list(predictors =  matrix(ncol = 0, nrow = 0))
    }

  } else {
    re_vars <- character(0)
    form_W <- list(predictors =  matrix(ncol = 0, nrow = 0))
  }

  # 4. Re model matrix
  map = list()
  if (length(re_vars) > 0) {
    Z_design <- Matrix::sparse.model.matrix(eval(parse(text = paste0('~ 0 + ', paste0(re_vars, collapse = ' + ')))),
                             data = X)
    has_re = 1
    map_re <- as.numeric(as.factor(unlist(lapply(re_vars, function(x) rep(x, length(unique(X[[x]])))))))

  } else {
    n = length(unique(c(D[,1], D[,2])))
    Z_design <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(n, 1))
    has_re = 0
    map_re = 0
  }

  # 5. Response
  if (is.null(Y_diss)) {
    if (family %in% c('binomial')) {
      Y_pair <- make_y_df(com = Y, D = D, method = method, binary = binary, num_den = TRUE)
      Y_diss <- Y_pair[,3]
      Y_den <- Y_pair[,4]
      D <- Y_pair[,1:2]
      map = list(log_scale = factor(NA))

    } else {
      Y_pair <- make_y_df(com = Y, D = D, method = method, binary = binary, num_den = FALSE)
      Y_diss <- Y_pair[,3]
      Y_den <- 0
      D <- Y_pair[,1:2]
    }
  }

  if(!is.null(scale_diss)){
    Y_diss <- scale_dist(Y_diss, scale_diss)
  }

  #### construct model ####
  data <- list(
    Y = Y_diss,
    Y_den = Y_den,
    D = as.matrix(D - 1), # remove 1 to match c++
    X = as.matrix(form_X$predictors),
    X_pair =  as.matrix(form_X_pair$predictors),
    W = as.matrix(form_W$predictors),
    Z = Z_design,
    has_random = has_re,
    map_re = map_re - 1,  # remove 1 to match c++ indexing
    mono = as.numeric(mono),
    mono_pair = as.numeric(mono_pair),
    link = link_num,
    family = family_num)

  parameters <- list(
    intercept = 0,
    beta = rep(0, ncol(form_X$predictors)),
    beta_p = rep(0, ncol(form_X_pair$predictors)),
    lambda = rep(0, ncol(form_W$predictors)),
    log_sigma_re = rep(0, length(re_vars)),
    u = rep(0, ncol(Z_design)*has_re),
    log_scale = 0)

  control_def <- list(rel.tol = 1e-10,
                      eval.max = 1000,
                      iter.max = 1000,
                      trace = 0)

  if (!is.null(control)) {
    for (nm in names(control)) {
      control_def[[nm]] <- control[[nm]]
    }
  }

  #### NO BAYESIAN BOOTSTRAPPING ####
  if (bboot == FALSE) {
    data[['weights']] <- rep(1, length(Y_diss))

    obj <- TMB::MakeADFun(
      data = data,
      parameters = parameters,
      DLL = 'gdmmTMB',
      map = map,
      random = 'u',
      silent = !trace
    )

    # not in use (using exp transformation in c++)
    n_par <- length(obj$par)
    lower_bounds <- rep(-Inf, n_par)
    upper_bounds <- rep(Inf, n_par)

    if (mono) {
      #lower_bounds[names(obj$par) == 'beta'] <- 0
    }

    opt <- nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      lower = lower_bounds,
      upper = upper_bounds,
      control = control_def)

      print(opt$message)

      fit <- list(Y = Y,
                  Y_diss = Y_diss,
                  Y_den = Y_den,
                  X = X,
                  X_pair = X_pair,
                  D = D,
                  form_X = form_X,
                  form_W = form_W,
                  form_X_pair = form_X_pair,
                  Z_design = Z_design,
                  map_re = map_re,
                  re_vars = re_vars,
                  diss_formula = diss_formula,
                  uniq_formula = uniq_formula,
                  pair_formula = pair_formula,
                  link = link,
                  family = family,
                  obj = obj,
                  opt = opt,
                  boot = FALSE,
                  mono = mono,
                  mono_pair = mono_pair,
                  call = match.call())

      out <- new_gdmm(fit)
      return(out)

  #### BAYESIAN BOOTSTRAPPING ####
  } else if (bboot == TRUE) {
    if (is.null(n_cores))  n_cores = parallel::detectCores(logical = TRUE) - 2

    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))


    cat('Running bayesian bootstrapping on', n_cores, 'cores ...')

    results <- foreach(i = 1:n_boot, .combine = rbind,
                       .packages = c('TMB', 'gtools')) %dopar% {
      # calculate random weights
      n <- length(unique(c(D[,1], D[,2])))
      weights_sample <-  c(gtools::rdirichlet(1,rep(1,n)))
      weights_pair <- (weights_sample[D[,1]] * weights_sample[D[,2]])
      data[['weights']] <- weights_pair/sum(weights_pair)*length(Y_diss) # standardise weights
      # objective
      obj <- MakeADFun(
        data = data,
        parameters = parameters,
        DLL = 'gdmmTMB',
        map = map,
        random = 'u',
        silent = !trace
      )

      n_par <- length(obj$par)
      lower_bounds <- rep(-Inf, n_par)
      upper_bounds <- rep( Inf, n_par)

      if (mono) {
        # lower_bounds[names(obj$par) == 'beta'] <- 0
      }

      opt <- nlminb(
        start = obj$par,
        objective = obj$fn,
        gradient = obj$gr,
        lower = lower_bounds,
        upper = upper_bounds,
        control = control_def)

      c(sdreport(obj)$value, u_re_ = obj$report()$u_0, obj$report()$sigma_re, logLikelihood = -opt$objective)
    }

    cat('done\n')
    fit <- list(Y = Y,
                Y_diss = Y_diss,
                Y_den = Y_den,
                X = X,
                X_pair = X_pair,
                D = D,
                form_X = form_X,
                form_W = form_W,
                form_X_pair = form_X_pair,
                Z_design = Z_design,
                map_re = map_re,
                re_vars = re_vars,
                diss_formula = diss_formula,
                uniq_formula = uniq_formula,
                pair_formula = pair_formula,
                link = link,
                family = family,
                boot_samples = results,
                n_boot = n_boot,
                boot = TRUE,
                mono = mono,
                mono_pair = mono_pair,
                call = match.call())
    out <- new_bbgdmm(fit)
    return(out)
  }
}


#' Print method for gdmm objects
#'
#' @description
#' Prints basic information from a fitted gdmm object
#'
#' @param m A fitted model object of class \code{gdmm}
#' @seealso \code{\link{print.bbgdmm}} \code{\link{gdmm}}
#'
#' @method print gdmm

print.gdmm <- function(m){
  print_title('Generalized dissimilarity mixed model (GDMM)', symb = '-')  # call

  cat('convergence:', ifelse(m$opt$convergence == 0, 'succesful convergence\n', 'optimization has not reached succesful convergence\n'))
  cat('nlminb message: ')
  cat(m$opt$message, '\n\n')

  cat('function call:\n')
  print(m$call)
  cat('\n')

  print_title2(' Details ', symb = '-')
  if (m$mono) cat('- monotonic dissimilarity effects (mono = TRUE) \n')
  cat('- distribution (family):', m$family, '\n')
  cat('- link:', m$link, '\n\n')

  print_title2(' Predictors ', symb = '-')

  cat('- dissimilarity gradient(s): ')
  if (!is.null(m$diss_formula)) {
    print(m$diss_formula)
  } else {
    cat('no dissimilarity gradients included \n')
  }

  cat('- effect(s) on uniqueness: ')
  if (!is.null(m$uniq_formula)) {
    print(lme4::nobars(m$uniq_formula))
  } else {
    cat('no predictors on uniqueness included \n')
  }

  cat('- random effects: ')
  cat(ifelse(length(m$re_vars) > 0, paste(m$re_vars,sep = ' ,'), 'no random effects included \n'))

}


#' Print method for bbgdmm objects
#'
#' @description
#' Prints basic information from a fitted bbgdmm object
#'
#' @param m A fitted model object of class \code{bbgdmm}
#' @seealso \code{\link{print.gdmm}} \code{\link{gdmm}}
#'
#' @method print bbgdmm

print.bbgdmm <- function(m){
  print_title('Generalized dissimilarity mixed model (GDMM)', symb = '-')  # call

  cat('Bayesian Bootstrapping with', m$n_boot, 'random samples\n\n')

  cat('function call:\n')
  print(m$call)
  cat('\n')

  print_title2(' Details ', symb = '-')
  if (m$mono) cat('- monotonic dissimilarity effects (mono = TRUE) \n')
  cat('- distribution (family):', m$family, '\n')
  cat('- link:', m$link, '\n\n')

  print_title2(' Predictors ', symb = '-')

  cat('- dissimilarity gradient(s): ')
  if (!is.null(m$diss_formula)) {
    print(m$diss_formula)
  } else {
    cat('none included \n')
  }

  cat('- effect(s) on uniqueness: ')
  if (!is.null(m$uniq_formula)) {
    print(lme4::nobars(m$uniq_formula))
  } else {
    cat('none included \n')
  }

  cat('- random effects: ')
  cat(ifelse(length(m$re_vars) > 0, paste(m$re_vars,sep = ' ,'), 'none included \n'))

  cat('\n')
  print_title2('', symb = '-')
}



