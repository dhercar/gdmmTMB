#' Summary method for gdmm objects
#'
#' S3 method for summarizing gdmm model fits.
#'
#' @param object An object of class \code{gdmm}.
#' @param ... additional arguments (not used)
#' @return An object of class \code{summary.gdmm}.
#' @method summary gdmm
#' @export
summary.gdmm <- function(object, ...) {

    # Get coefficient table with p-values
    sdr <- TMB::sdreport(object$obj)
    coef_table <- TMB::summary.sdreport(sdr, select = c('report'), p.value = TRUE)

    # Add significance stars
    p_values <- coef_table[, "Pr(>|z^2|)"]
    sig_stars <- ifelse(p_values < 0.001, "***",
                        ifelse(p_values < 0.01, "**",
                               ifelse(p_values < 0.05, "*",
                                      ifelse(p_values < 0.1, ".", " "))))

    # Create data.frame with lm-style column names and add significance stars
    result <- as.data.frame(coef_table)


    # Calculate AIC
    m_AIC <- AICc(log_likelihood = -object$opt$objective,
                  n = length(object$Y_diss),
                  k = length(object$opt$par))
    m_BIC <- BIC2(log_likelihood = -object$opt$objective,
                  n = length(object$Y_diss),
                  k = length(object$opt$par))

    out <- list(call = object$call,
                mono = object$mono,
                mono_pair = object$mono_pair,
                table = coef_table,
                p_values = p_values,
                names_beta =  paste0('diss: ', colnames(object$form_X$predictors)),
                names_beta_p = paste0('diss(p): ', colnames(object$form_X_pair$predictors)),
                names_lambda =  paste0('uniq: ', colnames(object$form_W$predictors)),
                sig = sig_stars,
                AIC = m_AIC$AIC,
                AICc = m_AIC$AICc,
                BIC = m_BIC,
                logL = -object$opt$objective)
    class(out) <- 'summary.gdmm'
    return(out)
}

#' Summary method for bbgdmm objects
#'
#' S3 method for summarizing bbgdmm model fits.
#'
#' @param object An object of class \code{bbgdmm}.
#' @param quantiles A numeric vector of quantiles to display.
#' @param null_value Value used for pseudo p-value tests.
#' @param ... additional arguments (not used)
#'
#' @return An object of class \code{summary.bbgdmm}.
#' @method summary bbgdmm
#' @importFrom stats sd
#' @export
summary.bbgdmm <- function(object,
                         quantiles = c(0.025, 0.5, 0.975),
                         null_value = 0, ...) {

  # Basic statistics
  samples <- object$boot_samples[, (colnames(object$boot_samples) != 'logLikelihood') & (substr(colnames(object$boot_samples), 1,5) != 'u_re_')]

  logL <- mean(object$boot_samples[,'logLikelihood'])

  # Calculate AIC
  m_AIC <- AICc(log_likelihood = logL,
                n = length(object$Y_diss),
                k = ncol(samples))
  m_BIC <- BIC2(log_likelihood = logL,
                n = length(object$Y_diss),
                k = ncol(samples))

  means <- apply(samples, 2, mean)
  sds <- apply(samples, 2, stats::sd)

  # Quantiles
  CI <- t(apply(samples,2, function(x) quantile(x, probs = quantiles)))

  # Pseudo Z
  pseudo_z <- (means - null_value) / sds

  # Pseudo P
  pseudo_p <- apply(samples, 2, function(x) {
    n <- length(x)
    obs_mean <- mean(x)
    side <- sign(obs_mean - null_value)
    return(sum( (x - null_value)*side < null_value) / n)
  })


  sig_stars <- ifelse(pseudo_p < 0.001, "***",
                      ifelse(pseudo_p < 0.01, "**",
                             ifelse(pseudo_p < 0.05, "*",
                                    ifelse(pseudo_p < 0.1, ".", " "))))

  out <- list(call = object$call,
              mono = object$mono,
              mono_pair = object$mono_pair,
              estimate = means,
              sds = sds,
              pseudo_zval = pseudo_z,
              pseudo_pval = pseudo_p,
              CI = CI,
              n_boot = object$n_boot,
              quantiles = quantiles,
              AIC = m_AIC$AIC,
              AICc = m_AIC$AICc,
              BIC = m_BIC,
              names_beta =  paste0('diss: ', colnames(object$form_X$predictors)),
              names_beta_p = paste0('diss(p): ', colnames(object$form_X_pair$predictors)),
              names_lambda =  paste0('uniq: ', colnames(object$form_W$predictors)),
              sig = sig_stars,
              logL = logL)

  class(out) <- 'summary.bbgdmm'
  return(out)
}

#' Print method for summary.gdmm
#'
#' Internal print method for objects of class \code{summary.gdmm}.
#'
#' @param x An object of class \code{summary.gdmm}.
#' @param ... Additional parameters (not used)
#' @method print summary.gdmm
#' @export
#' @noRd
print.summary.gdmm <- function(x, ...) {
  print_title('GDMM SUMMARY', symb = '-')  # call
  cat("call:\n\n")
  call_str <- deparse(x$call)
  cat(" ", call_str, "\n\n", sep = '')

  # table
  cat("coeff. summary table:\n\n")
  rownames(x$table)[rownames(x$table) == 'e_beta'] <- x$names_beta
  rownames(x$table)[rownames(x$table) == 'e_beta_p'] <- x$names_beta_p
  rownames(x$table)[rownames(x$table) == 'lambda'] <- x$names_lambda
  rownames(x$table)[rownames(x$table) == 'intercept'] <- '(Intercept)'
  x$table <- data.frame(x$table, x$sig)
  names(x$table) <-  c("Estimate", "Std.Err.", "z value", "Pr(>|z^2|)", " ")
  print(x$table, digits = 3)

  if (x$mono) message('\nCAUTION: P-values for dissimilarity components (diss) should not be used for hypothesis testing when "mono = TRUE"')
  if (x$mono_pair) message('\nCAUTION: P-values for dissimilarity components (diss(p)) should not be used for hypothesis testing when "mono_pair = TRUE"')

  cat("---\n")
  cat("signif. codes: '***' <0.001 '**' <0.01 '*' <0.05 '.' <0.1 \n")
  cat("---\n\n")
  # logL, AIC, AICc, and BIC
  cat("Marginal log-likelihood: ", x$logL, "\nAIC: ", x$AIC, ", AICc: ", x$AICc, ", BIC: ", x$BIC, "\n")
  cat('\n')
  print_title2('', symb = '-')
}

#' Print method for summary.bbgdmm
#'
#' Internal print method for objects of class \code{summary.bbgdmm}.
#'
#' @param x An object of class \code{summary.bbgdmm}.
#' @param ... Additional parameters (not used)
#' @method print summary.bbgdmm
#' @export
#' @noRd
print.summary.bbgdmm <- function(x, ...) {
  print_title('BBGDMM SUMMARY', symb = '-')  # call
  cat("call:\n\n")
  call_str <- deparse(x$call)
  cat(" ", call_str, "\n\n", sep = '')

  # table
  print_title2(" Coeff. Table ", symb = '-')
  rownames(x$CI)[rownames(x$CI) == 'e_beta'] <- x$names_beta
  rownames(x$CI)[rownames(x$CI) == 'e_beta_p'] <- x$names_beta_p
  rownames(x$CI)[rownames(x$CI)  == 'lambda'] <- x$names_lambda
  rownames(x$CI)[rownames(x$CI)  == 'intercept'] <- '(Intercept)'
  colnames(x$CI) <-  paste0(x$quantiles*100, '%')

  table <- data.frame(x$estimate,
                      x$sds,
                      x$CI,
                      x$pseudo_zval,
                      ifelse(x$pseudo_pval == 0, paste('<', 1/x$n_boot), x$pseudo_pval),
                      x$sig)
  names(table) <- c('Estimate', 'Std.Err.', paste0(x$quantiles*100,'%'), 'pseudo-Z', 'pseudo-Pval', ' ')
  print(table, digits = 4)

  if (x$mono) message('\nCAUTION: P-values for dissimilarity components (diss) should not be used for hypothesis testing when "mono = TRUE"')
  if (x$mono_pair) message('\nCAUTION: P-values for dissimilarity components (diss(p)) should not be used for hypothesis testing when "mono_pair = TRUE"')

  cat("---\n")
  cat("signif. codes: '***' <0.001 '**' <0.01 '*' <0.05 '.' <0.1 \n")
  cat("---\n\n")
  cat("Marginal log-likelihood (average): ", x$logL, "\nAIC: ", x$AIC, ", AICc: ", x$AICc, ", BIC: ", x$BIC, "\n")
  cat('\n')
  print_title2('', symb = '-')
}

