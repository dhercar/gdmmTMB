# Predict method for gdmm and bbgdmm objects

Generate predictions from fitted generalized dissimilarity mixed models.
This method can predict dissimilarity or uniqueness for new data, with
optional confidence intervals. Works with both regular gdmm and
bootstrap bbgdmm objects.

## Usage

``` r
# S3 method for class 'gdmm'
predict(
  object,
  new_X = NULL,
  new_W = NULL,
  new_X_pair = NULL,
  new_re = NULL,
  new_D = NULL,
  re_sd = logical(0),
  component = "dissimilarity",
  scale_uniq = TRUE,
  type = "response",
  CI = TRUE,
  CI_quant = c(0.95, 0.5),
  n_sim = NULL,
  n_cores = 1,
  ...
)

# S3 method for class 'bbgdmm'
predict(
  object,
  new_X = NULL,
  new_W = NULL,
  new_X_pair = NULL,
  new_re = NULL,
  new_D = NULL,
  re_sd = logical(0),
  component = "dissimilarity",
  scale_uniq = TRUE,
  type = "response",
  CI = TRUE,
  CI_quant = c(0.95, 0.5),
  n_sim = NULL,
  n_cores = 1,
  ...
)
```

## Arguments

- object:

  A fitted model object of class 'gdmm' or 'bbgdmm'

- new_X:

  A data frame or matrix of new predictor variables for the
  dissimilarity component. If NULL (default), uses the original data
  from the fitted model.

- new_W:

  A data frame or matrix of new predictor variables for the uniqueness
  component.

- new_X_pair:

  A data frame or matrix of new pairwise predictor variables for the
  dissimilarity component, usually distances.

- new_re:

  A data frame of new random effect variables. If NULL, uses random
  effects in the fitted model

- new_D:

  Matrix of pairwise combination indices. Required when new_X_pair is
  provided.

- re_sd:

  A character indicating which random effects should be treated as
  having standard deviation (used for simulation).

- component:

  Character string specifying what to predict. Either "dissimilarity"
  (default) for pairwise dissimilarities or "uniqueness" for site-level
  uniqueness values.

- scale_uniq:

  Logical indicating whether uniqueness values should be transformed
  into proportions as in LCBD. Default is TRUE.

- type:

  Character string specifying the type of prediction. Either "response"
  (default) for predictions on the response scale, or "link" for
  predictions on the linear predictor scale.

- CI:

  Logical indicating whether to compute confidence intervals. Default is
  TRUE.

- CI_quant:

  Numeric vector of confidence levels to compute. Default is
  `c(0.95, 0.5)` for 95% and 50% confidence intervals.

- n_sim:

  Integer specifying the number of simulations for confidence intervals.
  If NULL, uses 1000 for 'gdmm' objects and the number of bootstrap
  samples for 'bbgdmm' objects.

- n_cores:

  Integer specifying the number of cores to use for parallel
  computation. Currently not implemented. Default is 1.

- ...:

  Additional parameters (not used).

## Value

If `CI = FALSE`, returns a numeric vector of predictions. If
`CI = TRUE`, returns a matrix with columns for the mean prediction and
confidence interval bounds, named according to the confidence levels
specified in `CI_quant`.

## Details

This function generates predictions from fitted gdmm models by:

- Computing expected values using model coefficients

- Optionally generating confidence intervals through bootstrapping

For 'gdmm' objects, confidence intervals are computed using parametric
bootstrapping.

When `component = "uniqueness"`, the function computes site-level
uniqueness values based on the full dissimilarity matrix. When
`component = "dissimilarity"`, it returns pairwise dissimilarity
predictions. If site-level random effects are included via
`uniq_formula` when fitting the model (e.g.,`uniq_formula = (1|site)`),
confidence intervals can incorporate site-level uncertainty via
`re_sd = c('site')`.

## Functions

- `predict(bbgdmm)`: Method for bootstrap gdmm objects
