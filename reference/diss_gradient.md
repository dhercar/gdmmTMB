# Predicted Effects of Environmental Gradients on Dissimilarity

Computes the predicted effect of one or more covariates on compositional
dissimilarity from a fitted `gdmm` or `bbgdmm` model. Optionally,
returns uncertainty bands (credible or confidence intervals) using
either bootstrapped samples or parametric bootstrapping.

## Usage

``` r
diss_gradient(
  m,
  var = "all",
  n = 100,
  CI = TRUE,
  n_sim = NULL,
  CI_quant = c(0.95)
)
```

## Arguments

- m:

  A fitted model object of class `gdmm` or `bbgdmm`, representing a
  generalized dissimilarity (mixed) model.

- var:

  Character vector of variable names to evaluate. Use `"all"` (default)
  to compute for all predictors in the model.

- n:

  Integer. Number of points to evaluate along each gradient. Defaults to
  100.

- CI:

  Logical. Whether to compute confidence intervals. Default is `TRUE`.

- n_sim:

  Number of draws used to obtain CI. If `NULL`, uses `m$n_boot` for
  `bbgdmm` or 1000 for `gdmm`.

- CI_quant:

  Numeric vector specifying the quantile coverage (e.g., `0.95` for 95%
  CI). Default is `0.95`.

## Value

A named list of data frames, one for each variable in `var`. Each data
frame contains:

- `var`: Variable names.

- `f_x`: Predicted value after fitted non-linear transformation (e.g.,
  monotonic I-splines).

- `x`: The values of the focal gradient.

- `CI lower / CI upper`: (Optional) Lower and upper bounds of the
  confidence interval.

## Details

The function generates predicted dissimilarities along a gradient of
each specified predictor, while holding all other predictors constant.
If `CI = TRUE`, it uses bootstrap samples (for `bbgdmm`) or simulates
from the joint covariance matrix (for `gdmm`) to construct confidence or
credible intervals.

This approach is useful for visualizing directional compositional
changes along individual environmental gradients.
