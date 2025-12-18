# Fit Generalized Dissimilarity Mixed Model (GDMM).

Fits a generalized dissimilarity mixed model (gdmm) for ecological
community data. The model can estimate the effect of predictors on
community dissimilarity and uniqueness simultaneously.

The function supports two options for estimation:

- Hierarchical structure via site-level random effects in `uniq_formula`
  (e.g., `~ ... + (1|site)`).

- Bayesian bootstrapping (`bboot = TRUE`) to obtain site-level
  uncertainty in parameter estimates.

Dissimilarity gradients are specified in `diss_formula`, which may
include non-linear transformations such as monotonic I-splines (`isp(x)`
with `mono = TRUE`) to allow flexible monotonic responses.

Direct effects on community uniqueness are supported via `uniq_formula`.

## Usage

``` r
gdmm(
  Y = NULL,
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
  family = "normal",
  link = NULL,
  scale_diss = NULL,
  binary = FALSE,
  method = "bray",
  control = NULL,
  trace = FALSE,
  bboot = FALSE,
  n_boot = 1000,
  n_cores = NULL
)
```

## Arguments

- Y:

  Response matrix (e.g., site x species matrix). Either Y or Y_diss must
  be provided.

- Y_diss:

  Pre-calculated dissimilarity values between site pairs. If provided, D
  must also be supplied. Either Y or Y_diss must be provided.

- Y_den:

  Denominator values for Jaccard, Sørensenensen or Bray-Curtis
  dissimilarities when using Y_diss. Required when family = "binomial"
  and Y_diss is provided.

- D:

  Data frame or matrix with two columns specifying site pairs (indices)
  corresponding to Y_diss values. Only required when Y_diss is provided.

- X:

  Data frame containing predictor variables as site-level values.

- X_pair:

  Data frame containing predictor variables as pairwise distances.

- diss_formula:

  Formula specifying predictors for dissimilarity gradients. All
  variables must be numeric.
  [`isp()`](https://dhercar.github.io/gdmmTMB/reference/isp.md) can be
  used in combination with `mono = TRUE` to fit monotonic I-splines
  (e.g., `~ isp(elevation) + isp(temperature)`).

- uniq_formula:

  Formula specifying predictors and random effects for uniqueness. Uses
  lme4-style syntax for random effects (e.g., `~ treatment + (1|site)`).

- pair_formula:

  Formula specifying predictors for dissimilarity gradients. Variables
  must be numeric.

- mono:

  Logical. If `TRUE`, enforces monotonic (non-decreasing) dissimilarity
  effects on predictors included in `diss_formula`. Default is `FALSE`.

- mono_pair:

  Logical. If `TRUE`, enforces monotonic (non-decreasing) dissimilarity
  effects on predictors included in `pair_formula`. Default is `FALSE`.

- family:

  Distribution family for the response. One of `"normal"`, `"binomial"`,
  or `"beta"`. Default is `"normal"`.

- link:

  Link function. If `NULL` (default), automatically chosen based on
  family: `"identity"` for normal, `"logit"` for binomial/beta.

- scale_diss:

  Numeric vector of length 2 specifying the range of the re-scaling of
  dissimilarity values. Useful when using a beta distribution if some
  dissimilarities are exactly 0 and/or 1.

- binary:

  Logical. Whether to treat response as binary data before calculating
  dissimilarities from Y. Default is `FALSE`.

- method:

  Dissimilarity method applied to Y. See
  [[`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)](https://rdrr.io/cran/vegan/man/vegdist.html)
  for a list of compatible methods. Default is `"bray"`.

- control:

  List of control parameters passed to nlminb optimizer (e.g.,
  `control = list(rel.tol = 1e-8, iter.max = 500)`).

- trace:

  Logical. If `TRUE`, prints information during optimization. Default is
  `FALSE`.

- bboot:

  Logical. If `TRUE`, performs Bayesian bootstrapping. Default is
  `FALSE`.

- n_boot:

  Integer. Number of bootstrap samples when `bboot = TRUE`. Default is
  1000.

- n_cores:

  Integer. Number of cores for parallel processing during bootstrapping.
  If `NULL`, uses `detectCores() - 2`.

## Value

For `bboot = FALSE`: An object of class "gdmm" containing model fit,
parameters, and diagnostics. For `bboot = TRUE`: An object of class
"bbgdmm" containing bootstrap samples and model information.
