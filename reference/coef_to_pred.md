# Convert model coefficients to predictions

Internal function that transforms model parameters and new data into
dissimilarity or uniqueness predictions.

## Usage

``` r
coef_to_pred(
  object,
  intercept,
  beta,
  beta_p,
  lambda,
  u,
  new_W,
  new_X,
  new_X_pair,
  new_re,
  D,
  n,
  component,
  scale_uniq,
  type,
  sigma
)
```

## Arguments

- object:

  Model object containing fitted model components

- intercept:

  Numeric, model intercept

- beta:

  Numeric vector, coefficients for dissimilarity component

- beta_p:

  Numeric vector, coefficients for dissimilarity component (pairwise)

- lambda:

  Numeric vector, coefficients for uniqueness component

- u:

  Numeric vector, random effect values

- new_W:

  Matrix/data.frame, new data for uniqueness component

- new_X:

  Matrix/data.frame, new data for dissimilarity component

- new_X_pair:

  Matrix/data.frame, new data for dissimilarity component (pairwise)

- new_re:

  Matrix/data.frame, new random effect data

- D:

  Matrix, pairwise combination indices

- n:

  Integer, number of observations

- component:

  Character, "dissimilarity" or "uniqueness"

- scale_uniq:

  Logical, whether to scale uniqueness values

- type:

  Character, "response" or "link" scale

- sigma:

  Numeric, random effect standard deviations

## Value

Numeric vector of predictions (dissimilarities or uniqueness values)
