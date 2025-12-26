# Compute Sum of Squares or LCBD from a Symmetric Dissimilarity Matrix

Computes the per-site sum of squares (SS) from a symmetric dissimilarity
matrix using double-centering. Optionally returns LCBD values (relative
contributions to total SS).

## Usage

``` r
SS_calc(x, LCBD = FALSE)
```

## Arguments

- x:

  A numeric symmetric matrix (e.g., a squared Euclidean dissimilarity
  matrix).

- LCBD:

  Logical; if `TRUE`, returns normalized values summing to 1 (LCBD).

## Value

A numeric vector of length `nrow(x)`. If `LCBD = TRUE`, the values sum
to 1.

## Details

The function double-centers the dissimilarity matrix `x` and extracts
the diagonal elements. If `LCBD = TRUE`, the diagonal is normalized to
sum to 1.

## See also

[`LCBD.comp`](http://adeverse.github.io/adespatial/reference/LCBD.comp.md)
