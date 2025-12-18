# Ispline Basis with Degree 2 by Default

A thin wrapper around
[`splines2::isp()`](https://wwenjie.org/splines2/reference/iSpline.html).
Uses `degree = 2` unless otherwise specified. All other arguments are
passed through.

## Usage

``` r
isp(x, degree = 2, ...)
```

## Arguments

- x:

  A numeric vector.

- degree:

  Spline degree; defaults to 2 here (instead of 3).

- ...:

  additional parameters passed to
  [`splines2::isp`](https://wwenjie.org/splines2/reference/iSpline.html)
