# Scale a vector into an specified

Given a numeric vector `x`, this function re-scales `x` into the
interval \[`min`, `max`\].

## Usage

``` r
scale_dist(x, new = c(0.01, 0.99))
```

## Arguments

- x:

  A numeric vector

- new:

  Numeric scalar: the upper nad lower bound of the target interval.
  Default is \[`0.01`, `0.99`\].

## Value

A numeric vector of the same length as `x`, with values in the interval
\[`min(new)`, `max(new)`\].
