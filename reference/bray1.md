# Bray-Curtis dissimilarity components

Calculate the numerator and denominator components of the Bray-Curtis
dissimilarity index. The Bray-Curtis dissimilarity is calculated as
\\\frac{\sum\_{i} \|x_i - y_i\|}{\sum\_{i} (x_i + y_i)}\\

## Usage

``` r
bray1(x, y)
```

## Arguments

- x:

  numeric vector of the first sample to be compared

- y:

  numeric vector of the second sample to be compared

## Value

A named numeric vector with two elements:

- num - numerator: sum of absolute differences between x and y

- den - denominator: sum of x and y values
