# Summary method for bbgdmm objects

S3 method for summarizing bbgdmm model fits.

## Usage

``` r
# S3 method for class 'bbgdmm'
summary(object, quantiles = c(0.025, 0.5, 0.975), null_value = 0, ...)
```

## Arguments

- object:

  An object of class `bbgdmm`.

- quantiles:

  A numeric vector of quantiles to display.

- null_value:

  Value used for pseudo p-value tests.

- ...:

  additional arguments (not used)

## Value

An object of class `summary.bbgdmm`.
