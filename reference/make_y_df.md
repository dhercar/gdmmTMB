# Format community data

Formats species x site table into a long data frame of community
dissimilarities

## Usage

``` r
make_y_df(
  com = NULL,
  D = NULL,
  id = NULL,
  method = "bray",
  num_den = FALSE,
  trans = NULL,
  drop_empty_rows = FALSE,
  drop_empty_cols = TRUE,
  binary = FALSE,
  check_trans = TRUE,
  na.rm = FALSE,
  na.replace = NULL
)
```

## Arguments

- com:

  community matrix (sites x species)

- D:

  Optional. Numeric matrix with two columns

- id:

  Optional. Unique id for each sample

- method:

  Method used to calculate community dissimilarity between pairs of
  sites. Either a distance metric supported by vegan, 'abcd' to obtain
  each component of the binary contingecy matrix, 'decomp1' or 'decomp2'
  for Pierre Legendre's or Andres Baselga's turnover and nestedness
  decomposition

- num_den:

  For some dissimilarity indices, whether function should return the
  numerator and denominator in two different columns

- trans:

  Transformation to be applied to raw data before calculating distance
  matrix. Either a function or 'binary' to transform to
  presence-absence)

- drop_empty_rows:

  Indicates whether empty rows (i.e., sites containing 0s only) should
  be removed

- drop_empty_cols:

  Indicates whether empty columns (i.e., species without observations)
  should be removed

- binary:

  Indicates whether abundance data should be transform into presence
  absence

- check_trans:

  Check for non-finite values in transformation

- na.rm:

  Remove rows with NA?

- na.replace:

  Replace NAs with specific value

## Value

A data frame with the following columns

- s1, s2 - A combination of sites (rows).

- Additional columns indicating the dissimilarity between sites.
