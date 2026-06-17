# Bootstrap a quantile

Estimates a quantile by bootstrap resampling a numeric vector.

## Usage

``` r
boot_quantile(v, q, B = 1000)
```

## Arguments

- v:

  A numeric vector.

- q:

  Numeric. Quantile to estimate, between 0 and 1.

- B:

  Integer. Number of bootstrap replicates.

## Value

A numeric vector of bootstrap quantile estimates.
