# Build a bootstrap predictive envelope

Computes bootstrap summaries for selected quantiles of a numeric vector.

## Usage

``` r
predictive_envelope_full(
  v,
  q_levels = c(0.05, 0.25, 0.5, 0.75, 0.95),
  B = 1000
)
```

## Arguments

- v:

  Numeric vector of metric values.

- q_levels:

  Numeric vector of quantile levels to evaluate.

- B:

  Integer. Number of bootstrap replicates.

## Value

A data frame with columns `q`, `mean`, `sd`, `lo95`, and `hi95`.
