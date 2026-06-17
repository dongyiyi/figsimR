# Build predictive envelopes for model metrics

Computes bootstrap predictive envelopes for selected metrics from one
model output.

## Usage

``` r
envelope_for_model(
  df_one_model,
  metric_cols,
  q_levels = c(0.05, 0.25, 0.5, 0.75, 0.95),
  B = 1000
)
```

## Arguments

- df_one_model:

  A data frame containing model metric values.

- metric_cols:

  Character vector of metric columns to include.

- q_levels:

  Numeric vector of quantile levels used to build envelopes.

- B:

  Integer. Number of bootstrap replicates.

## Value

A data frame with predictive-envelope summaries for each metric and
quantile.
