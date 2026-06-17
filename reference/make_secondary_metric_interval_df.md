# Create interval summaries for secondary metrics

Creates simulated predictive intervals and standardized observed
deviations for secondary community metrics.

## Usage

``` r
make_secondary_metric_interval_df(sim_df, obs_df)
```

## Arguments

- sim_df:

  A data frame of simulated secondary metrics, usually generated from
  repeated simulations or bootstrap resampling.

- obs_df:

  A one-row data frame of observed secondary metrics.

## Value

A data frame containing simulated means, standard deviations, 80 and
interval-status labels.
