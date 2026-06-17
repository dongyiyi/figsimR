# Resample and Compute Community Metrics from Simulated Data

This function performs repeated subsampling of simulated fig wasp data
to assess the variability in community-level metrics. It automatically
detects columns containing emergence counts (i.e., columns starting with
\`"emergence\_"\`) and calculates diversity or network statistics using
a user-defined \`calc_all_metrics()\` function.

## Usage

``` r
resample_metrics_from_simulation(
  sim_data,
  n_reps = 500,
  sample_size = 200,
  seed = 42
)
```

## Arguments

- sim_data:

  A data frame returned by the simulator, typically
  `sim_result$summary`.

- n_reps:

  Number of bootstrap replicates to perform. Default is 500.

- sample_size:

  Number of figs to draw per replicate. Default is 200.

- seed:

  Random seed for reproducibility. Default is 42.

## Value

A data frame with one row per replicate, each row containing calculated
community metrics.

## Details

Columns with names starting with `"emergence_"` are automatically
extracted as wasp emergence data. If the metric computation fails for a
given replicate, it will be skipped with a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
sim_result <- simulate_figwasp_community(num_figs = 1000, ...)
summary_df <- sim_result$summary
metrics_df <- resample_metrics_from_simulation(summary_df)
} # }
```
