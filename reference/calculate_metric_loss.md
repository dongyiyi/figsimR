# Calculate Total Loss Between Simulated and Observed Metrics

This function compares the mean values of community metrics from
simulated data to a single row of observed values, computing the squared
error loss for each metric and returning both the per-metric loss and
the total loss.

## Usage

``` r
calculate_metric_loss(simulated_df, observed_df)
```

## Arguments

- simulated_df:

  A data frame of simulated metric results, typically generated from
  bootstrap resampling (e.g., via
  [`resample_metrics_from_simulation()`](https://dongyiyi.github.io/figsimR/reference/resample_metrics_from_simulation.md)).

- observed_df:

  A one-row data frame containing observed metric values. Column names
  should match those in `simulated_df`.

## Value

A list with two elements:

- `metric_table`:

  A data frame showing the simulated mean, observed value, and squared
  difference for each metric.

- `total_loss`:

  The total loss, defined as the sum of squared differences across all
  shared metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
sim_metrics <- resample_metrics_from_simulation(sim_data)
loss_result <- calculate_metric_loss(sim_metrics, observed_result_df)
print(loss_result$total_loss)
print(loss_result$metric_table)
} # }
```
