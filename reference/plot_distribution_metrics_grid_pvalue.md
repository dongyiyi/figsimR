# Plot Metric Distributions with Wilcoxon Test Significance

This function compares the distributions of community-level metrics
between observed and simulated datasets. It visualizes the density plots
for each numeric metric side-by-side and annotates each with the
Wilcoxon test p-value and corresponding significance level.

## Usage

``` r
plot_distribution_metrics_grid_pvalue(obs_metrics, sim_metrics, ncol = 3)
```

## Arguments

- obs_metrics:

  A data frame of observed community metrics (e.g., alpha, beta, network
  structure). Must include only numeric columns to be tested and
  plotted.

- sim_metrics:

  A data frame of simulated metrics from the same set of metrics as
  `obs_metrics`. The columns should have the same names and types as in
  `obs_metrics`.

- ncol:

  Integer. Number of columns to arrange plots in the output grid.
  Default is 3.

## Value

A grid of `ggplot2` density plots, one for each shared numeric metric.
Each plot includes a Wilcoxon rank-sum test result between the observed
and simulated values, reported as p-value and a significance label
(\*\*\* \< 0.001, \*\* \< 0.01, \* \< 0.05, ns otherwise).

## Details

This function is useful for model evaluation, allowing users to assess
whether the distributions of key diversity or network metrics produced
by a simulation (e.g., FIM) differ significantly from those observed in
empirical data. It uses nonparametric Wilcoxon rank-sum tests for
comparison.

Metrics that are not numeric or are not shared between `obs_metrics` and
`sim_metrics` will be ignored.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_distribution_metrics_grid_pvalue(observed_summary_df, simulated_summary_df, ncol = 4)
} # }
```
