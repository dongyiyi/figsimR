# Plot Metric Distributions with KS Test and Shared Legend

This function compares the distributions of community-level metrics
between observed and simulated datasets using kernel density plots. Each
metric is tested with a Kolmogorov-Smirnov (KS) test, and the resulting
p-value is displayed as the title of each plot. A shared legend is
displayed above all subplots.

## Usage

``` r
plot_distribution_metrics_grid_pvalue_legend(
  obs_metrics,
  sim_metrics,
  ncol = 3
)
```

## Arguments

- obs_metrics:

  A data frame containing observed metric values. Only numeric columns
  will be used.

- sim_metrics:

  A data frame containing simulated metric values. Must have the same
  structure as `obs_metrics`.

- ncol:

  Number of columns in the output plot grid. Default is 3.

## Value

A grid of density plots (one per metric) comparing observed and
simulated distributions, with a shared legend on top. Each plot shows a
KS test p-value between observed and simulated distributions.

## Details

This function helps evaluate model performance by comparing the
distribution of each community-level metric (e.g., alpha diversity,
richness, network structure) between observed and simulated datasets. A
Kolmogorov-Smirnov test is performed for each metric, and its p-value is
shown in the plot title.

The shared legend distinguishes between observed and simulated data and
is extracted from the first plot. All plots are styled consistently for
integration into figure panels or reports.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_distribution_metrics_grid_pvalue_legend(observed_df, simulated_df, ncol = 4)
} # }
```
