# Bootstrap Resampling of Observed Community Metrics

This function resamples the observed fig wasp community dataset multiple
times to estimate the distribution of community metrics, such as
diversity and network structure. It is typically used to compare
observed patterns against simulated outputs.

## Usage

``` r
resample_observed_all_metrics(
  observed_df,
  wasp_cols,
  n_draws = 500,
  sample_n = 200,
  seed = NULL,
  calc_func = calc_all_metrics
)
```

## Arguments

- observed_df:

  A data.frame containing fig-level community observations, with one row
  per fig and one column per wasp species (typically emergence counts).

- wasp_cols:

  A character vector indicating which columns in `observed_df`
  correspond to wasp species. These columns will be used to compute
  metrics.

- n_draws:

  Integer. Number of bootstrap replicates to generate. Default = 500.

- sample_n:

  Integer. Number of figs to sample per bootstrap iteration. Default =
  200.

- seed:

  Optional integer. If provided, sets the random seed for
  reproducibility.

- calc_func:

  A function to calculate community metrics. Must accept a data.frame of
  species abundance (rows = figs, columns = species). Default =
  [`calc_all_metrics`](https://dongyiyi.github.io/figsimR/reference/calc_all_metrics.md).

## Value

A data.frame containing `n_draws` rows, each row representing a
bootstrap replicate. Each row includes all metrics returned by
`calc_func`, along with a column `source = "Observed"`.

## Details

In each iteration, the function randomly samples `sample_n` rows (with
replacement) from the observed dataset, calculates community metrics
using `calc_func`, and stores the results.

This resampling approach estimates the empirical distribution of
community metrics under the observed data, providing a null reference
for comparison with simulation outputs.

## Examples

``` r
# Example with toy data
set.seed(123)
df <- data.frame(
  fig_id = 1:300,
  A = rpois(300, 5),
  B = rpois(300, 3),
  C = rpois(300, 1)
)

observed_metrics <- resample_observed_all_metrics(
  observed_df = df,
  wasp_cols = c("A", "B", "C"),
  n_draws = 100,
  sample_n = 50,
  seed = 42
)

head(observed_metrics)
#>   mean_richness mean_shannon mean_simpson mean_evenness mean_bray_curtis
#> 1          2.48    0.7670380    0.4879917     0.8812799        0.3182976
#> 2          2.62    0.8186713    0.5128465     0.8684918        0.3266789
#> 3          2.50    0.7768747    0.4896270     0.8669182        0.3284078
#> 4          2.58    0.8105405    0.5072708     0.8753333        0.2994410
#> 5          2.56    0.7811148    0.4859879     0.8438040        0.3155640
#> 6          2.58    0.7912086    0.4940710     0.8611107        0.3041828
#>   mean_jaccard nestedness connectance links_per_species modularity   source
#> 1    0.2253061   26.13252   0.8266667          23.05645 0.08243366 Observed
#> 2    0.1700680   22.14729   0.8733333          23.96947 0.06444846 Observed
#> 3    0.2201361   24.09408   0.8333333          23.11200 0.08355200 Observed
#> 4    0.1889796   29.72050   0.8600000          23.68992 0.06715342 Observed
#> 5    0.2068027   27.43282   0.8533333          23.35938 0.07638550 Observed
#> 6    0.2065306   27.43386   0.8600000          23.43411 0.06715342 Observed
```
