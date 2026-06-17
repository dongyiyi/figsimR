# Calculate emergent secondary community metrics

Calculates four emergent secondary metrics from a fig-level community
matrix: mean species richness per fig, mean Pielou's evenness per fig,
mean pairwise Bray-Curtis dissimilarity, and multivariate dispersion.

## Usage

``` r
calc_secondary_metrics(comm_mat)
```

## Arguments

- comm_mat:

  A matrix or data frame with figs as rows and wasp species as columns.
  Values should be species abundances.

## Value

A one-row data frame with columns `mean_richness`,
`mean_pielou_evenness`, `mean_bray_curtis`, and
`multivariate_dispersion`.
