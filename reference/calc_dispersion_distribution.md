# Calculate distance-to-centroid distribution

Calculates distances from each fig-level community to the multivariate
centroid using Bray-Curtis dissimilarity and
[`vegan::betadisper()`](https://vegandevs.github.io/vegan/reference/betadisper.html).

## Usage

``` r
calc_dispersion_distribution(comm_mat, dataset = "Observed")
```

## Arguments

- comm_mat:

  A matrix or data frame with figs as rows and wasp species as columns.

- dataset:

  Character. Dataset label.

## Value

A data frame with columns `dataset`, `fig_id`, and
`distance_to_centroid`.
