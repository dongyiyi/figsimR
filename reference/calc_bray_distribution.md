# Calculate pairwise Bray-Curtis dissimilarity distribution

Calculates pairwise Bray-Curtis dissimilarities among figs from a
community matrix.

## Usage

``` r
calc_bray_distribution(comm_mat, dataset = "Observed", max_pairs = 20000)
```

## Arguments

- comm_mat:

  A matrix or data frame with figs as rows and wasp species as columns.

- dataset:

  Character. Dataset label.

- max_pairs:

  Integer. Maximum number of pairwise distances to retain. If the number
  of pairwise distances exceeds this value, distances are randomly
  subsampled.

## Value

A data frame with columns `dataset` and `bray_curtis`.
