# Calculate species richness per fig

Calculates the number of species present in each fig from a community
matrix.

## Usage

``` r
calc_richness_per_fig(comm_mat, dataset_name)
```

## Arguments

- comm_mat:

  A matrix or data frame with figs as rows and wasp species as columns.

- dataset_name:

  Character. Label for the dataset, such as \`"Observed"\` or
  \`"Simulated"\`.

## Value

A data frame with columns \`dataset\`, \`fig_id\`, and \`richness\`.
