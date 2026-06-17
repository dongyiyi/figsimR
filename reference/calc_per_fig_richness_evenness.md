# Calculate fig-level richness and Pielou's evenness

Calculates species richness and Pielou's evenness for each fig.

## Usage

``` r
calc_per_fig_richness_evenness(comm_mat, dataset = "Observed")
```

## Arguments

- comm_mat:

  A matrix or data frame with figs as rows and wasp species as columns.

- dataset:

  Character. Dataset label, such as `"Observed"` or `"Simulated"`.

## Value

A data frame with columns `dataset`, `fig_id`, `richness`, and
`pielou_evenness`.
