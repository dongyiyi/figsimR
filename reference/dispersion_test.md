# Compare multivariate dispersion between observed and simulated communities

Compares observed and simulated community matrices using multivariate
dispersion based on
[`vegan::betadisper()`](https://vegandevs.github.io/vegan/reference/betadisper.html).

## Usage

``` r
dispersion_test(
  obs_comm_mat,
  sim_comm_mat,
  transform = c("hellinger", "none"),
  permutations = 999,
  seed = 1
)
```

## Arguments

- obs_comm_mat:

  Observed community matrix with figs as rows and species as columns.

- sim_comm_mat:

  Simulated community matrix with figs as rows and species as columns.

- transform:

  Character. Transformation applied before distance calculation. Options
  are `"hellinger"` and `"none"`.

- permutations:

  Integer. Number of permutations used in
  [`vegan::permutest()`](https://vegandevs.github.io/vegan/reference/anova.cca.html).

- seed:

  Integer. Random seed for the permutation test.

## Value

A list containing the observed-to-simulated dispersion ratio `rho`,
permutation-test `p` value, group means, the `betadisper` object,
permutation result, distances, and group labels.
