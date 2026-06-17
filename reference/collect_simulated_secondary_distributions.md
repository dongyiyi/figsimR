# Collect simulated distributions of secondary metrics

Runs repeated simulations and collects fig-level richness/evenness,
pairwise Bray-Curtis dissimilarities, and distance-to-centroid values
from sampled simulated figs.

## Usage

``` r
collect_simulated_secondary_distributions(
  n_draws = 500,
  sample_n = 200,
  num_figs = 1000,
  wasp_cols,
  species_names_out,
  simulator_func = simulate_figwasp_community,
  max_bray_pairs_per_draw = 5000,
  ...
)
```

## Arguments

- n_draws:

  Integer. Number of repeated simulation draws.

- sample_n:

  Integer. Number of retained figs sampled from each simulation.

- num_figs:

  Integer. Number of figs simulated in each simulation run.

- wasp_cols:

  Character vector of column names corresponding to wasp emergence
  abundances in the simulator output.

- species_names_out:

  Character vector of species names to assign to the sampled community
  matrix.

- simulator_func:

  Simulation function. Default is `simulate_figwasp_community`.

- max_bray_pairs_per_draw:

  Integer. Maximum number of Bray-Curtis pairwise distances retained per
  draw.

- ...:

  Additional arguments passed to `simulator_func`.

## Value

A list with three data frames: `per_fig`, `bray`, and `dispersion`.
