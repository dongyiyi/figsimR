# Run Simulations and Bootstrap Community Metrics Across Interaction Weights

This function automates the process of running the fig wasp community
simulator with multiple values of `interaction_weight`, and performs
bootstrap resampling on each result to assess variability in community
metrics.

## Usage

``` r
run_sim_and_bootstrap_multiple_weights(
  weights = c(0.1, 0.3, 0.5, 0.7, 0.9),
  n_reps = 500,
  sample_size = 200,
  seed = 42
)
```

## Arguments

- weights:

  A numeric vector of `interaction_weight` values to evaluate. Default:
  c(0.1, 0.3, 0.5, 0.7, 0.9).

- n_reps:

  Number of bootstrap replicates per weight. Default: 500.

- sample_size:

  Number of figs to sample per replicate. Default: 200.

- seed:

  Random seed for reproducibility. Default: 42.

## Value

A named list of data frames. Each element corresponds to one interaction
weight and contains bootstrap resampled community metrics computed from
the simulation.

## Details

This function requires that global objects such as
`lhs_best_param_list`, `observed_interaction_matrix`, and
[`simulate_figwasp_community()`](https://dongyiyi.github.io/figsimR/reference/simulate_figwasp_community.md)
be pre-defined. The function uses
[`resample_metrics_from_simulation()`](https://dongyiyi.github.io/figsimR/reference/resample_metrics_from_simulation.md)
internally to compute alpha diversity or other community structure
metrics for each simulation run.

## Examples

``` r
if (FALSE) { # \dontrun{
results_list <- run_sim_and_bootstrap_multiple_weights(
  weights = seq(0.1, 0.9, by = 0.2),
  n_reps = 300,
  sample_size = 150
)
} # }
```
