# Bootstrap Resampling of Simulated Community Metrics

This function performs multiple independent runs of a fig wasp community
simulation model and computes community metrics for a fixed number of
sampled figs from each run. It is used to generate the expected
distribution of metrics under a given simulation configuration.

## Usage

``` r
resample_simulated_metrics(
  n_draws = 500,
  sample_n = 200,
  num_figs = 1000,
  simulator_func = simulate_figwasp_community,
  calc_func = calc_all_metrics,
  wasp_cols,
  ...
)
```

## Arguments

- n_draws:

  Integer. Number of simulation replicates to generate. Default = 500.

- sample_n:

  Integer. Number of figs to sample per replicate for metric
  calculation. Default = 200.

- num_figs:

  Integer. Number of figs to simulate per call to the simulator. Default
  = 1000.

- simulator_func:

  A function that simulates fig wasp communities and returns a list with
  at least a `summary` data.frame including wasp emergence and an
  `is_dropped` column. Default =
  [`simulate_figwasp_community`](https://dongyiyi.github.io/figsimR/reference/simulate_figwasp_community.md).

- calc_func:

  A function to calculate community metrics from a sampled data.frame.
  Must accept a data.frame of species abundance (rows = figs, columns =
  species). Default =
  [`calc_all_metrics`](https://dongyiyi.github.io/figsimR/reference/calc_all_metrics.md).

- wasp_cols:

  A character vector indicating which columns in the simulation output
  correspond to wasp species.

- ...:

  Additional arguments passed to `simulator_func`, such as parameter
  lists or switches.

## Value

A data.frame where each row represents a replicate and includes all
community metrics returned by `calc_func`. A column
`source = "Simulated"` is added to distinguish simulated outputs.

## Details

In each replicate, the simulation function is called with `num_figs`
figs. Only figs where `is_dropped == 0` are retained, and a sample of
size `sample_n` is randomly drawn for metric calculation. If a replicate
contains fewer than `sample_n` valid figs, that iteration is skipped.

This function is commonly used to generate the expected distribution of
metrics under the Full Intrinsic Model (FIM) or other simulation
variants, for comparison with observed data.

## Examples

``` r
if (FALSE) { # \dontrun{
sim_metrics <- resample_simulated_metrics(
  n_draws = 100,
  sample_n = 100,
  num_figs = 500,
  simulator_func = simulate_figwasp_community,
  calc_func = calc_all_metrics,
  wasp_cols = c("emergence_Ceratosolen_sp",
                "emergence_Sycophaga_testacea",
                "emergence_Apocrypta_sp"),
  fecundity_mean = parameter_list_default$fecundity_mean,
  fecundity_dispersion = parameter_list_default$fecundity_dispersion,
  entry_mu = parameter_list_default$entry_mu,
  entry_size = parameter_list_default$entry_size,
  entry_priority = parameter_list_default$entry_priority,
  species_roles = parameter_list_default$species_roles,
  max_entry_table = parameter_list_default$max_entry_table,
  egg_success_prob = parameter_list_default$egg_success_prob,
  egg_success_prob_by_phase = parameter_list_default$egg_success_prob_by_phase,
  layer_preference = parameter_list_default$layer_preference
)
head(sim_metrics)
} # }
```
