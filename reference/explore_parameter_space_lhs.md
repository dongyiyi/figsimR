# Explore Parameter Space Using Latin Hypercube Sampling

This function systematically explores the parameter space of the fig
wasp community simulator using Latin Hypercube Sampling (LHS), evaluates
simulation output against observed metrics, and returns the top
parameter combinations that minimize loss.

## Usage

``` r
explore_parameter_space_lhs(
  new_param_ranges,
  observed_summary,
  num_figs = 1000,
  n_draws = 500,
  sample_n = 200,
  wasp_cols,
  parameter_list_template,
  n_samples = 5000,
  top_k = 5,
  n_cores = 4
)
```

## Arguments

- new_param_ranges:

  A named list of parameter ranges, where each element is a numeric
  vector of length 2 specifying the lower and upper bounds for a
  parameter. Names should follow the format: `type_species` (e.g.,
  `entry_mu_PollA`) or `species_phaseX` for phase-specific values (e.g.,
  `Apocrypta_sp_phase1`).

- observed_summary:

  A data.frame of observed community metrics (e.g., richness, Shannon,
  Bray-Curtis), typically from
  [`resample_observed_all_metrics`](https://dongyiyi.github.io/figsimR/reference/resample_observed_all_metrics.md).

- num_figs:

  Integer. Number of figs to simulate in each replicate. Default = 1000.

- n_draws:

  Integer. Number of bootstrap replicates per parameter combination.
  Default = 500.

- sample_n:

  Integer. Number of figs sampled per replicate for calculating
  community metrics. Default = 200.

- wasp_cols:

  Character vector. Names of columns in simulation output corresponding
  to fig wasp species.

- parameter_list_template:

  A full template of simulation parameters (typically based on
  `parameter_list_default`), used to insert sampled parameter values.

- n_samples:

  Integer. Number of LHS parameter combinations to evaluate. Default =
  5000.

- top_k:

  Integer. Number of top-performing parameter sets (with lowest loss) to
  return. Default = 5.

- n_cores:

  Integer. Number of parallel cores to use. Default = 4.

## Value

A data.frame of the top `top_k` parameter combinations, sorted by loss,
with one row per combination and columns for each parameter and
corresponding loss value.

## Details

This function generates `n_samples` unique parameter sets using LHS
within the user-defined bounds, runs a stochastic simulation with each
set, and computes a loss value based on the squared difference between
simulated and observed mean metrics. Simulation is parallelized for
efficiency.

The `parameter_list_template` must include all required simulator
inputs. Sampled values are inserted by name-matching into their
corresponding list elements. Phase-specific parameters (e.g., egg
success probabilities by phase) must be named accordingly.

## Examples

``` r
if (FALSE) { # \dontrun{
best_params <- explore_parameter_space_lhs(
  new_param_ranges = list(
    fecundity_mean_PollA = c(5, 20),
    entry_mu_PollA = c(0.1, 2.0),
    Apocrypta_sp_phase1 = c(0.05, 0.3)
  ),
  observed_summary = resample_observed_all_metrics(...),
  num_figs = 1000,
  n_draws = 100,
  sample_n = 200,
  wasp_cols = c("PollA", "GallerA", "ParaA"),
  parameter_list_template = parameter_list_default,
  n_samples = 1000,
  top_k = 3,
  n_cores = 4
)
} # }
```
