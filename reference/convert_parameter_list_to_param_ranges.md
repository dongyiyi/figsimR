# Convert a Parameter List into a Param Range List for Sampling

Converts a parameter list into a list of value ranges suitable for Latin
Hypercube Sampling or Bayesian optimization.

## Usage

``` r
convert_parameter_list_to_param_ranges(
  param_list,
  prop = 0.3,
  floor_list = list(default = 0.01),
  ceiling_list = list(default = Inf, egg_success_prob = 1, layer = 1, parasitism_prob =
    1)
)
```

## Arguments

- param_list:

  A named list of species parameters used in the fig wasp simulator.
  Must include at least: `entry_mu`, `entry_size`, `fecundity_mean`,
  `fecundity_dispersion`, and `egg_success_prob`. Optional entries
  include `layer_preference`, `parasitism_prob`, and
  `egg_success_prob_by_phase`.

- prop:

  A numeric scalar (e.g., 0.3) indicating the percentage range around
  each parameter's base value to use for lower and upper bounds.

- floor_list:

  A named list setting minimum values for each parameter type. Example
  values: `default = 0.01` `egg_success_prob = 0.01` `layer = 0.01`

- ceiling_list:

  A named list setting maximum values for each parameter type. Example
  values: `default = Inf`; `egg_success_prob = 1`; `layer = 1`;
  `parasitism_prob = 1`.

## Value

A named list of numeric vectors of length 2 (lower, upper bounds),
including: `entry_mu_*`, `entry_size_*`, `fecundity_mean_*`,
`egg_success_prob_*` If applicable: `layer_core_raw_*`,
`parasitism_prob_*`, and `interaction_weight`

## Details

This function prevents zero or NA values by assigning default values.
Bounds are adjusted based on proportional range `prop`, and clipped to
valid floor and ceiling values per parameter type.

## Examples

``` r
data(parameter_list_default)
param_input <- parameter_list_default
param_input[c("interaction_matrix", "interaction_weight")] <- NULL
ranges <- convert_parameter_list_to_param_ranges(param_input)
str(ranges[1:3])
#> List of 3
#>  $ entry_mu_Sycophaga_testacea      : num [1:2] 4.2 7.8
#>  $ entry_size_Sycophaga_testacea    : num [1:2] 6.3 11.7
#>  $ fecundity_mean_Sycophaga_testacea: num [1:2] 15.4 28.6
```
