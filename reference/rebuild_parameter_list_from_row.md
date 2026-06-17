# Rebuild Parameter List from a Row of Flat Parameter Values

This function reconstructs a nested `parameter_list` object by inserting
flat parameter values (e.g., from optimization or Latin Hypercube
Sampling results) into the appropriate fields of a template list.

## Usage

``` r
rebuild_parameter_list_from_row(row, parameter_list_template)
```

## Arguments

- row:

  A named vector or single-row data.frame containing flat parameter
  values. Names must follow specific naming patterns (see Details).

- parameter_list_template:

  A complete simulation parameter list (e.g., `parameter_list_default`)
  to be used as the base for reconstruction.

## Value

A fully structured `parameter_list` list, where values in `row` are
inserted into their corresponding nested positions.

## Details

This function supports multiple naming conventions for parameter fields:

- **`egg_success_prob_by_phase_Species_phaseX`**: assigns phase-specific
  egg success probability (e.g.,
  `egg_success_prob_by_phase_Apocrypta_sp_phase2 = 0.6`).

- **`layer_core_raw_Species`**, `layer_mid_raw_Species`,
  `layer_outer_raw_Species`: assigns spatial layer preferences for
  species (e.g., `layer_mid_raw_PollA = 0.3`).

- **`egg_success_prob_Species`**: assigns non-phase-specific egg success
  probability.

- **`prefix_species`**: general pattern for parameters such as
  `entry_mu`, `entry_size`, `fecundity_mean`, etc.

- **`prefix_subprefix_species`**: supports two-part prefixes for
  compatibility with list structures such as
  `entry_priority$submodel$species` (less common).

The function parses each parameter name, extracts the relevant field and
species identity, and updates the `parameter_list_template` accordingly.
Unrecognized fields are ignored.

## Examples

``` r
if (FALSE) { # \dontrun{
row <- list(
  entry_mu_PollA = 0.5,
  fecundity_mean_PollA = 12,
  egg_success_prob_by_phase_Apocrypta_sp_phase2 = 0.6,
  layer_mid_raw_PollA = 0.3
)
param_list <- rebuild_parameter_list_from_row(row, parameter_list_default)
} # }
```
