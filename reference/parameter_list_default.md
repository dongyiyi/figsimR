# Default Parameter List for Simulating Fig Wasp Communities

This named list contains all required input parameters for running
[`simulate_figwasp_community`](https://dongyiyi.github.io/figsimR/reference/simulate_figwasp_community.md)
using the six fig wasp species observed on *Ficus racemosa*. It includes
species-specific trait values, community entry priority, host–parasitoid
relationships, and ecological interaction settings.

## Usage

``` r
data(parameter_list_default)
```

## Format

A named list with the following components:

- `entry_mu`:

  Named numeric vector of mean number of individuals of each species
  entering a fig.

- `entry_size`:

  Named numeric vector controlling entry variation; used as size in NB
  or sdlog in lognormal.

- `fecundity_mean`:

  Named numeric vector of the mean number of eggs per wasp (i.e., per
  female individual).

- `fecundity_dispersion`:

  Named numeric vector of fecundity dispersion (size parameter in
  negative binomial).

- `egg_success_prob`:

  Named numeric vector of baseline oviposition success probabilities
  (0–1).

- `egg_success_prob_by_phase`:

  Nested list of success probabilities for species-by-phase
  combinations. Overrides `egg_success_prob` if enabled.

- `layer_preference`:

  Named list giving relative probabilities of oviposition in fig layers:
  `core`, `mid`, `outer`. Sum should be 1. Used when
  `use_layering = TRUE`.

- `max_entry_table`:

  Named numeric vector of maximum number of individuals of each species
  allowed to enter a fig, adjusted by fig size.

- `parasitism_prob`:

  Named numeric vector of supplemental parasitism probabilities (0–1),
  used only if `use_supplemental_parasitism = TRUE`.

- `entry_priority`:

  Named list of phases. Each element is a character vector of species
  names representing entry order. Entry is processed by phase.

- `interaction_matrix`:

  Square matrix of interaction coefficients among species. Used to
  simulate facilitation or inhibition. Default is all zero.

- `interaction_weight`:

  Scalar controlling strength of interaction effects. Higher values
  amplify the influence of `interaction_matrix`.

- `species_roles`:

  A list defining ecological roles and trophic relationships. Contains:

  `guild`

  :   Named character vector assigning each species to a guild:
      `"pollinator"`, `"galler"`, or `"parasitoid"`.

  `hosts`

  :   Named list of host species for each parasitoid. Empty character
      vector if no host.

  `parasitoid`

  :   Named list of parasitoids attacking each host species. Empty if
      species is not parasitized.

## Details

All vectors and lists are named using standardized species names from
[`species_list`](https://dongyiyi.github.io/figsimR/reference/species_list.md).

## See also

[`simulate_figwasp_community`](https://dongyiyi.github.io/figsimR/reference/simulate_figwasp_community.md),
[`species_list`](https://dongyiyi.github.io/figsimR/reference/species_list.md)

## Examples

``` r
data(parameter_list_default)
names(parameter_list_default)
#>  [1] "entry_mu"                  "entry_size"               
#>  [3] "fecundity_mean"            "fecundity_dispersion"     
#>  [5] "egg_success_prob"          "egg_success_prob_by_phase"
#>  [7] "layer_preference"          "max_entry_table"          
#>  [9] "parasitism_prob"           "entry_priority"           
#> [11] "interaction_matrix"        "interaction_weight"       
#> [13] "species_roles"            
parameter_list_default$entry_priority
#> $phase1
#> [1] "Sycophaga_testacea" "Apocrypta_sp"      
#> 
#> $phase2
#> [1] "Apocrypta_sp"        "Sycophaga_mayri"     "Ceratosolen_sp"     
#> [4] "Sycophaga_agraensis" "Apocrypta_westwoodi"
#> 
#> $phase3
#> [1] "Sycophaga_agraensis" "Apocrypta_westwoodi"
#> 
parameter_list_default$species_roles$guild
#>      Ceratosolen_sp     Sycophaga_mayri  Sycophaga_testacea        Apocrypta_sp 
#>        "pollinator"            "galler"            "galler"        "parasitoid" 
#> Apocrypta_westwoodi Sycophaga_agraensis 
#>        "parasitoid"        "parasitoid" 
```
