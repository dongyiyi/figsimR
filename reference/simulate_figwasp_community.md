# Simulate Fig Wasp Community Assembly (Observed on *Ficus racemosa*)

This function simulates fig wasp community structure within a set of fig
fruits, incorporating priority effects, host–parasitoid interactions,
resource constraints, pollination, oviposition success, and fig abortion
mechanisms.

## Usage

``` r
simulate_figwasp_community(
  num_figs = 1000,
  fig_diameter_mean = 2.5,
  fig_diameter_sd = 1.2,
  k = 300,
  alpha = 1.3,
  max_entry_table = NULL,
  fecundity_mean,
  fecundity_dispersion,
  entry_mu,
  entry_size,
  entry_priority,
  species_roles,
  entry_distribution = "lognormal",
  interaction_matrix = NULL,
  interaction_weight = 0,
  seed = NULL,
  egg_success_prob = NULL,
  egg_success_prob_by_phase = NULL,
  parasitism_prob = NULL,
  layer_preference = NULL,
  p_pollination_per_ovule = 0.98,
  p_no_entry = 0.002,
  host_sanction = 0.8,
  fig_diameter_min = 1.3,
  fig_diameter_max = 4,
  record_individual = FALSE,
  use_layering = TRUE,
  use_layer_preference = TRUE,
  enable_drop = TRUE,
  drop_cancels_emergence = FALSE,
  use_supplemental_parasitism = FALSE,
  use_flower_limit = TRUE,
  use_egg_success_by_phase = TRUE,
  use_sink_strength = FALSE,
  sink_w_gall = 1,
  sink_w_seed = 1.5,
  sink_linear_coef = 1,
  sink_min_prop = 0.2,
  sink_max_prop = 0.95
)
```

## Arguments

- num_figs:

  Integer. Number of fig fruits to simulate (i.e., number of independent
  communities, default = 1000). Each fig will be treated as a discrete
  community.

- fig_diameter_mean:

  Numeric. Mean fig diameter (in cm, default = 2.5). It is used to
  determine ovule number. Affects resource availability.

- fig_diameter_sd:

  Numeric. Standard deviation of fig diameter under the normal
  distribution used to generate fig diameters, default = 1.2).

- k:

  Numeric. Scaling constant for estimating flower number from fig
  diameter (default = 300). Controls ovule-based competition.

- alpha:

  Numeric. Exponent in the diameter–flower power-law used to estimate
  flower number from fig diameter; controls shape of seed-yield curve
  (default = 1.3).

- max_entry_table:

  Named numeric vector giving the baseline maximum number of arrivals or
  oviposition-attempting wasps per species. Names must match
  `species_roles$guild`. If `NULL`, the function uses the default *Ficus
  racemosa* values when the supplied species set matches the worked
  example; otherwise, a generic value is assigned to all species, with a
  warning recommending system-specific values.

- fecundity_mean:

  The mean potential number of eggs per individual wasp.

- fecundity_dispersion:

  Named numeric vectors for fecundity distribution. The dispersion
  parameter for the negative binomial distribution modeling individual
  fecundity.

- entry_mu:

  Named numeric vectors controlling the expected number of wasp arrival
  or oviposition-attempt events per fig. For pollinators, this
  corresponds to physical entry into the fig cavity; for non-pollinating
  fig wasps, it represents arrival at the fig and oviposition attempts
  from outside the fig wall. `entry_mu` should not be interpreted as
  physical entry into the fig cavity for all species.

- entry_size:

  Named numeric vectors controlling entry distribution (the dispersion
  of the arrival or oviposition-attempt distribution). See
  `entry_distribution`. For `entry_distribution = "nb"`: used as NB size
  parameter. For `lognormal`: used to derive
  `sdlog = 1 / sqrt(entry_size)`. In the negative binomial option, this
  corresponds to the "size" parameter; smaller values produce stronger
  aggregation and greater among-fig variation. In the lognormal option,
  it controls the log-scale variance of the arrival process. This
  parameter is not a sample size and does not represent the number of
  eggs.

- entry_priority:

  Named list of phases and species.

- species_roles:

  List with \$guild, \$hosts, \$parasitoid.

- entry_distribution:

  "nb" or "lognormal". See `entry_size`.

- interaction_matrix:

  Optional numeric matrix (pairwise interactions).

- interaction_weight:

  Numeric. Strength of interaction effects.

- seed:

  RNG seed (optional).

- egg_success_prob:

  Global probability that a single oviposition attempt results in a
  successful egg. Applies when `use_egg_success_by_phase = FALSE`.

- egg_success_prob_by_phase:

  Named list of species-specific, phase-specific egg-success
  probabilities. Each element is a named numeric vector whose names
  correspond to phases in `entry_priority`. When
  `use_egg_success_by_phase = TRUE`, a matching species-phase value
  overrides the baseline value in `egg_success_prob`. If no
  phase-specific value is provided, the baseline `egg_success_prob` is
  used; if that is also missing, the default success probability is 1.

- parasitism_prob:

  Named numeric vector for supplemental parasitism.

- layer_preference:

  Named list specifying species-specific ovary-layer preference or
  accessibility. Each species can be assigned either a numeric vector
  giving relative probabilities for `core`, `mid`, or `outer` ovary
  layers; Or `none` to indicate no species-specific ovary-layer
  preference. If `NULL`, figsimR uses the default *Ficus
  racemosa*-specific layer preference when the supplied species set
  matches the worked example; otherwise, species are assigned `none` and
  users are encouraged to provide system-specific values.

- p_pollination_per_ovule:

  Numeric \[0–1\]. Probability that each unoccupied ovule becomes a
  seed.

- p_no_entry:

  Numeric \[0–1\]. Probability a fig receives no entries.

- host_sanction:

  Numeric \[0–1\]. Legacy overuse threshold for drop logic.

- fig_diameter_min:

  Numeric. Minimum fig diameter used to cap values (optional).
  Truncation bounds (if used). Truncate simulated diameters (and thus
  flower counts) to a biologically reasonable range.

- fig_diameter_max:

  Numeric. Maximum fig diameter used to cap values (optional).
  Truncation bounds (if used). Truncate simulated diameters (and thus
  flower counts) to a biologically reasonable range.

- record_individual:

  Logical. If `TRUE`, per-individual oviposition recorded.

- use_layering:

  Logical switches for layering. If `TRUE`, fig flowers are divided into
  spatial layers (core, mid, outer).

- use_layer_preference:

  Logical. If `TRUE`, apply species-specific ovary-layer preference or
  accessibility values supplied through `layer_preference`. If `FALSE`,
  ovary-layer resources are still used when `use_layering = TRUE`, but
  no species-specific layer preference is applied.

- enable_drop:

  Logical. Enable legacy drop (host sanction) mechanism. If flower use
  proportion exceeds host_sanction, fig may abort. In the output
  "drop_prob" column: Ranges from 0 (never drop) to 1 (always drop);
  "is_dropped" column: no drop (0) and drop (1). When sink strength is
  used (`use_sink_strength = TRUE`) and `enable_drop = FALSE` (host
  sanction is disabled); the columns "is_dropped" and "drop_prob" in the
  outputs are not meaningful.

- drop_cancels_emergence:

  Logical. If `FALSE` (figs are dropped and no wasps emerge), dropped
  figs yield zero output (no wasps emerge at the end of phase).

- use_supplemental_parasitism:

  Logical. Whether to activate supplemental parasitism. Default FALSE.

- use_flower_limit:

  Logical. If TRUE, clamps fig diameter between `fig_diameter_min` and
  `fig_diameter_max` before calculating flower number.

- use_egg_success_by_phase:

  Logical. Use egg success probabilities by phase. Default TRUE.

- use_sink_strength:

  Logical. If `TRUE`, compute sink-strength metrics but do NOT alter
  legacy drop decision: `use_sink_strength = TRUE` and
  `enable_drop = FALSE` and the columns "is_dropped" and "drop_prob" are
  not meaningful. The fig is dropped when the value of "sink_prop" falls
  outside the thresholds defined by `sink_min_prop` and `sink_max_prop`;
  The columns "sink_below_min" and "sink_above_max" are logical: TRUE =
  drop; FALSE = no drop. Adds columns: `sink_strength`, `sink_prop`,
  `sink_below_min`, `sink_above_max`. Sink strength threshold
  calculation:
  `sink_strength = sink_linear_coef * (sink_w_gall * galled_total + sink_w_seed * seed_count)`.

- sink_w_gall:

  Numeric in \[0, +). Per-ovule sink weight for galled ovules (default
  1.0).

- sink_w_seed:

  Numeric in \[0, +). Per-ovule sink weight for seeds (default 1.5).

- sink_linear_coef:

  Numeric. Global multiplier on sink strength (default 1.0).

- sink_min_prop:

  Numeric in \[0,1\]. Reference lower bound of sink proportion (default
  0.20).

- sink_max_prop:

  Numeric in \[0,1\]. Reference upper bound of sink proportion (default
  0.95).

## Value

A list containing:

- `summary`: a compact data frame containing fig-level variables, total
  entry, egg, and emergence counts for each species, seed and ovule
  outcomes, sink-strength metrics, drop status, and the entry matrix.

- `original_summary`: the complete simulation output, including
  layer-specific resource-use and egg-count columns and other diagnostic
  variables.

- `individual_eggs`: a list of per-individual egg counts, returned only
  when `record_individual = TRUE`.

## Details

Each simulation models community assembly in a number of fig fruits,
where:

- Flower number is determined by fig diameter and stochastic
  heterogeneity.

- Wasp species arrive in `entry_priority` phases (priority effects).

- Pollinators and gallers lay eggs into limited ovules.

- Parasitoids require host presence and attack host eggs.

- The ovary-layer preference filter can be bypassed with `"none"`; Or
  layer preference can restrict oviposition to `core`, `mid`, or `outer`
  layers.

- Figs may abort (drop) if flower usage exceeds `host_sanction` or
  `sink_strength`.

Species roles, hosts, and parasitoid links must be defined in
`species_roles`. Default species and parameter values are provided in
[`species_list`](https://dongyiyi.github.io/figsimR/reference/species_list.md)
and
[`parameter_list_default`](https://dongyiyi.github.io/figsimR/reference/parameter_list_default.md).

## Examples

``` r
data(species_list)
data(parameter_list_default)
set.seed(123)
sim_result <- simulate_figwasp_community(
  num_figs = 10,
  fecundity_mean = parameter_list_default$fecundity_mean,
  fecundity_dispersion = parameter_list_default$fecundity_dispersion,
  entry_mu = parameter_list_default$entry_mu,
  entry_size = parameter_list_default$entry_size,
  entry_priority = parameter_list_default$entry_priority,
  species_roles = parameter_list_default$species_roles,
  max_entry_table = parameter_list_default$max_entry_table,
  interaction_matrix = parameter_list_default$interaction_matrix,
  interaction_weight = parameter_list_default$interaction_weight,
  egg_success_prob = parameter_list_default$egg_success_prob,
  egg_success_prob_by_phase = parameter_list_default$egg_success_prob_by_phase,
  parasitism_prob = parameter_list_default$parasitism_prob,
  layer_preference = parameter_list_default$layer_preference,
  record_individual = TRUE,
  use_sink_strength = TRUE,
  sink_w_gall = 1.0,
  sink_w_seed = 1.5,
  sink_linear_coef = 1.0,
  sink_min_prop = 0.20,
  sink_max_prop = 0.95
)
head(sim_result$summary)
#>   fig_id fig_diameter flower_count resource_use richness_skipped
#> 1      1     1.827429         1661          319                0
#> 2      2     2.223787         1791          590                0
#> 3      3     4.000000         3876         1361                0
#> 4      4     2.584610         2061          746                0
#> 5      5     2.655145         1828          859                0
#> 6      6     4.000000         4663          253                0
#>   entry_Ceratosolen_sp eggs_Ceratosolen_sp entry_Sycophaga_mayri
#> 1                    2                 243                     2
#> 2                    6                 350                     8
#> 3                   16                1062                    15
#> 4                    9                 563                     4
#> 5                   10                 674                     6
#> 6                    2                  77                     8
#>   eggs_Sycophaga_mayri entry_Sycophaga_testacea eggs_Sycophaga_testacea
#> 1                   29                        4                      47
#> 2                  240                        1                       0
#> 3                  195                        5                     104
#> 4                  114                        5                      69
#> 5                  158                        2                      27
#> 6                  152                        2                      24
#>   entry_Apocrypta_sp eggs_Apocrypta_sp entry_Apocrypta_westwoodi
#> 1                  6                 9                         8
#> 2                 11                 0                        10
#> 3                 13                 9                         9
#> 4                 23                47                         8
#> 5                 21                19                         7
#> 6                  7                20                         8
#>   eggs_Apocrypta_westwoodi entry_Sycophaga_agraensis eggs_Sycophaga_agraensis
#> 1                       21                         4                       10
#> 2                        2                         5                        7
#> 3                       26                         7                        2
#> 4                        2                        10                       10
#> 5                       14                         9                       17
#> 6                       13                         9                       14
#>   emergence_Ceratosolen_sp emergence_Sycophaga_mayri
#> 1                      233                         8
#> 2                      343                       238
#> 3                     1060                       169
#> 4                      553                       112
#> 5                      657                       144
#> 6                       63                       139
#>   emergence_Sycophaga_testacea emergence_Apocrypta_sp
#> 1                           38                      9
#> 2                            0                      0
#> 3                           95                      9
#> 4                           22                     47
#> 5                            8                     19
#> 6                            4                     20
#>   emergence_Apocrypta_westwoodi emergence_Sycophaga_agraensis
#> 1                            21                            10
#> 2                             2                             7
#> 3                            26                             2
#> 4                             2                            10
#> 5                            14                            17
#> 6                            13                            14
#>   unoccupied_flowers seeds failed_ovules used_flowers resource_ratio
#> 1               1342  1318            24         1661     0.19205298
#> 2               1201  1175            26         1791     0.32942490
#> 3               2515  2462            53         3876     0.35113519
#> 4               1315  1291            24         2061     0.36196021
#> 5                969   951            18         1828     0.46991247
#> 6               4410  4325            85         4663     0.05425692
#>   sink_strength sink_prop sink_below_min sink_above_max    drop_prob is_dropped
#> 1        2296.0 0.9215332          FALSE          FALSE 0.0022841599          0
#> 2        2352.5 0.8756747          FALSE          FALSE 0.0089620754          0
#> 3        5054.0 0.8692810          FALSE          FALSE 0.0111109823          0
#> 4        2682.5 0.8677018          FALSE          FALSE 0.0123655548          0
#> 5        2285.5 0.8335157          FALSE          FALSE 0.0355411744          0
#> 6        6740.5 0.9636858          FALSE           TRUE 0.0005768041          0
#>          entry_matrix
#> 1    2, 2, 4, 6, 8, 4
#> 2  6, 8, 1, 11, 10, 5
#> 3 16, 15, 5, 13, 9, 7
#> 4  9, 4, 5, 23, 8, 10
#> 5  10, 6, 2, 21, 7, 9
#> 6    2, 8, 2, 7, 8, 9
```
