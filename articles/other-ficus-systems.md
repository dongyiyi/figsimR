# Applying figsimR to other Ficus systems

## Is figsimR restricted to Ficus racemosa?

No. The default worked example is parameterized for `Ficus racemosa`,
but the number of simulated species is not hard-coded.

## What users need to provide

To adapt `figsimR` to another system, users need:

1.  Species names
2.  Guild assignments
3.  Host–parasitoid links
4.  Entry or oviposition schedules
5.  Fecundity parameters
6.  Layer preference or `"none"`
7.  Optional fig-retention rules

## Two cases: simulation for other “fig” systems

To demonstrate this flexibility, here, we tested the package using two
additional user-defined examples: a minimal two-species community and a
ten-species community that includes multiple pollinating fig wasp
species. In both cases, the species names and interaction structures
were user-defined. These examples show that figsimR can accommodate both
simplified and more complex fig-fig wasp communities.

### function for parameter check

Let define a parameter check function first

``` r

validate_figsim_parameter_list <- function(parameter_list) {

  required_components <- c(
    "entry_mu",
    "entry_size",
    "fecundity_mean",
    "fecundity_dispersion",
    "egg_success_prob",
    "egg_success_prob_by_phase",
    "layer_preference",
    "max_entry_table",
    "parasitism_prob",
    "entry_priority",
    "interaction_matrix",
    "interaction_weight",
    "species_roles"
  )

  missing_components <- setdiff(required_components, names(parameter_list))
  if (length(missing_components) > 0) {
    stop(
      "Missing required components: ",
      paste(missing_components, collapse = ", ")
    )
  }

  species_roles <- parameter_list$species_roles
  species_names <- names(species_roles$guild)

  if (length(species_names) == 0) {
    stop("species_roles$guild must be a named vector.")
  }

  # Check species-indexed numeric vectors
  species_vectors <- c(
    "entry_mu",
    "entry_size",
    "fecundity_mean",
    "fecundity_dispersion",
    "egg_success_prob",
    "max_entry_table"
  )

  for (v in species_vectors) {
    missing_species <- setdiff(species_names, names(parameter_list[[v]]))
    extra_species <- setdiff(names(parameter_list[[v]]), species_names)

    if (length(missing_species) > 0) {
      stop(v, " is missing species: ", paste(missing_species, collapse = ", "))
    }

    if (length(extra_species) > 0) {
      warning(v, " contains extra species: ", paste(extra_species, collapse = ", "))
    }
  }

  # Check hosts and parasitoid lists
  if (!all(species_names %in% names(species_roles$hosts))) {
    stop("species_roles$hosts must contain all species.")
  }

  if (!all(species_names %in% names(species_roles$parasitoid))) {
    stop("species_roles$parasitoid must contain all species.")
  }

  all_hosts <- unique(unlist(species_roles$hosts))
  all_parasitoids <- unique(unlist(species_roles$parasitoid))

  invalid_hosts <- setdiff(all_hosts, species_names)
  invalid_parasitoids <- setdiff(all_parasitoids, species_names)

  if (length(invalid_hosts) > 0) {
    stop("Unknown host species: ", paste(invalid_hosts, collapse = ", "))
  }

  if (length(invalid_parasitoids) > 0) {
    stop("Unknown parasitoid species: ", paste(invalid_parasitoids, collapse = ", "))
  }

  # Check entry priority
  entry_species <- unique(unlist(parameter_list$entry_priority))
  missing_from_priority <- setdiff(species_names, entry_species)
  invalid_priority_species <- setdiff(entry_species, species_names)

  if (length(missing_from_priority) > 0) {
    warning(
      "These species are not included in entry_priority: ",
      paste(missing_from_priority, collapse = ", ")
    )
  }

  if (length(invalid_priority_species) > 0) {
    stop(
      "entry_priority contains unknown species: ",
      paste(invalid_priority_species, collapse = ", ")
    )
  }

  # Check interaction matrix
  interaction_matrix <- parameter_list$interaction_matrix

  if (!is.matrix(interaction_matrix)) {
    stop("interaction_matrix must be a matrix.")
  }

  if (!all(species_names %in% rownames(interaction_matrix)) ||
      !all(species_names %in% colnames(interaction_matrix))) {
    stop("interaction_matrix rownames and colnames must include all species.")
  }

  # Check layer preference
  if (!all(species_names %in% names(parameter_list$layer_preference))) {
    stop("layer_preference must contain all species.")
  }

  # Check egg_success_prob_by_phase
  if (!all(species_names %in% names(parameter_list$egg_success_prob_by_phase))) {
    stop("egg_success_prob_by_phase must contain all species.")
  }

  message("Parameter list passed all checks.")
  invisible(TRUE)
}
```

### ten-species community with two pollinators

``` r

library(figsimR)
```

    ## 
    ## Attaching package: 'figsimR'

    ## The following object is masked from 'package:base':
    ## 
    ##     %||%

``` r

set.seed(42)

# define parameter list
custom_species <- c(
  "Pollinator_A",
  "Pollinator_B",
  "Galler_A",
  "Galler_B",
  "Galler_C",
  "Parasitoid_A",
  "Parasitoid_B",
  "Parasitoid_C",
  "Hyperparasitoid_A",
  "Hyperparasitoid_B"
)

custom_roles <- list(
  guild = c(
    Pollinator_A      = "pollinator",
    Pollinator_B      = "pollinator",
    Galler_A          = "galler",
    Galler_B          = "galler",
    Galler_C          = "galler",
    Parasitoid_A      = "parasitoid",
    Parasitoid_B      = "parasitoid",
    Parasitoid_C      = "parasitoid",
    Hyperparasitoid_A = "parasitoid",
    Hyperparasitoid_B = "parasitoid"
  ),

  hosts = list(
    Pollinator_A      = character(0),
    Pollinator_B      = character(0),
    Galler_A          = character(0),
    Galler_B          = character(0),
    Galler_C          = character(0),
    Parasitoid_A      = c("Galler_A"),
    Parasitoid_B      = c("Galler_B", "Galler_C"),
    Parasitoid_C      = c("Pollinator_A"),
    Hyperparasitoid_A = c("Parasitoid_A"),
    Hyperparasitoid_B = c("Parasitoid_B", "Parasitoid_C")
  ),

  parasitoid = list(
    Pollinator_A      = c("Parasitoid_C"),
    Pollinator_B      = character(0),
    Galler_A          = c("Parasitoid_A"),
    Galler_B          = c("Parasitoid_B"),
    Galler_C          = c("Parasitoid_B"),
    Parasitoid_A      = c("Hyperparasitoid_A"),
    Parasitoid_B      = c("Hyperparasitoid_B"),
    Parasitoid_C      = c("Hyperparasitoid_B"),
    Hyperparasitoid_A = character(0),
    Hyperparasitoid_B = character(0)
  )
)

custom_parameter_list <- list(
  entry_mu = c(
    Pollinator_A      = 8,
    Pollinator_B      = 4,
    Galler_A          = 5,
    Galler_B          = 4,
    Galler_C          = 4,
    Parasitoid_A      = 3,
    Parasitoid_B      = 3,
    Parasitoid_C      = 2,
    Hyperparasitoid_A = 1.5,
    Hyperparasitoid_B = 1.5
  ),

  entry_size = c(
    Pollinator_A      = 6,
    Pollinator_B      = 5,
    Galler_A          = 5,
    Galler_B          = 5,
    Galler_C          = 5,
    Parasitoid_A      = 4,
    Parasitoid_B      = 4,
    Parasitoid_C      = 4,
    Hyperparasitoid_A = 3,
    Hyperparasitoid_B = 3
  ),

  fecundity_mean = c(
    Pollinator_A      = 80,
    Pollinator_B      = 70,
    Galler_A          = 25,
    Galler_B          = 22,
    Galler_C          = 20,
    Parasitoid_A      = 8,
    Parasitoid_B      = 8,
    Parasitoid_C      = 7,
    Hyperparasitoid_A = 5,
    Hyperparasitoid_B = 5
  ),

  fecundity_dispersion = c(
    Pollinator_A      = 1.2,
    Pollinator_B      = 1.2,
    Galler_A          = 1.5,
    Galler_B          = 1.5,
    Galler_C          = 1.5,
    Parasitoid_A      = 2,
    Parasitoid_B      = 2,
    Parasitoid_C      = 2,
    Hyperparasitoid_A = 2,
    Hyperparasitoid_B = 2
  ),

  egg_success_prob = c(
    Pollinator_A      = 0.9,
    Pollinator_B      = 0.85,
    Galler_A          = 0.75,
    Galler_B          = 0.7,
    Galler_C          = 0.7,
    Parasitoid_A      = 0.6,
    Parasitoid_B      = 0.6,
    Parasitoid_C      = 0.55,
    Hyperparasitoid_A = 0.5,
    Hyperparasitoid_B = 0.5
  ),

  egg_success_prob_by_phase = list(
    Pollinator_A      = c(phase2 = 0.9),
    Pollinator_B      = c(phase2 = 0.8),
    Galler_A          = c(phase1 = 0.7),
    Galler_B          = c(phase1 = 0.7, phase2 = 0.6),
    Galler_C          = c(phase2 = 0.7),
    Parasitoid_A      = c(phase2 = 0.6),
    Parasitoid_B      = c(phase2 = 0.6, phase3 = 0.6),
    Parasitoid_C      = c(phase2 = 0.5),
    Hyperparasitoid_A = c(phase3 = 0.5),
    Hyperparasitoid_B = c(phase3 = 0.5)
  ),

  layer_preference = list(
    Pollinator_A      = c(core = 0.7, mid = 0.2, outer = 0.1),
    Pollinator_B      = c(core = 0.6, mid = 0.3, outer = 0.1),
    Galler_A          = c(core = 0.5, mid = 0.3, outer = 0.2),
    Galler_B          = c(core = 0.4, mid = 0.4, outer = 0.2),
    Galler_C          = c(core = 0.3, mid = 0.4, outer = 0.3),
    Parasitoid_A      = c(core = 0.2, mid = 0.4, outer = 0.4),
    Parasitoid_B      = c(core = 0.2, mid = 0.3, outer = 0.5),
    Parasitoid_C      = c(core = 0.3, mid = 0.3, outer = 0.4),
    Hyperparasitoid_A = c(core = 0.1, mid = 0.3, outer = 0.6),
    Hyperparasitoid_B = c(core = 0.1, mid = 0.3, outer = 0.6)
  ),

  max_entry_table = setNames(rep(20, length(custom_species)), custom_species),

  parasitism_prob = c(
    Parasitoid_A      = 0.85,
    Parasitoid_B      = 0.80,
    Parasitoid_C      = 0.75,
    Hyperparasitoid_A = 0.60,
    Hyperparasitoid_B = 0.60
  ),

  entry_priority = list(
    phase1 = c("Galler_A", "Galler_B"),
    phase2 = c(
      "Pollinator_A",
      "Pollinator_B",
      "Galler_B",
      "Galler_C",
      "Parasitoid_A",
      "Parasitoid_B",
      "Parasitoid_C"
    ),
    phase3 = c(
      "Parasitoid_B",
      "Hyperparasitoid_A",
      "Hyperparasitoid_B"
    )
  ),

  interaction_matrix = matrix(
    0,
    nrow = length(custom_species),
    ncol = length(custom_species),
    dimnames = list(custom_species, custom_species)
  ),

  interaction_weight = 0,

  species_roles = custom_roles
)
```

### check parameter for the community with ten wasp species

``` r

validate_figsim_parameter_list(custom_parameter_list)
```

    ## Parameter list passed all checks.

Everthing looks great, let’s simulate the community!

``` r

set.seed(123)

custom_sim <- do.call(
  simulate_figwasp_community,
  c(
    list(
      num_figs = 50,
      seed = 123,
      use_layering = TRUE,
      use_layer_preference = TRUE,
      use_sink_strength = TRUE,
      enable_drop = FALSE
    ),
    custom_parameter_list
  )
)

custom_summary <- custom_sim$summary

head(custom_summary)
```

    ##   fig_id fig_diameter flower_count resource_use richness_skipped
    ## 1      1     1.827429         1356          705                0
    ## 2      2     2.223787         1643          111                0
    ## 3      3     4.000000         4602          678                0
    ## 4      4     2.584610         2138          695                0
    ## 5      5     2.655145         1617          265                0
    ## 6      6     4.000000         2412          538                0
    ##   entry_Pollinator_A eggs_Pollinator_A entry_Pollinator_B eggs_Pollinator_B
    ## 1                  6               473                  1                92
    ## 2                  0                 0                  3                10
    ## 3                  5               289                  1                96
    ## 4                  5               287                  3               156
    ## 5                  3               212                  0                 0
    ## 6                  3               216                  3               212
    ##   entry_Galler_A eggs_Galler_A entry_Galler_B eggs_Galler_B entry_Galler_C
    ## 1              2            21              8           111              2
    ## 2              4            23              5            78              0
    ## 3              6           146              4            85              4
    ## 4              5           123              4           103              1
    ## 5              1             0              5            24              3
    ## 6              1            22              0             0              6
    ##   eggs_Galler_C entry_Parasitoid_A eggs_Parasitoid_A entry_Parasitoid_B
    ## 1             8                  0                 0                  0
    ## 2             0                  0                 0                  2
    ## 3            62                  3                37                  2
    ## 4            26                  2                 1                  5
    ## 5            29                  2                 0                  2
    ## 6            88                  3                 0                  1
    ##   eggs_Parasitoid_B entry_Parasitoid_C eggs_Parasitoid_C
    ## 1                 0                  1                 0
    ## 2                 0                  0                 0
    ## 3                 6                  2                 7
    ## 4                21                  0                 0
    ## 5                 6                  2                 0
    ## 6                10                  1                 0
    ##   entry_Hyperparasitoid_A eggs_Hyperparasitoid_A entry_Hyperparasitoid_B
    ## 1                       0                      0                       1
    ## 2                       0                      0                       0
    ## 3                       0                      0                       1
    ## 4                       3                      0                       0
    ## 5                       0                      0                       1
    ## 6                       0                      0                       4
    ##   eggs_Hyperparasitoid_B emergence_Pollinator_A emergence_Pollinator_B
    ## 1                      0                    473                     92
    ## 2                      0                      0                     10
    ## 3                      3                    282                     96
    ## 4                      0                    287                    156
    ## 5                      0                    212                      0
    ## 6                      9                    216                    212
    ##   emergence_Galler_A emergence_Galler_B emergence_Galler_C
    ## 1                 21                111                  8
    ## 2                 23                 78                  0
    ## 3                109                 81                 60
    ## 4                122                 82                 26
    ## 5                  0                 18                 29
    ## 6                 22                  0                 78
    ##   emergence_Parasitoid_A emergence_Parasitoid_B emergence_Parasitoid_C
    ## 1                      0                      0                      0
    ## 2                      0                      0                      0
    ## 3                     37                      6                      4
    ## 4                      1                     21                      0
    ## 5                      0                      6                      0
    ## 6                      0                      1                      0
    ##   emergence_Hyperparasitoid_A emergence_Hyperparasitoid_B unoccupied_flowers
    ## 1                           0                           0                651
    ## 2                           0                           0               1532
    ## 3                           0                           3               3924
    ## 4                           0                           0               1443
    ## 5                           0                           0               1352
    ## 6                           0                           9               1874
    ##   seeds failed_ovules used_flowers resource_ratio sink_strength sink_prop
    ## 1   632            19         1356     0.51991150        1653.0 0.8126844
    ## 2  1489            43         1643     0.06755934        2344.5 0.9513086
    ## 3  3841            83         4602     0.14732725        6439.5 0.9328553
    ## 4  1423            20         2138     0.32507016        2829.5 0.8822887
    ## 5  1330            22         1617     0.16388374        2260.0 0.9317666
    ## 6  1832            42         2412     0.22305141        3286.0 0.9082366
    ##   sink_below_min sink_above_max drop_prob is_dropped
    ## 1          FALSE          FALSE         0          0
    ## 2          FALSE           TRUE         0          0
    ## 3          FALSE          FALSE         0          0
    ## 4          FALSE          FALSE         0          0
    ## 5          FALSE          FALSE         0          0
    ## 6          FALSE          FALSE         0          0
    ##                   entry_matrix
    ## 1 6, 1, 2, 8, 2, 0, 0, 1, 0, 1
    ## 2 0, 3, 4, 5, 0, 0, 2, 0, 0, 0
    ## 3 5, 1, 6, 4, 4, 3, 2, 2, 0, 1
    ## 4 5, 3, 5, 4, 1, 2, 5, 0, 3, 0
    ## 5 3, 0, 1, 5, 3, 2, 2, 2, 0, 1
    ## 6 3, 3, 1, 0, 6, 3, 1, 1, 0, 4

## two-species community

Here, we have a minimal two-species example.

``` r

two_species <- c("Pollinator_A", "Galler_A")

two_roles <- list(
  guild = c(
    Pollinator_A = "pollinator",
    Galler_A     = "galler"
  ),
  hosts = list(
    Pollinator_A = character(0),
    Galler_A     = character(0)
  ),
  parasitoid = list(
    Pollinator_A = character(0),
    Galler_A     = character(0)
  )
)

two_parameter_list <- list(
  entry_mu = c(
    Pollinator_A = 8,
    Galler_A     = 5
  ),

  entry_size = c(
    Pollinator_A = 6,
    Galler_A     = 5
  ),

  fecundity_mean = c(
    Pollinator_A = 80,
    Galler_A     = 25
  ),

  fecundity_dispersion = c(
    Pollinator_A = 1.2,
    Galler_A     = 1.5
  ),

  egg_success_prob = c(
    Pollinator_A = 0.9,
    Galler_A     = 0.7
  ),

  egg_success_prob_by_phase = list(
    Pollinator_A = c(phase2 = 0.9),
    Galler_A     = c(phase1 = 0.7)
  ),

  layer_preference = list(
    Pollinator_A = c(core = 0.7, mid = 0.2, outer = 0.1),
    Galler_A     = c(core = 0.4, mid = 0.4, outer = 0.2)
  ),

  max_entry_table = c(
    Pollinator_A = 20,
    Galler_A     = 20
  ),

  parasitism_prob = c(),

  entry_priority = list(
    phase1 = c("Galler_A"),
    phase2 = c("Pollinator_A")
  ),

  interaction_matrix = matrix(
    0,
    nrow = length(two_species),
    ncol = length(two_species),
    dimnames = list(two_species, two_species)
  ),

  interaction_weight = 0,

  species_roles = two_roles
)

validate_figsim_parameter_list(two_parameter_list)
```

    ## Parameter list passed all checks.

``` r

two_sim <- do.call(
  simulate_figwasp_community,
  c(
    list(
      num_figs = 100,
      seed = 123,
      use_layering = TRUE,
      use_layer_preference = TRUE,
      enable_drop = TRUE
    ),
    two_parameter_list
  )
)

two_summary <- two_sim$summary
head(two_summary)
```

    ##   fig_id fig_diameter flower_count resource_use richness_skipped
    ## 1      1     1.827429         1083          388                0
    ## 2      2     2.223787         1668           73                0
    ## 3      3     4.000000         4554          978                0
    ## 4      4     2.584610         1764          197                0
    ## 5      5     2.655145         2061          326                0
    ## 6      6     4.000000         3919          508                0
    ##   entry_Pollinator_A eggs_Pollinator_A entry_Galler_A eggs_Galler_A
    ## 1                  5               388              1             0
    ## 2                  4                73              1             0
    ## 3                 10               945              3            33
    ## 4                  4               194              1             3
    ## 5                  6               246              4            80
    ## 6                  6               441              6            67
    ##   emergence_Pollinator_A emergence_Galler_A unoccupied_flowers seeds
    ## 1                    388                  0                695   682
    ## 2                     73                  0               1595  1561
    ## 3                    945                 33               3576  3508
    ## 4                    194                  3               1567  1525
    ## 5                    246                 80               1735  1703
    ## 6                    441                 67               3411  3351
    ##   failed_ovules used_flowers resource_ratio sink_strength sink_prop
    ## 1            13         1083     0.35826408            NA        NA
    ## 2            34         1668     0.04376499            NA        NA
    ## 3            68         4554     0.21475626            NA        NA
    ## 4            42         1764     0.11167800            NA        NA
    ## 5            32         2061     0.15817564            NA        NA
    ## 6            60         3919     0.12962490            NA        NA
    ##   sink_below_min sink_above_max    drop_prob is_dropped entry_matrix
    ## 1             NA             NA 0.0119222005          0         5, 1
    ## 2             NA             NA 0.0005193827          0         4, 1
    ## 3             NA             NA 0.0028646583          0        10, 3
    ## 4             NA             NA 0.0010237896          0         4, 1
    ## 5             NA             NA 0.0016288618          0         6, 4
    ## 6             NA             NA 0.0012248015          0         6, 6

## Practical limitations

Use this package to other Ficus-fig wasp systems, the prerequisite is
that you need to prepare the input data following the package
guidelines, especially species roles and guilds, basic biological
information, such as expected oviposition or egg-laying rates,
interactions among species, host-parasitoid relationships. In fact, for
many fig wasp communities, this information is scattered across
different publications or remains unclear, which may limit practical
application. However, this is a limitation of available biological
information rather than a restriction of the package.
