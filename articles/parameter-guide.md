# Parameter guide

## Overview

This guide explains the major biological and statistical parameters used
by
[`simulate_figwasp_community()`](https://dongyiyi.github.io/figsimR/reference/simulate_figwasp_community.md).

``` r

library(figsimR)
```

    ## 
    ## Attaching package: 'figsimR'

    ## The following object is masked from 'package:base':
    ## 
    ##     %||%

``` r

library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r

library(forcats)
library(ggplot2)
library(scales)
library(DiagrammeR)
library(vegan)
```

    ## Loading required package: permute

``` r

library(patchwork)
library(purrr)
```

    ## 
    ## Attaching package: 'purrr'

    ## The following object is masked from 'package:scales':
    ## 
    ##     discard

    ## The following object is masked from 'package:figsimR':
    ## 
    ##     %||%

``` r

library(stringr)
library(tibble)
library(readr)
```

    ## 
    ## Attaching package: 'readr'

    ## The following object is masked from 'package:scales':
    ## 
    ##     col_factor

``` r

library(ggpubr)
```

    ## Registered S3 methods overwritten by 'car':
    ##   method       from
    ##   hist.boot    FSA 
    ##   confint.boot FSA

``` r

library(Hmsc)       # For Hierarchical Modelling of Species Communities (HMSC)
```

    ## Loading required package: coda

``` r

library(future)      # for parallel backend
library(progressr)   # to display progress bar during long computation
library(tidyr)

# Set custom font (Arial) for figure output
# This section ensures consistent font rendering across figures
# NOTE: This step may be platform-specific (Windows). Users on macOS/Linux may need to adjust font path or use a different font.
library(extrafont)
```

    ## Registering fonts with R

``` r

library(showtext)
```

    ## Loading required package: sysfonts

    ## Loading required package: showtextdb

    ## 
    ## Attaching package: 'showtextdb'

    ## The following object is masked from 'package:extrafont':
    ## 
    ##     font_install

``` r

###################### Data inputs & simulations ##################

set.seed(42)

data("observed_data")
data("species_list")
data("parameter_list_default")

# Prepare parameter ranges for model resampling or optimization
new_param_ranges <- parameter_list_default

# Remove entries not relevant for optimization (e.g., fixed interaction matrix)
new_param_ranges[c("interaction_matrix", "interaction_weight")] <- NULL

# Convert the cleaned parameter list into a structured format suitable for optimization
# This typically standardizes into a list of vectors/ranges for each parameter
new_param_ranges <- convert_parameter_list_to_param_ranges(new_param_ranges)

# Check the resulting structure of the parameter ranges
str(new_param_ranges)
```

    ## List of 60
    ##  $ entry_mu_Sycophaga_testacea                         : num [1:2] 4.2 7.8
    ##  $ entry_size_Sycophaga_testacea                       : num [1:2] 6.3 11.7
    ##  $ fecundity_mean_Sycophaga_testacea                   : num [1:2] 15.4 28.6
    ##  $ fecundity_dispersion_Sycophaga_testacea             : num [1:2] 1.05 1.95
    ##  $ egg_success_prob_Sycophaga_testacea                 : num [1:2] 0.56 1
    ##  $ layer_core_raw_Sycophaga_testacea                   : num [1:2] 0.49 0.91
    ##  $ layer_mid_raw_Sycophaga_testacea                    : num [1:2] 0.175 0.325
    ##  $ layer_outer_raw_Sycophaga_testacea                  : num [1:2] 0.035 0.065
    ##  $ entry_mu_Apocrypta_sp                               : num [1:2] 7 13
    ##  $ entry_size_Apocrypta_sp                             : num [1:2] 5.6 10.4
    ##  $ fecundity_mean_Apocrypta_sp                         : num [1:2] 2.8 5.2
    ##  $ fecundity_dispersion_Apocrypta_sp                   : num [1:2] 1.54 2.86
    ##  $ egg_success_prob_Apocrypta_sp                       : num [1:2] 0.42 0.78
    ##  $ layer_core_raw_Apocrypta_sp                         : num [1:2] 0.035 0.065
    ##  $ layer_mid_raw_Apocrypta_sp                          : num [1:2] 0.07 0.13
    ##  $ layer_outer_raw_Apocrypta_sp                        : num [1:2] 0.595 1
    ##  $ parasitism_prob_Apocrypta_sp                        : num [1:2] 0.595 1
    ##  $ entry_mu_Sycophaga_mayri                            : num [1:2] 6.3 11.7
    ##  $ entry_size_Sycophaga_mayri                          : num [1:2] 6.3 11.7
    ##  $ fecundity_mean_Sycophaga_mayri                      : num [1:2] 14.7 27.3
    ##  $ fecundity_dispersion_Sycophaga_mayri                : num [1:2] 1.05 1.95
    ##  $ egg_success_prob_Sycophaga_mayri                    : num [1:2] 0.42 0.78
    ##  $ layer_core_raw_Sycophaga_mayri                      : num [1:2] 0.28 0.52
    ##  $ layer_mid_raw_Sycophaga_mayri                       : num [1:2] 0.28 0.52
    ##  $ layer_outer_raw_Sycophaga_mayri                     : num [1:2] 0.014 0.026
    ##  $ entry_mu_Ceratosolen_sp                             : num [1:2] 9.1 16.9
    ##  $ entry_size_Ceratosolen_sp                           : num [1:2] 5.6 10.4
    ##  $ fecundity_mean_Ceratosolen_sp                       : num [1:2] 56 104
    ##  $ fecundity_dispersion_Ceratosolen_sp                 : num [1:2] 0.7 1.3
    ##  $ egg_success_prob_Ceratosolen_sp                     : num [1:2] 0.63 1
    ##  $ layer_core_raw_Ceratosolen_sp                       : num [1:2] 0.49 0.91
    ##  $ layer_mid_raw_Ceratosolen_sp                        : num [1:2] 0.07 0.13
    ##  $ layer_outer_raw_Ceratosolen_sp                      : num [1:2] 0.01 0.013
    ##  $ entry_mu_Sycophaga_agraensis                        : num [1:2] 4.2 7.8
    ##  $ entry_size_Sycophaga_agraensis                      : num [1:2] 4.9 9.1
    ##  $ fecundity_mean_Sycophaga_agraensis                  : num [1:2] 2.1 3.9
    ##  $ fecundity_dispersion_Sycophaga_agraensis            : num [1:2] 1.05 1.95
    ##  $ egg_success_prob_Sycophaga_agraensis                : num [1:2] 0.28 0.52
    ##  $ layer_core_raw_Sycophaga_agraensis                  : num [1:2] 0.01 0.013
    ##  $ layer_mid_raw_Sycophaga_agraensis                   : num [1:2] 0.35 0.65
    ##  $ layer_outer_raw_Sycophaga_agraensis                 : num [1:2] 0.343 0.637
    ##  $ parasitism_prob_Sycophaga_agraensis                 : num [1:2] 0.63 1
    ##  $ entry_mu_Apocrypta_westwoodi                        : num [1:2] 4.2 7.8
    ##  $ entry_size_Apocrypta_westwoodi                      : num [1:2] 4.9 9.1
    ##  $ fecundity_mean_Apocrypta_westwoodi                  : num [1:2] 3.5 6.5
    ##  $ fecundity_dispersion_Apocrypta_westwoodi            : num [1:2] 0.7 1.3
    ##  $ egg_success_prob_Apocrypta_westwoodi                : num [1:2] 0.42 0.78
    ##  $ layer_core_raw_Apocrypta_westwoodi                  : num [1:2] 0.01 0.013
    ##  $ layer_mid_raw_Apocrypta_westwoodi                   : num [1:2] 0.35 0.65
    ##  $ layer_outer_raw_Apocrypta_westwoodi                 : num [1:2] 0.343 0.637
    ##  $ parasitism_prob_Apocrypta_westwoodi                 : num [1:2] 0.42 0.78
    ##  $ egg_success_prob_by_phase_Apocrypta_sp_phase1       : num [1:2] 0.28 0.52
    ##  $ egg_success_prob_by_phase_Apocrypta_sp_phase2       : num [1:2] 0.28 0.52
    ##  $ egg_success_prob_by_phase_Apocrypta_westwoodi_phase2: num [1:2] 0.07 0.13
    ##  $ egg_success_prob_by_phase_Apocrypta_westwoodi_phase3: num [1:2] 0.49 0.91
    ##  $ egg_success_prob_by_phase_Sycophaga_agraensis_phase2: num [1:2] 0.07 0.13
    ##  $ egg_success_prob_by_phase_Sycophaga_agraensis_phase3: num [1:2] 0.56 1
    ##  $ egg_success_prob_by_phase_Sycophaga_mayri_phase2    : num [1:2] 0.63 1
    ##  $ egg_success_prob_by_phase_Ceratosolen_sp_phase2     : num [1:2] 0.56 1
    ##  $ egg_success_prob_by_phase_Sycophaga_testacea_phase1 : num [1:2] 0.56 1

``` r

obs_comm_mat <- observed_data

# resample observed data and calculate metrics
file_path <- system.file("extdata", "resample_observed_metrics_df_results.rds", package = "figsimR")

if (file_path != "" && file.exists(file_path)) {
  message("Loading precomputed resampled simulated result...")
  resample_observed_metrics_df_results <- readRDS(file_path)
} else {
  message("Running resampling may take several minutes...")
  resample_observed_metrics_df_results  <- resample_observed_all_metrics(
  observed_df = observed_data,       # empirical data
  seed = 42,                         # set seed for reproducibility
  wasp_cols = species_list,          # columns corresponding to species counts
  n_draws = 500,                     # number of bootstrap replicates
  sample_n = 200,                    # number of figs sampled per replicate
  calc_func = calc_all_metrics       # metric calculation function
  )
}
```

    ## Loading precomputed resampled simulated result...

``` r

head(resample_observed_metrics_df_results)
```

    ##   mean_richness mean_shannon mean_simpson mean_evenness mean_bray_curtis
    ## 1         4.785    0.8048667    0.4238703     0.5216083        0.6382065
    ## 2         4.815    0.8081919    0.4212343     0.5249506        0.6505153
    ## 3         4.775    0.7528049    0.3957639     0.4889228        0.6051254
    ## 4         4.755    0.7880249    0.4160750     0.5102444        0.6386990
    ## 5         4.870    0.8167150    0.4188874     0.5198021        0.6279106
    ## 6         4.940    0.8147954    0.4223479     0.5186080        0.6317665
    ##   mean_jaccard nestedness connectance links_per_species modularity   source
    ## 1    0.3047061   29.04625   0.8055556          84.05434 0.06025393 Observed
    ## 2    0.3101264   34.97029   0.8065327          84.17134 0.05209199 Observed
    ## 3    0.3024372   32.24482   0.8038721          84.16754 0.06811655 Observed
    ## 4    0.3071446   32.74933   0.8045685          83.50999 0.06081097 Observed
    ## 5    0.2942915   33.18542   0.8116667          85.50308 0.06205649 Observed
    ## 6    0.2763337   33.63853   0.8274707          86.35830 0.05504669 Observed

``` r

# Summary statistics across all bootstrap replicates
summary(resample_observed_metrics_df_results )
```

    ##  mean_richness    mean_shannon     mean_simpson    mean_evenness   
    ##  Min.   :4.650   Min.   :0.7161   Min.   :0.3751   Min.   :0.4660  
    ##  1st Qu.:4.810   1st Qu.:0.7899   1st Qu.:0.4107   1st Qu.:0.5068  
    ##  Median :4.860   Median :0.8055   Median :0.4202   Median :0.5161  
    ##  Mean   :4.859   Mean   :0.8053   Mean   :0.4197   Mean   :0.5157  
    ##  3rd Qu.:4.910   3rd Qu.:0.8222   3rd Qu.:0.4282   3rd Qu.:0.5254  
    ##  Max.   :5.070   Max.   :0.8761   Max.   :0.4622   Max.   :0.5597  
    ##  mean_bray_curtis  mean_jaccard      nestedness     connectance    
    ##  Min.   :0.5887   Min.   :0.2507   Min.   :28.03   Min.   :0.7825  
    ##  1st Qu.:0.6181   1st Qu.:0.2816   1st Qu.:32.47   1st Qu.:0.8066  
    ##  Median :0.6270   Median :0.2933   Median :34.02   Median :0.8149  
    ##  Mean   :0.6267   Mean   :0.2932   Mean   :33.96   Mean   :0.8151  
    ##  3rd Qu.:0.6360   3rd Qu.:0.3036   3rd Qu.:35.45   3rd Qu.:0.8240  
    ##  Max.   :0.6576   Max.   :0.3403   Max.   :39.58   Max.   :0.8492  
    ##  links_per_species   modularity            source   
    ##  Min.   :81.75     Min.   :0.03600   Length   :500  
    ##  1st Qu.:84.32     1st Qu.:0.05366   N.unique :  1  
    ##  Median :85.11     Median :0.06003   N.blank  :  0  
    ##  Mean   :85.12     Mean   :0.05869   Min.nchar:  8  
    ##  3rd Qu.:85.94     3rd Qu.:0.06479   Max.nchar:  8  
    ##  Max.   :88.11     Max.   :0.07826

``` r

# Optionally extract a cleaner subset of core diversity/network metrics LHS optimization
observed_result_df_clean <- resample_observed_metrics_df_results[, 1:10]
head(resample_observed_metrics_df_results)
```

    ##   mean_richness mean_shannon mean_simpson mean_evenness mean_bray_curtis
    ## 1         4.785    0.8048667    0.4238703     0.5216083        0.6382065
    ## 2         4.815    0.8081919    0.4212343     0.5249506        0.6505153
    ## 3         4.775    0.7528049    0.3957639     0.4889228        0.6051254
    ## 4         4.755    0.7880249    0.4160750     0.5102444        0.6386990
    ## 5         4.870    0.8167150    0.4188874     0.5198021        0.6279106
    ## 6         4.940    0.8147954    0.4223479     0.5186080        0.6317665
    ##   mean_jaccard nestedness connectance links_per_species modularity   source
    ## 1    0.3047061   29.04625   0.8055556          84.05434 0.06025393 Observed
    ## 2    0.3101264   34.97029   0.8065327          84.17134 0.05209199 Observed
    ## 3    0.3024372   32.24482   0.8038721          84.16754 0.06811655 Observed
    ## 4    0.3071446   32.74933   0.8045685          83.50999 0.06081097 Observed
    ## 5    0.2942915   33.18542   0.8116667          85.50308 0.06205649 Observed
    ## 6    0.2763337   33.63853   0.8274707          86.35830 0.05504669 Observed

``` r

file_path <- system.file("extdata", "lhs_optimize_result.rds", package = "figsimR")

# Load precomputed result if available; otherwise run the full optimization (may take ~1 hour)
if (file_path != "" && file.exists(file_path)) {
  message("Loading precomputed LHS optimization result...")
  lhs_optimize_result <- readRDS(file_path)
} else {
  message("Running LHS optimization from scratch (this may take several hours; under testing: running 50 minetes with 62 cores)...")
  lhs_optimize_result <- explore_parameter_space_lhs(
    new_param_ranges = new_param_ranges,                  # parameter ranges to explore
    observed_summary = observed_result_df_clean,          # summary statistics from observed data
    num_figs = 1000,                                      # number of figs to simulate per run
    n_draws = 500,                                        # bootstrap iterations for metric calculation
    sample_n = 200,                                       # number of figs drawn per bootstrap sample
    n_samples = 50000,                                    # total number of parameter combinations to test
    wasp_cols = species_list,                             # columns to evaluate in diversity metrics
    n_cores = 7,                                          # number of CPU cores to use in parallel
    parameter_list_template = parameter_list_default      # base parameter list structure
  )
}
```

    ## Loading precomputed LHS optimization result...

``` r

# Extract the single best-performing parameter row (minimum loss)
# Rebuild the full parameter list object from the best result row
lhs_optimize_parameter <- rebuild_parameter_list_from_row(lhs_optimize_result[1, ],
                                                          parameter_list_default)
print(lhs_optimize_parameter)
```

    ## $entry_mu
    ##  Sycophaga_testacea        Apocrypta_sp     Sycophaga_mayri      Ceratosolen_sp 
    ##            6.866784           12.231390            7.971764           11.245201 
    ## Sycophaga_agraensis Apocrypta_westwoodi 
    ##            5.345928            6.661942 
    ## 
    ## $entry_size
    ##  Sycophaga_testacea        Apocrypta_sp     Sycophaga_mayri      Ceratosolen_sp 
    ##            8.600065            9.419983           10.433115            8.888246 
    ## Sycophaga_agraensis Apocrypta_westwoodi 
    ##            7.889184            8.440404 
    ## 
    ## $fecundity_mean
    ##  Sycophaga_testacea        Apocrypta_sp     Sycophaga_mayri      Ceratosolen_sp 
    ##           21.299431            3.399390           19.494681           90.098359 
    ## Sycophaga_agraensis Apocrypta_westwoodi 
    ##            3.370820            4.275399 
    ## 
    ## $fecundity_dispersion
    ##  Sycophaga_testacea        Apocrypta_sp     Sycophaga_mayri      Ceratosolen_sp 
    ##            1.729786            1.962193            1.467205            1.292465 
    ## Sycophaga_agraensis Apocrypta_westwoodi 
    ##            1.747348            1.298638 
    ## 
    ## $egg_success_prob
    ##  Sycophaga_testacea        Apocrypta_sp     Sycophaga_mayri      Ceratosolen_sp 
    ##           0.8614453           0.4386163           0.5173069           0.6489697 
    ## Sycophaga_agraensis Apocrypta_westwoodi 
    ##           0.3797118           0.4249933 
    ## 
    ## $egg_success_prob_by_phase
    ## $egg_success_prob_by_phase$Apocrypta_sp
    ##    phase1    phase2 
    ## 0.4298199 0.3413492 
    ## 
    ## $egg_success_prob_by_phase$Apocrypta_westwoodi
    ##    phase2    phase3 
    ## 0.0947316 0.6109436 
    ## 
    ## $egg_success_prob_by_phase$Sycophaga_agraensis
    ##    phase2    phase3 
    ## 0.1265416 0.6976915 
    ## 
    ## $egg_success_prob_by_phase$Sycophaga_mayri
    ##    phase2 
    ## 0.7154422 
    ## 
    ## $egg_success_prob_by_phase$Ceratosolen_sp
    ##    phase2 
    ## 0.6625057 
    ## 
    ## $egg_success_prob_by_phase$Sycophaga_testacea
    ##    phase1 
    ## 0.6818591 
    ## 
    ## 
    ## $layer_preference
    ## $layer_preference$Ceratosolen_sp
    ##       core        mid      outer 
    ## 0.65321622 0.12383243 0.01016421 
    ## 
    ## $layer_preference$Sycophaga_mayri
    ##       core        mid      outer 
    ## 0.38104585 0.46268065 0.02136136 
    ## 
    ## $layer_preference$Sycophaga_testacea
    ##       core        mid      outer 
    ## 0.88792815 0.21047071 0.05052812 
    ## 
    ## $layer_preference$Apocrypta_sp
    ##       core        mid      outer 
    ## 0.04705676 0.07985802 0.92596364 
    ## 
    ## $layer_preference$Apocrypta_westwoodi
    ##       core        mid      outer 
    ## 0.01189974 0.59900773 0.59935123 
    ## 
    ## $layer_preference$Sycophaga_agraensis
    ##       core        mid      outer 
    ## 0.01058241 0.39322237 0.56093442 
    ## 
    ## 
    ## $max_entry_table
    ##  Sycophaga_testacea        Apocrypta_sp     Sycophaga_mayri      Ceratosolen_sp 
    ##                  20                  20                  20                  20 
    ## Sycophaga_agraensis Apocrypta_westwoodi 
    ##                  20                  20 
    ## 
    ## $parasitism_prob
    ##        Apocrypta_sp Sycophaga_agraensis Apocrypta_westwoodi 
    ##           0.9456121           0.8170753           0.4585077 
    ## 
    ## $entry_priority
    ## $entry_priority$phase1
    ## [1] "Sycophaga_testacea" "Apocrypta_sp"      
    ## 
    ## $entry_priority$phase2
    ## [1] "Apocrypta_sp"        "Sycophaga_mayri"     "Ceratosolen_sp"     
    ## [4] "Sycophaga_agraensis" "Apocrypta_westwoodi"
    ## 
    ## $entry_priority$phase3
    ## [1] "Sycophaga_agraensis" "Apocrypta_westwoodi"
    ## 
    ## 
    ## $interaction_matrix
    ##                     Sycophaga_testacea Apocrypta_sp Sycophaga_mayri
    ## Sycophaga_testacea                   0            0               0
    ## Apocrypta_sp                         0            0               0
    ## Sycophaga_mayri                      0            0               0
    ## Ceratosolen_sp                       0            0               0
    ## Sycophaga_agraensis                  0            0               0
    ## Apocrypta_westwoodi                  0            0               0
    ##                     Ceratosolen_sp Sycophaga_agraensis Apocrypta_westwoodi
    ## Sycophaga_testacea               0                   0                   0
    ## Apocrypta_sp                     0                   0                   0
    ## Sycophaga_mayri                  0                   0                   0
    ## Ceratosolen_sp                   0                   0                   0
    ## Sycophaga_agraensis              0                   0                   0
    ## Apocrypta_westwoodi              0                   0                   0
    ## 
    ## $interaction_weight
    ## [1] 0
    ## 
    ## $species_roles
    ## $species_roles$guild
    ##      Ceratosolen_sp     Sycophaga_mayri  Sycophaga_testacea        Apocrypta_sp 
    ##        "pollinator"            "galler"            "galler"        "parasitoid" 
    ## Apocrypta_westwoodi Sycophaga_agraensis 
    ##        "parasitoid"        "parasitoid" 
    ## 
    ## $species_roles$hosts
    ## $species_roles$hosts$Ceratosolen_sp
    ## character(0)
    ## 
    ## $species_roles$hosts$Sycophaga_mayri
    ## character(0)
    ## 
    ## $species_roles$hosts$Sycophaga_testacea
    ## character(0)
    ## 
    ## $species_roles$hosts$Apocrypta_sp
    ## [1] "Sycophaga_testacea"
    ## 
    ## $species_roles$hosts$Apocrypta_westwoodi
    ## [1] "Sycophaga_mayri"
    ## 
    ## $species_roles$hosts$Sycophaga_agraensis
    ## [1] "Ceratosolen_sp"
    ## 
    ## 
    ## $species_roles$parasitoid
    ## $species_roles$parasitoid$Ceratosolen_sp
    ## [1] "Sycophaga_agraensis"
    ## 
    ## $species_roles$parasitoid$Sycophaga_mayri
    ## [1] "Apocrypta_westwoodi"
    ## 
    ## $species_roles$parasitoid$Sycophaga_testacea
    ## [1] "Apocrypta_sp"
    ## 
    ## $species_roles$parasitoid$Apocrypta_sp
    ## character(0)
    ## 
    ## $species_roles$parasitoid$Apocrypta_westwoodi
    ## character(0)
    ## 
    ## $species_roles$parasitoid$Sycophaga_agraensis
    ## character(0)

``` r

# simulate communities with the optimized parameters, keep consistency of the other parameters between explore_parameter_space_lhs() and simulate_figwasp_community()
sim_df <- simulate_figwasp_community(
  num_figs = 1000,

  # default fig / resource settings used inside LHS unless otherwise specified
  fig_diameter_mean = 2.5,
  fig_diameter_sd = 1.2,
  k = 300,
  alpha = 1.3,
  fig_diameter_min = 1.3,
  fig_diameter_max = 4.0,
  use_flower_limit = TRUE,

  # optimized or template parameters
  fecundity_mean = lhs_optimize_parameter$fecundity_mean,
  fecundity_dispersion = lhs_optimize_parameter$fecundity_dispersion,
  entry_mu = lhs_optimize_parameter$entry_mu,
  entry_size = lhs_optimize_parameter$entry_size,
  entry_priority = lhs_optimize_parameter$entry_priority,
  max_entry_table = lhs_optimize_parameter$max_entry_table,
  species_roles = lhs_optimize_parameter$species_roles,
  interaction_matrix = lhs_optimize_parameter$interaction_matrix,
  interaction_weight = lhs_optimize_parameter$interaction_weight,
  egg_success_prob = lhs_optimize_parameter$egg_success_prob,
  egg_success_prob_by_phase = lhs_optimize_parameter$egg_success_prob_by_phase,
  layer_preference = lhs_optimize_parameter$layer_preference,

  # important: match the LHS call
  parasitism_prob = NULL,
  enable_drop = TRUE,
  use_supplemental_parasitism = FALSE,
  drop_cancels_emergence = FALSE,
  use_layering = TRUE,

  # implicit defaults used by simulate_figwasp_community()
  use_layer_preference = TRUE,
  use_egg_success_by_phase = TRUE,
  entry_distribution = "lognormal",
  record_individual = FALSE,
  use_sink_strength = FALSE,
  p_pollination_per_ovule = 0.98,
  p_no_entry = 0.002,
  host_sanction = 0.8,

  seed = 42
)

# check simulation results
head(sim_df$summary)
```

    ##   fig_id fig_diameter flower_count resource_use richness_skipped
    ## 1      1     4.000000         5660          286                0
    ## 2      2     1.822362         1432          427                0
    ## 3      3     2.935754         2923          827                0
    ## 4      4     3.259435         2955          621                0
    ## 5      5     2.985122         1908          172                0
    ## 6      6     2.372651         1360          588                0
    ##   entry_Ceratosolen_sp eggs_Ceratosolen_sp entry_Sycophaga_mayri
    ## 1                    7                 192                     3
    ## 2                    4                 328                     3
    ## 3                   11                 611                    11
    ## 4                   16                 413                     5
    ## 5                    1                   0                     5
    ## 6                    6                 537                     2
    ##   eggs_Sycophaga_mayri entry_Sycophaga_testacea eggs_Sycophaga_testacea
    ## 1                   31                        1                      63
    ## 2                   52                        4                      47
    ## 3                  173                        4                      43
    ## 4                  106                        6                     102
    ## 5                   95                        4                      77
    ## 6                    3                        5                      48
    ##   entry_Apocrypta_sp eggs_Apocrypta_sp entry_Apocrypta_westwoodi
    ## 1                 20                22                         5
    ## 2                 13                 4                         8
    ## 3                 12                15                        11
    ## 4                 15                37                         8
    ## 5                 17                 7                         4
    ## 6                 22                12                        12
    ##   eggs_Apocrypta_westwoodi entry_Sycophaga_agraensis eggs_Sycophaga_agraensis
    ## 1                       14                         3                        0
    ## 2                       13                         5                        8
    ## 3                        5                         8                       19
    ## 4                        1                         9                       23
    ## 5                        5                         4                        0
    ## 6                        0                         4                        6
    ##   emergence_Ceratosolen_sp emergence_Sycophaga_mayri
    ## 1                      192                        17
    ## 2                      320                        39
    ## 3                      592                       168
    ## 4                      390                       105
    ## 5                        0                        90
    ## 6                      531                         3
    ##   emergence_Sycophaga_testacea emergence_Apocrypta_sp
    ## 1                           41                     22
    ## 2                           43                      4
    ## 3                           28                     15
    ## 4                           65                     37
    ## 5                           70                      7
    ## 6                           36                     12
    ##   emergence_Apocrypta_westwoodi emergence_Sycophaga_agraensis
    ## 1                            14                             0
    ## 2                            13                             8
    ## 3                             5                            19
    ## 4                             1                            23
    ## 5                             5                             0
    ## 6                             0                             6
    ##   unoccupied_flowers seeds failed_ovules used_flowers resource_ratio
    ## 1               5374  5257           117         5660     0.05053004
    ## 2               1005   975            30         1432     0.29818436
    ## 3               2096  2059            37         2923     0.28292850
    ## 4               2334  2295            39         2955     0.21015228
    ## 5               1736  1691            45         1908     0.09014675
    ## 6                772   761            11         1360     0.43235294
    ##   sink_strength sink_prop sink_below_min sink_above_max    drop_prob is_dropped
    ## 1            NA        NA             NA             NA 0.0005557147          0
    ## 2            NA        NA             NA             NA 0.0065732208          0
    ## 3            NA        NA             NA             NA 0.0056484198          0
    ## 4            NA        NA             NA             NA 0.0027361129          0
    ## 5            NA        NA             NA             NA 0.0008256344          0
    ## 6            NA        NA             NA             NA 0.0246872665          0
    ##           entry_matrix
    ## 1    7, 3, 1, 20, 5, 3
    ## 2    4, 3, 4, 13, 8, 5
    ## 3 11, 11, 4, 12, 11, 8
    ## 4   16, 5, 6, 15, 8, 9
    ## 5    1, 5, 4, 17, 4, 4
    ## 6   6, 2, 5, 22, 12, 4

``` r

# extract communities
sim_df_comm <- sim_df$summary
sim_comm_mat <- sim_df_comm[,c(18:23)]

# clean headers and reorder
colnames(sim_comm_mat) <- gsub("^emergence_", "", colnames(sim_comm_mat))
sim_comm_mat <- sim_comm_mat[, match(colnames(obs_comm_mat), colnames(sim_comm_mat))]
head(obs_comm_mat)
```

    ##   Ceratosolen_sp Sycophaga_testacea Apocrypta_sp Sycophaga_mayri
    ## 1           1163                  3            0               2
    ## 2           1510                109            2              70
    ## 3            744                102            0              39
    ## 4           1494                  9            0               8
    ## 5           1224                 10            0              42
    ## 6            478                 77            1              97
    ##   Sycophaga_agraensis Apocrypta_westwoodi
    ## 1                  32                 189
    ## 2                  22                  45
    ## 3                  62                  12
    ## 4                  54                  72
    ## 5                  14                 103
    ## 6                  39                  32

``` r

head(sim_comm_mat)
```

    ##   Ceratosolen_sp Sycophaga_testacea Apocrypta_sp Sycophaga_mayri
    ## 1            192                 41           22              17
    ## 2            320                 43            4              39
    ## 3            592                 28           15             168
    ## 4            390                 65           37             105
    ## 5              0                 70            7              90
    ## 6            531                 36           12               3
    ##   Sycophaga_agraensis Apocrypta_westwoodi
    ## 1                   0                  14
    ## 2                   8                  13
    ## 3                  19                   5
    ## 4                  23                   1
    ## 5                   0                   5
    ## 6                   6                   0

``` r

stopifnot(identical(colnames(obs_comm_mat), colnames(sim_comm_mat)))
```

## key parameters

For full parameter, please see the help document:
?simulate_figwasp_community()

### alpha

Exponent in the power-law relationship between fig diameter and flower
number. Empirically estimated around 1.2-1.4 for many figs.

### drop_cancels_emergence

Logical. If FALSE, dropped figs produce zero emergence and seeds. If
TRUE, their output remains.

### egg_success_prob

Global probability that a single oviposition attempt results in a
successful egg. Applies when use_egg_success_by_phase = FALSE.

### egg_success_prob_by_phase

Overrides global success probabilities by entry phase. A species-level
list of phase-specific reproductive-success probabilities. For each
species and developmental phase, the supplied value is used as the
probability that an arriving or oviposition-attempting female
contributes eggs during that phase. If no phase-specific value is
supplied, the model uses the baseline egg_success_prob. Realized egg
numbers are subsequently determined by fecundity, host availability,
ovary-layer accessibility, and resource constraints. Format:
species_name = list(phase1 = 0.5, phase2 = 0.7,…).

### enable_drop

Enables fig abortion due to host sanctions. If flower use proportion
exceeds host_sanction, fig may abort. In the output “drop_prob” column:
Ranges from 0 (never drop) to 1 (always drop); “is_dropped” column: no
drop (0) and drop (1).

### entry_distribution

Either “nb” (negative binomial) or “lognormal”. Defines how per-species
entry numbers are simulated: “nb”: entries are drawn from NB( mu =
entry_mu, size = entry_size ). “lognormal”: log(entry_mu) - 0.5 used as
meanlog; sdlog = 1 / sqrt(entry_size); result used as Poisson rate.

### entry_mu

Mean number of individuals of each species attempting to “enter”
(arrival or oviposition-attempt events) each fig. Acts as expected value
for entry distribution.

### entry_priority

Defines the temporal entry order of species. Each element is a phase
(e.g., phase1), with a vector of species names. Earlier phases enter
first and may influence resource or host availability for later phases.

### entry_size

Controls variation in entry distribution (the dispersion of the arrival
or oviposition-attempt distribution): For entry_distribution = “nb”:
used as NB size parameter. For “lognormal”: used to derive sdlog = 1 /
sqrt(entry_size).

### fecundity_dispersion

Dispersion parameter for fecundity (size in NB distribution). Controls
variability in egg output across individuals.

### fecundity_mean

Mean number of eggs laid per female individual, per species. Used as mu
in a negative binomial distribution.

### fig_diameter_max

aximum fig diameter used to cap values (optional).

### fig_diameter_mean

Mean fig diameter (in cm). Used to generate per-fig flower number via
the formula: flower_count = k, diameter^alpha, Gamma(shape = 20, scale =
0.1). Larger diameters yield more ovules.

### fig_diameter_min

Truncate simulated diameters (and thus flower counts) to a biologically
reasonable range.

### fig_diameter_sd

Standard deviation of fig diameter. Each fig’s diameter is drawn from a
normal distribution with this mean and SD.

### host_sanction

Threshold parameter used in the optional host-sanction/fig-abortion
module. Threshold of flower use (e.g., 0.8). If exceeded, the fig may
abort (drop) due to excessive exploitation by wasps. When enable_drop =
TRUE, fig drop probability is calculated as a logistic function of
resource_ratio relative to this threshold. This parameter does not
directly represent the probability of abortion, and it is retained as an
optional alternative fig-retention rule.

### k

Scaling constant for estimating flower number from fig diameter.
Represents average number of flowers per cm^ alpha.

### layer_preference

For each species, a named numeric vector giving probabilities of
oviposition in core, mid, and outer ovary layers. Sum must be 1. Or
“none” is not given species-specific layer preferences; their
oviposition attempts are allocated according to remaining ovary-layer
availability. Used only if use_layering = TRUE.

### species_list

Character vector of species names (used to extract relevant columns).

### species_roles

Defines the ecological roles and interactions of all species. Contains
three components: guild: character vector assigning each species to
“pollinator”, “galler”, or “parasitoid”. hosts: list mapping parasitoids
to their hosts. parasitoid: list mapping hosts to their attacking
parasitoids. This structure defines trophic constraints (e.g.,
parasitoids cannot oviposit unless their hosts are present).

### use_egg_success_by_phase

Use egg success probabilities by phase. Default TRUE.

### use_flower_limit

If TRUE, clamps fig diameter between fig_diameter_min and
fig_diameter_max before calculating flower number.

### use_layer_preference

If TRUE, species will preferentially oviposit in specific layers as
defined in layer_preference.

### use_layering species

If TRUE, fig flowers are divided into spatial layers (core, mid, outer).

### use_supplemental_parasitism

Whether to activate supplemental parasitism. Default FALSE.

### use_sink_strength

If TRUE, compute sink-strength metrics but do NOT alter legacy drop
decision. Adds columns: sink_strength, sink_prop, sink_below_min,
sink_above_max.

### sink_w_gall

Per-ovule sink weight for galled ovules (default 1.0).

### sink_w_seed

Per-ovule sink weight for seeds (default 1.5).

### sink_linear_coef

Global multiplier on sink strength (default 1.0).

### sink_min_prop

Reference lower bound of sink proportion (default 0.20).

### sink_max_prop

Reference upper bound of sink proportion (default 0.95).
