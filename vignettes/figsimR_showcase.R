## ----setup, include=FALSE-----------------------------------------------------
# Global chunk options for a fast, deterministic vignette
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  fig.width = 6,
  fig.height = 4,
  dpi = 120,
  cache = FALSE
)

# We'll keep heavy computations disabled by default in the vignette.
RUN_HEAVY <- FALSE

## ----setup-library, echo=FALSE------------------------------------------------
# set seed and load library in a separate, non-displayed chunk
set.seed(123)
library(figsimR)

## ----load-data----------------------------------------------------------------
# Example data shipped with the package
data(observed_data)          # 935 x 6 abundance matrix (one row per fig)
data(parameter_list_default) # Named lists of biological parameters
data(species_list)           # Character vector of the six focal species

# Inspect the first rows
head(observed_data)

## ----bootstrap-metrics--------------------------------------------------------
obs_metrics <- resample_observed_all_metrics(
  observed_df = observed_data,
  wasp_cols   = species_list,
  n_draws     = 50,     
  sample_n    = 100,    
  seed        = 42,
  calc_func   = calc_all_metrics
)

head(obs_metrics)

## ----param-ranges-------------------------------------------------------------
param_ranges <- parameter_list_default
param_ranges[c("interaction_matrix", "interaction_weight")] <- NULL
param_ranges <- convert_parameter_list_to_param_ranges(param_ranges)
str(param_ranges, max.level = 1)

## ----explore-params-----------------------------------------------------------
# Optional: progress bars (kept OFF in vignette)
# progressr::handlers(global = TRUE); progressr::handlers("txtprogressbar")

# load a precomputed optimization (ship it in inst/extdata/)
lhs_file <- system.file("extdata", "lhs_optimize_result.rds", package = "figsimR")

if (file.exists(lhs_file)) {
  lhs_opt <- readRDS(lhs_file)
} else if (RUN_HEAVY) {
  lhs_opt <- explore_parameter_space_lhs(
    new_param_ranges        = param_ranges,
    observed_summary        = obs_metrics[, 1:10], # numeric summary columns
    num_figs                = 1000,
    n_draws                 = 500,
    sample_n                = 200,
    wasp_cols               = species_list,
    parameter_list_template = parameter_list_default,
    n_samples               = 5000,
    top_k                   = 5,
    n_cores                 = max(1, parallel::detectCores() - 1)
  )
} else {
  
  lhs_opt <- explore_parameter_space_lhs(
    new_param_ranges        = param_ranges,
    observed_summary        = obs_metrics[, 1:10],
    num_figs                = 200,
    n_draws                 = 20,
    sample_n                = 100,
    wasp_cols               = species_list,
    parameter_list_template = parameter_list_default,
    n_samples               = 50,    
    top_k                   = 1,
    n_cores                 = 1
  )
}

best_row <- lhs_opt[1, , drop = FALSE]
best_params <- rebuild_parameter_list_from_row(best_row, parameter_list_default)

## ----simulate-bcm-------------------------------------------------------------
sim_out <- simulate_figwasp_community(
  num_figs                = if (RUN_HEAVY) 1000 else 300,
  fecundity_mean          = best_params$fecundity_mean,
  fecundity_dispersion    = best_params$fecundity_dispersion,
  entry_mu                = best_params$entry_mu,
  entry_size              = best_params$entry_size,
  entry_priority          = best_params$entry_priority,
  species_roles           = best_params$species_roles,
  max_entry_table         = best_params$max_entry_table,
  enable_drop             = TRUE,
  drop_cancels_emergence  = FALSE,
  entry_distribution      = "lognormal",
  interaction_matrix      = best_params$interaction_matrix,
  interaction_weight      = 0,
  egg_success_prob        = best_params$egg_success_prob,
  egg_success_prob_by_phase = best_params$egg_success_prob_by_phase,
  layer_preference        = best_params$layer_preference,
  use_layering            = TRUE,
  seed                    = 42
)

sim_df <- sim_out$summary
head(sim_df)

