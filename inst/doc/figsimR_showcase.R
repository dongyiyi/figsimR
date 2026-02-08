## ----setup_progressr, include=FALSE-------------------------------------------

safe_handlers <- function(...) {
  if (requireNamespace("progressr", quietly = TRUE)) {
    try(safe_handlers(...), silent = TRUE)
  }
}

knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

has_progressr <- requireNamespace("progressr", quietly = TRUE)
# Avoid parallel errors on CRAN/check
if (requireNamespace("future", quietly = TRUE)) {
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    future::plan("sequential")
  } else {
    future::plan("multisession")
  }
}

with_p <- if (has_progressr) progressr::with_progress else function(expr) force(expr)
p_factory <- if (has_progressr) progressr::progressor else function(...) { function(...) {} }

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

# loading packages
library(igraph)
library(ggraph)
library(patchwork)
library(dplyr)
library(stringr)
library(dplyr)
library(stringr)
library(FSA)        
library(ggplot2)
library(RColorBrewer)
library(purrr)
library(tidyr)
library(ggpubr)
library(vegan)     
library(Hmsc)       # For Hierarchical Modelling of Species Communities (HMSC)
library(reshape2)   
# Load required libraries for parallel computing and progress tracking
library(future)      # for parallel backend
library(progressr)   # to display progress bar during long computation
library(tidyverse)


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
  n_draws     = 50,     # small for vignette
  sample_n    = 100,    # small for vignette
  seed        = 42,
  calc_func   = calc_all_metrics
)

head(obs_metrics)

# Summary statistics across all bootstrap replicates
summary(obs_metrics)

# Optionally extract a cleaner subset of core diversity/network metrics
obs_metrics_mean <- obs_metrics[, 1:10]




## ----param-ranges-------------------------------------------------------------
# Prepare parameter ranges for model resampling or optimization
new_param_ranges <- parameter_list_default

# Remove entries not relevant for optimization (e.g., fixed interaction matrix)
new_param_ranges[c("interaction_matrix", "interaction_weight")] <- NULL

# Convert the cleaned parameter list into a structured format suitable for optimization
# This typically standardizes into a list of vectors/ranges for each parameter
new_param_ranges <- convert_parameter_list_to_param_ranges(new_param_ranges)

# Check the resulting structure of the parameter ranges
str(new_param_ranges)

# To avoid long computation time during review (~1 hour), a precomputed optimization result file 
# (`lhs_optimize_result.rds`) is provided. If this file is present in the working directory, it 
# will be loaded automatically. 
# 
# If you prefer to re-run the full Latin Hypercube Sampling optimization yourself, you may comment 
# out the conditional block below and run `explore_parameter_space_lhs()` directly.

# Load required libraries for parallel computing and progress tracking
#library(future)      # for parallel backend
#library(progressr)   # to display progress bar during long computation

# Set up global progress bar handler and seed for reproducibility
safe_handlers(global = TRUE)
safe_handlers("txtprogressbar")  # show text progress bar in console
set.seed(42)

# Load precomputed result if available; otherwise run the full optimization (may take ~1 hour), we run small dataset for testing
if (file.exists("lhs_optimize_result.rds")) {
  message("Loading precomputed LHS optimization result...")
  lhs_optimize_result <- readRDS("lhs_optimize_result.rds")
} else {
  message("Running LHS optimization from scratch (this may take ~1 hour)...")
  lhs_optimize_result <- explore_parameter_space_lhs(
    new_param_ranges = new_param_ranges,                  # parameter ranges to explore
    observed_summary = obs_metrics_mean,                  # summary statistics from observed data
    num_figs = 200,                                      # number of figs to simulate per run
    n_draws = 50,                                        # bootstrap iterations for metric calculation
    sample_n = 20,                                       # number of figs drawn per bootstrap sample
    n_samples = 50,                                      # total number of parameter combinations to test
    wasp_cols = species_list,                             # columns to evaluate in diversity metrics
    n_cores = 4,                                          # number of CPU cores to use in parallel
    parameter_list_template = parameter_list_default      # base parameter list structure
  )
}


# Extract the single best-performing parameter row (minimum loss)
best_lhs_optimize_result <- lhs_optimize_result[1, ]
print(best_lhs_optimize_result)

#saveRDS(lhs_optimize_result, file = "lhs_optimize_result.rds")

# Rebuild the full parameter list object from the best result row
# This list will be used in subsequent simulations
lhs_best_param_list <- rebuild_parameter_list_from_row(
  best_lhs_optimize_result,
  parameter_list_default
)

# Display the reconstructed best-fit parameter list
lhs_best_param_list




## ----explore-params-----------------------------------------------------------
# Optional: progress bars (kept OFF in vignette)
# safe_handlers(global = TRUE); safe_handlers("txtprogressbar")

# Try to load a precomputed optimization (ship it in )
lhs_file <- system.file("extdata", "lhs_optimize_result.rds", package = "figsimR")

if (file.exists(lhs_file)) {
  lhs_opt <- readRDS(lhs_file)
} else if (RUN_HEAVY) {
  lhs_opt <- explore_parameter_space_lhs(
    new_param_ranges        = new_param_ranges,
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
  # Tiny demo search to keep vignette fast; results not scientific
  lhs_opt <- explore_parameter_space_lhs(
    new_param_ranges        = new_param_ranges,
    observed_summary        = obs_metrics[, 1:10],
    num_figs                = 200,
    n_draws                 = 20,
    sample_n                = 100,
    wasp_cols               = species_list,
    parameter_list_template = parameter_list_default,
    n_samples               = 50,    # tiny demo
    top_k                   = 1,
    n_cores                 = 1
  )
}

best_row <- lhs_opt[1, , drop = FALSE]
best_params <- rebuild_parameter_list_from_row(best_row, parameter_list_default)

## ----simulate-fim-------------------------------------------------------------
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

# let plot it
summarize_simulated_metrics(sim_df, species_list, version_label = "lhs_best_param_list")



## ----quantify-gap-------------------------------------------------------------
sim_metrics <- resample_simulated_metrics(
  n_draws   = if (RUN_HEAVY) 500 else 50,
  sample_n  = if (RUN_HEAVY) 200 else 100,
  num_figs  = nrow(sim_df),
  simulator_func = function(...) list(summary = sim_df),  # reuse existing sim
  calc_func = calc_all_metrics,
  wasp_cols = paste0("emergence_", species_list)
)

# Combine tags and run ordination/dispersion tests (example using NMDS + betadisper)
obs_tag <- observed_data
sim_tag <- sim_df[, paste0("emergence_", species_list), drop = FALSE]

obs_tag$group <- "Observed"
sim_tag$group <- "Simulated"

# Ensure column names match for rbind
colnames(sim_tag) <- gsub("^emergence_", "", colnames(sim_tag))

combined <- rbind(obs_tag, sim_tag)

combined <- combined[rowSums(combined[ , species_list]) > 0, ]
sp_mat   <- vegan::decostand(combined[ , species_list], method = "hellinger")
nmds     <- vegan::metaMDS(sp_mat, distance = "euclidean", k = 2, trymax = 10)

# Dispersion comparison
dist_eu  <- vegan::vegdist(sp_mat, "euclidean")
bd       <- vegan::betadisper(dist_eu, group = factor(combined$group))
bd_perm  <- vegan::permutest(bd, permutations = 499)

bd_perm


## ----plot-variation-gap-pca, fig.cap="Figure 2: Visualization of the Variation Gap using Principal Component Analysis (PCA). The dense cloud of blue points represents the 'Theoretical Envelope' generated by the FIM. The widely dispersed red points are the empirical observations from the field, clearly falling outside the theoretical boundary."----

pca <- stats::prcomp(sp_mat, scale. = TRUE)
plot_data <- as.data.frame(pca$x[, 1:2])
plot_data$group <- combined$group


pca_summary <- summary(pca)
variance_explained <- pca_summary$importance[2, 1:2] * 100

# --- PCA ---

centroids <- aggregate(cbind(PC1, PC2) ~ group, data = plot_data, FUN = mean)

ggplot(plot_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.8, size = 2) +
  stat_ellipse(aes(group = group), type = "t", level = 0.95, linetype = "dashed", linewidth = 0.6) +
  geom_point(data = centroids, aes(x = PC1, y = PC2), size = 5, shape = 18) +
  scale_color_manual(values = c("Observed" = "#8da0cb", "Simulated" = "#fc8d62"),
                     name = "Community Source") +
  labs(title = "PCA of Simulated vs. Observed Communities",
       subtitle = "Demonstrating the Variation Gap",
       x = paste0("Principal Component 1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("Principal Component 2 (", round(variance_explained[2], 1), "%)")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  coord_fixed() 


## ----knockouts----------------------------------------------------------------

# module control options, more details: ?simulate_figwasp_community()
# use_layering = TRUE/FALSE
# use_layer_preference = TRUE/FALSE
# enable_drop = TRUE/FALSE
# drop_cancels_emergence = TRUE/FALSE
# use_supplemental_parasitism = TRUE/FALSE
# use_flower_limit = TRUE/FALSE
# use_egg_success_by_phase = TRUE/FALSE

# Example: remove spatial layering
sim_no_space <- simulate_figwasp_community(
  num_figs = if (RUN_HEAVY) 1000 else 300,
  fecundity_mean       = best_params$fecundity_mean,
  fecundity_dispersion = best_params$fecundity_dispersion,
  entry_mu             = best_params$entry_mu,
  entry_size           = best_params$entry_size,
  entry_priority       = best_params$entry_priority,
  species_roles        = best_params$species_roles,
  max_entry_table      = best_params$max_entry_table,
  enable_drop          = TRUE,
  drop_cancels_emergence = FALSE,
  entry_distribution   = "lognormal",
  interaction_matrix   = best_params$interaction_matrix,
  interaction_weight   = 0,
  egg_success_prob     = best_params$egg_success_prob,
  egg_success_prob_by_phase = best_params$egg_success_prob_by_phase,
  layer_preference     = best_params$layer_preference,
  use_layering         = FALSE,  # KO here
  seed                 = 43
)$summary

# Summarize metrics for FIM vs. KO
fim_metrics  <- resample_metrics_from_simulation(sim_df,       n_reps = 50, sample_size = 100, seed = 1)
ko_space   <- resample_metrics_from_simulation(sim_no_space, n_reps = 50, sample_size = 100, seed = 1)

# Loss vs empirical means (smaller = better to observed)
obs_means_row <- as.data.frame(t(colMeans(obs_metrics[ , setdiff(names(obs_metrics), "source") ], na.rm = TRUE)))
fim_loss      <- calculate_metric_loss(fim_metrics, obs_means_row)
ko_space_loss       <- calculate_metric_loss(ko_space,  obs_means_row)

list(FIM_total_loss = fim_loss$total_loss,
     KO_no_space_total_loss = ko_space_loss$total_loss)

# Example: stop use phase-specific success probabilities
sim_no_prob_by_phase <- simulate_figwasp_community(
  num_figs = if (RUN_HEAVY) 1000 else 300,
  fecundity_mean       = best_params$fecundity_mean,
  fecundity_dispersion = best_params$fecundity_dispersion,
  entry_mu             = best_params$entry_mu,
  entry_size           = best_params$entry_size,
  entry_priority       = best_params$entry_priority,
  species_roles        = best_params$species_roles,
  max_entry_table      = best_params$max_entry_table,
  enable_drop          = TRUE,
  drop_cancels_emergence = FALSE,
  entry_distribution   = "lognormal",
  interaction_matrix   = best_params$interaction_matrix,
  interaction_weight   = 0,
  egg_success_prob     = best_params$egg_success_prob,
  egg_success_prob_by_phase = NULL, # turn it off
  use_egg_success_by_phase = FALSE, # turn it off
  layer_preference     = best_params$layer_preference,
  use_layering         = TRUE,  
  seed                 = 43
)$summary

# Summarize metrics for KO prob_by_phase
ko_prob_by_phase   <- resample_metrics_from_simulation(sim_no_prob_by_phase, n_reps = 50, sample_size = 100, seed = 1)

# Loss vs empirical means (smaller = better to observed)
ko_prob_by_phase_loss <- calculate_metric_loss(ko_prob_by_phase,  obs_means_row)

list(FIM_total_loss = fim_loss$total_loss,
     KO_no_space_total_loss = ko_space_loss$total_loss,
     KO_prob_by_phase_total_loss = ko_prob_by_phase_loss$total_loss)



## ----plot-knockouts, fig.cap="Figure 3: Dissection of internal mechanisms via knockout experiments. The bar height represents the increase in model loss (i.e., deviation from empirical means) when a specific mechanism is removed, relative to the full FIM. A larger value indicates a greater contribution of that mechanism to maintaining the theoretical community structure."----

knockout_results <- data.frame(
  mechanism = c(
    "Full Intrinsic Model (FIM)",
    "KO: Spatial Layering",
    "KO: Phase-specific Egg Success"
  ),
  total_loss = c(
    fim_loss$total_loss,
    ko_space_loss$total_loss,
    ko_prob_by_phase_loss$total_loss
  )
)


ggplot(knockout_results, aes(x = reorder(mechanism, total_loss), y = total_loss)) +
  geom_col(fill = "#4682B4", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = round(total_loss, 2)), vjust = -0.5, size = 4.5) +
  labs(
    title = "Model Performance with and without Key Mechanisms",
    subtitle = "Comparing Total Loss against Observed Data",
    x = "Model Configuration",
    y = "Total Loss (vs. Observed Data)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )


## ----strict host-sanction threshold-------------------------------------------
# First, ensure the sanction mechanism *enable_sanctions <- TRUE *
# Run the simulation with this specific sanction rule
# Key change: set a strict host-sanction threshold at 10% of total flowers
# If ovule-occupying eggs (pollinators + gallers) > 10% of flowers -> fig drops -> zero emergence.
# We mark figs as dropped but keep eggs/emergence for inspection (drop_cancels_emergence = FALSE).
# For analysis, we create a masked copy where dropped figs' emergence counts are zeroed.

sim_sanction_10 <- simulate_figwasp_community(
  num_figs = if (RUN_HEAVY) 1000 else 300,
  fecundity_mean       = best_params$fecundity_mean,
  fecundity_dispersion = best_params$fecundity_dispersion,
  entry_mu             = best_params$entry_mu,
  entry_size           = best_params$entry_size,
  entry_priority       = best_params$entry_priority,
  species_roles        = best_params$species_roles,
  max_entry_table      = best_params$max_entry_table,
  enable_drop          = TRUE,
  drop_cancels_emergence = FALSE, # keep it as FALSE: it means "drop => zero emergence"
  entry_distribution   = "lognormal",
  interaction_matrix   = best_params$interaction_matrix,
  interaction_weight   = 0,
  host_sanction = 0.1, # the key change, default = 0.8
  egg_success_prob     = best_params$egg_success_prob,
  egg_success_prob_by_phase = best_params$egg_success_prob_by_phase,
  layer_preference     = best_params$layer_preference,
  use_layering         = TRUE,  # let's turn it on
  seed                 = 43
)

sim_sanction_10 <- sim_sanction_10$summary

# Quick check: how many figs dropped, and confirm dropped figs have zero emergence
summarize_simulated_metrics(sim_sanction_10, species_list, version_label = "sanction 10%")$drop_summary  

summarize_simulated_metrics(sim_sanction_10, species_list, version_label = "sanction 10%")$stats_table


## ----session-info-------------------------------------------------------------
sessionInfo()

