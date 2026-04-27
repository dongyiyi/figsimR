# figsimR

<!-- badges: start -->
[![R-CMD-check](https://github.com/dongyiyi/figsimR/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/dongyiyi/figsimR/actions/workflows/R-CMD-check.yml)
<!-- badges: end -->


figsimR is a mechanistic, agent-based simulation package in R that turns process-based hypotheses into executable, testable models of community assembly

# figsimR: Installation
```
# 1) GitHub (recommended)
# install.packages("remotes")     # if needed
remotes::install_github("dongyiyi/figsimR")
library(figsimR)

# (optional) build vignettes
# remotes::install_github("dongyiyi/figsimR", build_vignettes = TRUE)
# requires: rmarkdown, knitr, pandoc

# 2) Using pak (fast dependency solver)
# install.packages("pak")
pak::pak("dongyiyi/figsimR")
library(figsimR)

# 3) From a local tarball (offline)
# replace with your actual file path:
install.packages("path/to/figsimR_0.1.0.tar.gz", repos = NULL, type = "source")
library(figsimR)

# 4) From a local source directory (after git clone)
# git clone https://github.com/dongyiyi/figsimR.git
# then in R (point to the folder or tar.gz):
# install.packages("devtools")    # if needed
devtools::install_local("path/to/figsimR", upgrade = "never")
# or: devtools::install_local("path/to/figsimR_0.1.0.tar.gz")
library(figsimR)
```

```{r setup-library, echo=FALSE}
# set seed
set.seed(123)

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

```

# 1. Introduction: When theory meets variability

Ecologists have long recognized that natural communities exhibit large variability. `figsimR` operationalizes a mechanistic null workflow: we calibrate a baseline configuration modeling (BCM) to reproduce empirical means, then quantify the Variation Gap—the residual dispersion that not captured by the current BCM configuration. Finally, we use mechanism knockout experiments to dissect which internal processes maintain the theoretical structure.

This vignette shows the full workflow on a fig–wasp case study, with small defaults for speed. For the full reproduction (long runs), see the package website article “Reproducing the case study” (or run this vignette with `RUN_HEAVY <- TRUE`).

Framing & scope. In this vignette, the “Variation Gap” is defined strictly as model–data residual variation under the current BCM configuration. It is not a claim that ecological variability is inherently “extrinsic” or unknowable. Many well-studied drivers—e.g., environmental/landscape heterogeneity, dispersal, demographic/environmental stochasticity, evolutionary change—can be added as explicit modules and may reduce or eliminate the gap. Here we hold these out of scope to (i) delineate the explanatory envelope of intrinsic mechanisms, and (ii) provide a transparent baseline for subsequent extensions.

# 2. Data and quick look

We use an empirical dataset of fig wasp communities collected from *Ficus racemosa*. Each column is a species’ abundance; each row is one fig (community unit).

```{r load-data}
# Example data shipped with the package
data(observed_data)          # 935 x 6 abundance matrix (one row per fig)
data(parameter_list_default) # Named lists of biological parameters
data(species_list)           # Character vector of the six focal species

# Inspect the first rows
head(observed_data)
```

**Provenance.**
Aung, K. M., Chen, H. H., Segar, S. T., Miao, B. G., Peng, Y. Q., & Liu, C. (2022). *Journal of Animal Ecology*, 91, 1303–1315. https://doi.org/10.1111/1365-2656.13701
* **Background consensus.** The field has long developed rich frameworks for environmental heterogeneity, coexistence, and spatial processes. Our contribution here is procedural: to define a transparent baseline workflow that first saturates intrinsic structure, then diagnoses what remains to be explained.

# 3. Build the theoretical world (BCM)

## 3.1 Bootstrap empirical metrics (fast)

We summarize alpha, beta, and simple network metrics by bootstrapping figs. Defaults are small to keep this vignette snappy.


```{r bootstrap-metrics}
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



```

## 3.2 Parameter ranges for BCM calibration

We convert a parameter list into ranges for Latin hypercube exploration. We remove fixed items (interaction matrix/weight).

```{r param-ranges}
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



```

## 3.3 Explore parameter space (short demo vs. full)

**Important:** LHS over thousands of samples × 1000‑fig simulations is expensive. In the vignette we either (A) load a precomputed result, or (B) run a tiny demo search. For full runs, set `RUN_HEAVY <- TRUE` and increase `n_samples`, `n_draws`, `sample_n`, and `num_figs`.

```{r explore-params}
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
```

## 3.4 Simulate the BCM with best parameters

```{r simulate-bcm}
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


```

# 4. Quantify the model–data residual variance (‘Variation Gap’)

We compare Observed vs Simulated by resampling simulated figs to the same sample size and computing identical metrics.

```{r quantify-gap}
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

```
```{r plot-variation-gap-pca, fig.cap="Figure 2: Visualization of the Variation Gap using Principal Component Analysis (PCA). The dense cloud of blue points represents the 'Theoretical Envelope' generated by the BCM. The widely dispersed red points are the empirical observations from the field, clearly falling outside the theoretical boundary."}

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

```

The Variation Gap occurs when the simulated cloud matches means but shows less multivariate dispersion than what is seen (the intrinsic model understates actual variability). This result diagnoses a limit of the present model, not of ecological theory writ large. The magnitude and even the presence of the gap are model-, metric-, and dataset-dependent (e.g., choice of distance, sample size, site heterogeneity). As additional, well-supported processes are incorporated, the gap is expected to shrink.

# 5. Dissect internal mechanisms (knockouts)

We switch off individual modules (e.g., spatial layering, host sanctions) to quantify each mechanism’s contribution relative to the BCM.
We set drop_cancels_emergence = FALSE to retain auditability (dropped figs keep their raw emergence counts). For analysis, we apply a masking policy on a copy of the data (set emergence to 0 for dropped figs), keeping the raw outputs intact.

```{r knockouts}

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

# Summarize metrics for BCM vs. KO
bcm_metrics  <- resample_metrics_from_simulation(sim_df,       n_reps = 50, sample_size = 100, seed = 1)
ko_space   <- resample_metrics_from_simulation(sim_no_space, n_reps = 50, sample_size = 100, seed = 1)

# Loss vs empirical means (smaller = better to observed)
obs_means_row <- as.data.frame(t(colMeans(obs_metrics[ , setdiff(names(obs_metrics), "source") ], na.rm = TRUE)))
bcm_loss      <- calculate_metric_loss(bcm_metrics, obs_means_row)
ko_space_loss       <- calculate_metric_loss(ko_space,  obs_means_row)

list(BCM_total_loss = bcm_loss$total_loss,
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

list(BCM_total_loss = bcm_loss$total_loss,
     KO_no_space_total_loss = ko_space_loss$total_loss,
     KO_prob_by_phase_total_loss = ko_prob_by_phase_loss$total_loss)


```
```{r plot-knockouts, fig.cap="Figure 3: Dissection of internal mechanisms via knockout experiments. The bar height represents the increase in model loss (i.e., deviation from empirical means) when a specific mechanism is removed, relative to the full BCM. A larger value indicates a greater contribution of that mechanism to maintaining the theoretical community structure."}

knockout_results <- data.frame(
  mechanism = c(
    "Full Intrinsic Model (BCM)",
    "KO: Spatial Layering",
    "KO: Phase-specific Egg Success"
  ),
  total_loss = c(
    bcm_loss$total_loss,
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

```

# 6. Host Sanction Threshold

Example — Host-sanction threshold (extreme setting).
Here we illustrate an extreme host-sanction setting. We set the host-sanction threshold to 10% of a fig’s total flowers. If the number of ovules occupied by eggs of ovule-occupying guilds (pollinators + gallers) exceeds 10% of the fig’s flower count, the fig drops. When a fig drops, no fig wasps emerge. This example is for demonstration only; in practice, the threshold should be calibrated to empirical biology.

Notes: (i) Parasitoid eggs are not counted toward the ovule-occupation total because they do not occupy ovules directly. (ii) If your function still uses the legacy argument name host_sanction, use it; otherwise use the corrected host_sanction_threshold. (iii) A new sink theory was incorporated into figsimR_0.2.0.tar.gz, more details: Simon T Segar, Sotiria Boutsi, Daniel Souto-Vilarós, Martin Volf, Derek W Dunn, Astrid Cruaud, Rodrigo A S Pereira, Jean-Yves Rasplus, Finn Kjellberg, The diversity of Ficus, Annals of Botany, 2025;, mcaf280, https://doi.org/10.1093/aob/mcaf280

```{r strict host-sanction threshold }
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
  # new feature: sink metrics
  use_sink_strength = FALSE,
  sink_w_gall          = 1.0,
  sink_w_seed          = 1.5,
  sink_linear_coef     = 1.0,
  sink_min_prop        = 0.20,
  sink_max_prop        = 0.95
  seed                 = 43
)

sim_sanction_10 <- sim_sanction_10$summary

# Quick check: how many figs dropped, and confirm dropped figs have zero emergence
summarize_simulated_metrics(sim_sanction_10, species_list, version_label = "sanction 10%")$drop_summary  

summarize_simulated_metrics(sim_sanction_10, species_list, version_label = "sanction 10%")$stats_table

```

The bigger the loss jump when a module is removed, the more important that module is for reproducing community structure within the intrinsic envelope. We use ‘intrinsic envelope’ purely as a modeling term; it does not imply other processes are stochastic by nature—only that they are not yet modeled here.

# 7. Practical notes for users

* **Speed vs. rigor.** The vignette runs with small `n_draws`, `sample_n`, and `n_samples`. For real analyses, increase these substantially.

* **Reproducibility.** Use `set.seed()` everywhere. If you enable parallelism, document RNG streams.

* **Progress bars.** `progressr::handlers()` may conflict with test runners—keep them off by default in vignettes.

* **Precomputed objects.** For long workflows, ship `.rds` and load them when available (as shown above).

* **Metric sensitivity.** Gap estimates depend on the dispersion metric (e.g., Bray–Curtis vs. Jaccard) and resampling scheme. Report effect sizes with uncertainty and perform sensitivity analyses.

* **Extensibility.** Treat the BCM vignette as a baseline. If a project requires environmental or spatial structure, add those modules and re-estimate the gap.

# 8. Session info

```{r session-info}
sessionInfo()
