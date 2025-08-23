#' Run Simulations and Bootstrap Community Metrics Across Interaction Weights
#'
#' This function automates the process of running the fig wasp community simulator
#' with multiple values of \code{interaction_weight}, and performs bootstrap resampling
#' on each result to assess variability in community metrics.
#'
#' @param weights A numeric vector of \code{interaction_weight} values to evaluate. Default: c(0.1, 0.3, 0.5, 0.7, 0.9).
#' @param n_reps Number of bootstrap replicates per weight. Default: 500.
#' @param sample_size Number of figs to sample per replicate. Default: 200.
#' @param seed Random seed for reproducibility. Default: 42.
#'
#' @return A named list of data frames. Each element corresponds to one interaction weight and contains
#' bootstrap resampled community metrics computed from the simulation.
#'
#' @details
#' This function requires that global objects such as \code{lhs_best_param_list},
#' \code{observed_interaction_matrix}, and \code{simulate_figwasp_community()}
#' be pre-defined. The function uses \code{resample_metrics_from_simulation()} internally
#' to compute alpha diversity or other community structure metrics for each simulation run.
#'
#' @examples
#' \dontrun{
#' results_list <- run_sim_and_bootstrap_multiple_weights(
#'   weights = seq(0.1, 0.9, by = 0.2),
#'   n_reps = 300,
#'   sample_size = 150
#' )
#' }
#'
#' @export

run_sim_and_bootstrap_multiple_weights <- function(weights = c(0.1, 0.3, 0.5, 0.7, 0.9),
                                                   n_reps = 500,
                                                   sample_size = 200,
                                                   seed = 42) {

  results <- list()

  for (w in weights) {
    message("Running simulation for interaction_weight = ", w)

    sim_data <- simulate_figwasp_community(
      num_figs = 1000,
      fecundity_mean = lhs_best_param_list$fecundity_mean,
      fecundity_dispersion = lhs_best_param_list$fecundity_dispersion,
      entry_mu = lhs_best_param_list$entry_mu,
      entry_size = lhs_best_param_list$entry_size,
      entry_priority = lhs_best_param_list$entry_priority,
      species_roles = lhs_best_param_list$species_roles,
      max_entry_table = lhs_best_param_list$max_entry_table,
      enable_drop = TRUE,
      drop_cancels_emergence = FALSE,
      entry_distribution = "lognormal",
      interaction_matrix = observed_interaction_matrix,
      interaction_weight = w,
      seed = seed,
      egg_success_prob = NULL,
      egg_success_prob_by_phase = lhs_best_param_list$egg_success_prob_by_phase,
      host_sanction = 0.8,
      use_supplemental_parasitism = FALSE,
      parasitism_prob = lhs_best_param_list$parasitism_prob,
      record_individual = FALSE,
      use_layering = TRUE,
      layer_preference = lhs_best_param_list$layer_preference,
      p_pollination_per_ovule = 0.98
    )$summary


    boot_metrics <- resample_metrics_from_simulation(sim_data, n_reps = n_reps, sample_size = sample_size, seed = seed)


    name <- paste0("fim_bootstrap_metrics_", w)
    results[[name]] <- boot_metrics
  }

  return(results)
}
