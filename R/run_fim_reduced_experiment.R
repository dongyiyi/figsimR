#' Run a Reduced FIM Simulation and Bootstrap Key Metrics
#'
#' This function runs a reduced version of the Full Intrinsic Model (FIM)
#' simulation using configurable ecological mechanisms, and computes
#' bootstrapped summaries of key diversity and network metrics.
#'
#' @param version_label Character. Label to assign to the output (e.g., "NoLayering").
#' @param lhs_param_list List of input parameters (entry, fecundity, etc.).
#' @param species_list List of species to include in network summary.
#' @param num_figs Integer. Number of figs to simulate. Default is 1000.
#' @param host_sanction Numeric. Host sanction strength (0–1).
#' @param use_layering Logical. Whether to enable spatial layering. Default TRUE.
#' @param use_layer_preference Logical. Whether to apply layer preference. Default TRUE.
#' @param use_supplemental_parasitism Logical. Whether to activate supplemental parasitism. Default FALSE.
#' @param use_egg_success_by_phase Logical. Use egg success probabilities by phase. Default TRUE.
#' @param drop_cancels_emergence Logical. Whether fruit drop cancels emergence. Default FALSE.
#' @param seed Integer. Random seed for reproducibility.
#' @param resample_n Integer. Number of bootstrap iterations. Default 500.
#' @param resample_size Integer. Sample size per bootstrap iteration. Default 200.
#'
#' @return A data frame of 6 metrics (plus source label), including:
#' \itemize{
#'   \item Richness
#'   \item Pielou's Evenness
#'   \item Bray-Curtis Dissimilarity
#'   \item Nestedness
#'   \item Connectance
#'   \item Modularity
#' }
#'
#' @seealso \code{\link{simulate_figwasp_community}}, \code{\link{resample_metrics_from_simulation}}, \code{\link{summarize_simulated_metrics}}
#'
#' @export


run_fim_reduced_experiment <- function(
    version_label,
    lhs_param_list,
    species_list,
    num_figs = 1000,
    host_sanction = 0.8,
    use_layering = TRUE,
    use_layer_preference = TRUE,
    use_supplemental_parasitism = FALSE,
    use_egg_success_by_phase = TRUE,
    drop_cancels_emergence = FALSE,
    seed = 42,
    resample_n = 500,
    resample_size = 200
) {

  sim_result <- simulate_figwasp_community(
    num_figs = num_figs,
    fecundity_mean = lhs_param_list$fecundity_mean,
    fecundity_dispersion = lhs_param_list$fecundity_dispersion,
    entry_mu = lhs_param_list$entry_mu,
    entry_size = lhs_param_list$entry_size,
    entry_priority = lhs_param_list$entry_priority,
    species_roles = lhs_param_list$species_roles,
    max_entry_table = lhs_param_list$max_entry_table,
    enable_drop = TRUE,
    drop_cancels_emergence = drop_cancels_emergence,
    entry_distribution = "lognormal",
    interaction_matrix = lhs_param_list$interaction_matrix,
    interaction_weight = 0,
    seed = seed,
    egg_success_prob = NULL,
    egg_success_prob_by_phase = lhs_param_list$egg_success_prob_by_phase,
    host_sanction = host_sanction,
    use_supplemental_parasitism = use_supplemental_parasitism,
    parasitism_prob = lhs_param_list$parasitism_prob,
    record_individual = FALSE,
    use_egg_success_by_phase = use_egg_success_by_phase,
    use_layering = use_layering,
    use_layer_preference = use_layer_preference,
    layer_preference = lhs_param_list$layer_preference,
    p_pollination_per_ovule = 0.98
  )


  sim_summary <- sim_result$summary
  summarized <- summarize_simulated_metrics(sim_summary, species_list, version_label = version_label)


  bootstrap_result <- resample_metrics_from_simulation(sim_summary,
                                                       n_reps = resample_n,
                                                       sample_size = resample_size,
                                                       seed = seed)
  bootstrap_result$Source <- version_label


  reduced <- bootstrap_result[, c(1,4,5,7,8,10,11)]
  colnames(reduced) <- c("Richness", "Pielou's Evenness", "Bray-Curtis Dissimilarity",
                         "Nestedness", "Connectance", "Modularity", "Source")

  return(reduced)
}
