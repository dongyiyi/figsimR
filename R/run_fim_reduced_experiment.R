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
#' @param use_sink_strength Logical. If TRUE, compute sink-strength metrics but do NOT alter legacy
#' drop decision: use_sink_strength = TRUE and enable_drop = FALSE and the columns "is_dropped" and "drop_prob"
#' are not meaningful. The fig is dropped when the value of "sink_prop" falls outside the thresholds defined
#' by sink_min_prop and sink_max_prop; The columns "sink_below_min" and "sink_above_max" are
#' logical: TRUE = drop; FALSE = no drop. Adds columns: sink_strength, sink_prop, sink_below_min,
#' sink_above_max. Sink strength threshold calculation: sink_strength = sink_linear_coef * (sink_w_gall * galled_total + sink_w_seed * seed_count).
#' @param enable_drop Logical. Enable legacy drop (host sanction) mechanism. If flower use proportion exceeds
#' host_sanction, fig may abort. In the output "drop_prob" column: Ranges from 0 (never drop) to 1 (always drop);
#' "is_dropped" column: no drop (0) and drop (1). When sink strength is used (use_sink_strength = TRUE)
#' and enable_drop = FALSE (host sanction is disabled); the columns "is_dropped" and "drop_prob"
#' in the outputs are not meaningful.
#' @param sink_min_prop Numeric in [0,1]. Reference lower bound of sink proportion (default 0.20).
#' @param sink_max_prop Numeric in [0,1]. Reference upper bound of sink proportion (default 0.95).
#' @param sink_controls_drop Logical. If \code{TRUE}, sink-strength thresholds
#'   are used to mark figs as dropped after simulation.
#' @param sink_drop_condition Character. Sink-based drop rule; one of
#'   \code{"outside_range"}, \code{"below_min"}, or \code{"above_max"}.
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
    enable_drop = TRUE,
    use_egg_success_by_phase = TRUE,
    drop_cancels_emergence = FALSE,
    use_sink_strength = FALSE,
    sink_min_prop = 0.20,
    sink_max_prop = 0.95,
    sink_controls_drop = FALSE,
    sink_drop_condition = c("outside_range", "below_min", "above_max"),
    seed = 42,
    resample_n = 500,
    resample_size = 200
) {

  sink_drop_condition <- match.arg(sink_drop_condition)

  sim_result <- simulate_figwasp_community(
    num_figs = num_figs,

    fecundity_mean = lhs_param_list$fecundity_mean,
    fecundity_dispersion = lhs_param_list$fecundity_dispersion,
    entry_mu = lhs_param_list$entry_mu,
    entry_size = lhs_param_list$entry_size,
    entry_priority = lhs_param_list$entry_priority,
    species_roles = lhs_param_list$species_roles,
    max_entry_table = lhs_param_list$max_entry_table,

    # important: do not hard-code this as TRUE
    enable_drop = enable_drop,
    drop_cancels_emergence = drop_cancels_emergence,
    host_sanction = host_sanction,

    entry_distribution = "lognormal",
    interaction_matrix = lhs_param_list$interaction_matrix,
    interaction_weight = 0,
    seed = seed,

    egg_success_prob = NULL,
    egg_success_prob_by_phase = lhs_param_list$egg_success_prob_by_phase,
    use_egg_success_by_phase = use_egg_success_by_phase,

    use_supplemental_parasitism = use_supplemental_parasitism,
    parasitism_prob = lhs_param_list$parasitism_prob,

    record_individual = FALSE,
    use_layering = use_layering,
    use_layer_preference = use_layer_preference,
    layer_preference = lhs_param_list$layer_preference,
    p_pollination_per_ovule = 0.98,

    # sink-strength module
    use_sink_strength = use_sink_strength,
    sink_min_prop = sink_min_prop,
    sink_max_prop = sink_max_prop
  )

  sim_summary <- sim_result$summary

  # Optional sink-strength-based drop assignment.
  # This is applied after simulation because use_sink_strength is diagnostic-only
  # in the current simulate_figwasp_community() implementation.
  if (use_sink_strength && sink_controls_drop) {

    required_sink_cols <- c("sink_below_min", "sink_above_max")
    missing_sink_cols <- setdiff(required_sink_cols, names(sim_summary))

    if (length(missing_sink_cols) > 0) {
      stop(
        "sink_controls_drop = TRUE, but these sink columns are missing: ",
        paste(missing_sink_cols, collapse = ", ")
      )
    }

    sink_drop_flag <- switch(
      sink_drop_condition,
      below_min = sim_summary$sink_below_min,
      above_max = sim_summary$sink_above_max,
      outside_range = sim_summary$sink_below_min | sim_summary$sink_above_max
    )

    sim_summary$is_dropped <- as.integer(sink_drop_flag %in% TRUE)

    # This is deterministic threshold-based drop, not probabilistic host sanction.
    sim_summary$drop_prob <- NA_real_
    sim_summary$drop_mechanism <- ifelse(
      sim_summary$is_dropped == 1,
      "sink_strength",
      "none"
    )
  }

  # Keep full drop summary before filtering
  if (!"is_dropped" %in% names(sim_summary)) {
    sim_summary$is_dropped <- 0L
  }

  drop_summary_full <- table(
    factor(sim_summary$is_dropped, levels = c(0, 1))
  )

  # Remove dropped figs before metric calculation
  sim_summary_for_metrics <- sim_summary[
    is.na(sim_summary$is_dropped) | sim_summary$is_dropped == 0,
    ,
    drop = FALSE
  ]

  message(
    version_label, ": removed ",
    sum(sim_summary$is_dropped == 1, na.rm = TRUE),
    " dropped figs; retained ",
    nrow(sim_summary_for_metrics),
    " figs for metric calculation."
  )

  if (nrow(sim_summary_for_metrics) == 0) {
    stop(version_label, ": no figs retained after dropping.")
  }

  summarized <- summarize_simulated_metrics(
    sim_summary_for_metrics,
    species_list,
    version_label = version_label
  )

  # Preserve the full drop table, because summarize_simulated_metrics()
  # now only sees retained figs.
  summarized$drop_summary_full <- drop_summary_full

  bootstrap_result <- resample_metrics_from_simulation(
    sim_summary_for_metrics,
    n_reps = resample_n,
    sample_size = resample_size,
    seed = seed
  )

  bootstrap_result$Source <- version_label

  reduced <- bootstrap_result[, c(1, 4, 5, 7, 8, 10, 11)]

  colnames(reduced) <- c(
    "Richness",
    "Pielou's Evenness",
    "Bray-Curtis Dissimilarity",
    "Nestedness",
    "Connectance",
    "Modularity",
    "Source"
  )

  return(reduced)
}
