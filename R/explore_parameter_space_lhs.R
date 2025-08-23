#' Explore Parameter Space Using Latin Hypercube Sampling
#'
#' This function systematically explores the parameter space of the fig wasp community simulator
#' using Latin Hypercube Sampling (LHS), evaluates simulation output against observed metrics,
#' and returns the top parameter combinations that minimize loss.
#'
#' @param new_param_ranges A named list of parameter ranges, where each element is a numeric vector of length 2
#'   specifying the lower and upper bounds for a parameter. Names should follow the format: \code{type_species}
#'   (e.g., \code{entry_mu_PollA}) or \code{species_phaseX} for phase-specific values (e.g., \code{Apocrypta_sp_phase1}).
#' @param observed_summary A data.frame of observed community metrics (e.g., richness, Shannon, Bray-Curtis),
#'   typically from \code{\link{resample_observed_all_metrics}}.
#' @param num_figs Integer. Number of figs to simulate in each replicate. Default = 1000.
#' @param n_draws Integer. Number of bootstrap replicates per parameter combination. Default = 500.
#' @param sample_n Integer. Number of figs sampled per replicate for calculating community metrics. Default = 200.
#' @param wasp_cols Character vector. Names of columns in simulation output corresponding to fig wasp species.
#' @param parameter_list_template A full template of simulation parameters (typically based on
#'   \code{parameter_list_default}), used to insert sampled parameter values.
#' @param n_samples Integer. Number of LHS parameter combinations to evaluate. Default = 5000.
#' @param top_k Integer. Number of top-performing parameter sets (with lowest loss) to return. Default = 5.
#' @param n_cores Integer. Number of parallel cores to use. Default = 4.
#'
#' @return A data.frame of the top \code{top_k} parameter combinations, sorted by loss,
#'   with one row per combination and columns for each parameter and corresponding loss value.
#'
#' @details
#' This function generates \code{n_samples} unique parameter sets using LHS within the user-defined bounds,
#' runs a stochastic simulation with each set, and computes a loss value based on the squared difference
#' between simulated and observed mean metrics. Simulation is parallelized for efficiency.
#'
#' The \code{parameter_list_template} must include all required simulator inputs. Sampled values are inserted
#' by name-matching into their corresponding list elements. Phase-specific parameters (e.g., egg success probabilities
#' by phase) must be named accordingly.
#'
#' @examples
#' \dontrun{
#' best_params <- explore_parameter_space_lhs(
#'   new_param_ranges = list(
#'     fecundity_mean_PollA = c(5, 20),
#'     entry_mu_PollA = c(0.1, 2.0),
#'     Apocrypta_sp_phase1 = c(0.05, 0.3)
#'   ),
#'   observed_summary = resample_observed_all_metrics(...),
#'   num_figs = 1000,
#'   n_draws = 100,
#'   sample_n = 200,
#'   wasp_cols = c("PollA", "GallerA", "ParaA"),
#'   parameter_list_template = parameter_list_default,
#'   n_samples = 1000,
#'   top_k = 3,
#'   n_cores = 4
#' )
#' }
#'
#' @importFrom lhs randomLHS
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor with_progress
#' @importFrom dplyr arrange
#' @export

explore_parameter_space_lhs<- function(
    new_param_ranges,
    observed_summary,
    num_figs = 1000,
    n_draws = 500,
    sample_n = 200,
    wasp_cols,
    parameter_list_template,
    n_samples = 5000,
    top_k = 5,
    n_cores = 4
) {




  options(future.rng.onMisuse = "ignore")


  plan(multisession, workers = n_cores)
  #handlers(global = TRUE)


  param_names <- names(new_param_ranges)
  lower <- sapply(new_param_ranges, `[`, 1)
  upper <- sapply(new_param_ranges, `[`, 2)


  lhs_samples <- lhs::randomLHS(n_samples, length(param_names))
  sampled_params <- as.data.frame(t(t(lhs_samples) * (upper - lower) + lower))
  colnames(sampled_params) <- param_names


  loss_fn <- function(param_vec_named) {
    param_list <- parameter_list_template
    for (p in names(param_vec_named)) {
      val <- param_vec_named[[p]]
      if (grepl("_phase", p)) {
        parts <- strsplit(p, "_phase")[[1]]
        sp <- parts[1]
        ph <- gsub(".*phase", "", p)
        param_list$egg_success_prob_by_phase[[sp]][[paste0("phase", ph)]] <- val
      } else {
        parts <- strsplit(p, "_")[[1]]
        param_type <- parts[1]
        sp <- paste(parts[-1], collapse = "_")
        param_list[[param_type]][[sp]] <- val
      }
    }

    sim_df <- tryCatch({
      resample_simulated_metrics(
        n_draws = n_draws,
        sample_n = sample_n,
        num_figs = num_figs,
        wasp_cols = wasp_cols,
        fecundity_mean = param_list$fecundity_mean,
        fecundity_dispersion = param_list$fecundity_dispersion,
        entry_mu = param_list$entry_mu,
        entry_size = param_list$entry_size,
        entry_priority = param_list$entry_priority,
        species_roles = param_list$species_roles,
        interaction_matrix = param_list$interaction_matrix,
        interaction_weight = param_list$interaction_weight,
        egg_success_prob = param_list$egg_success_prob,
        egg_success_prob_by_phase = param_list$egg_success_prob_by_phase,
        parasitism_prob = NULL,
        layer_preference = param_list$layer_preference,
        enable_drop = TRUE,
        use_supplemental_parasitism = FALSE,
        drop_cancels_emergence = FALSE,
        use_layering = TRUE
      )
    }, error = function(e) return(NULL))

    if (is.null(sim_df)) return(1e6)

    sim_means <- colMeans(sim_df[, colnames(observed_summary)], na.rm = TRUE)
    obs_means <- colMeans(observed_summary[, colnames(observed_summary)], na.rm = TRUE)
    loss <- sum((sim_means - obs_means)^2)

    return(loss)
  }


  progressr::with_progress({
    p <- progressr::progressor(steps = n_samples)
    losses <- future_lapply(1:n_samples, future.seed = TRUE, FUN = function(i) {
      p()
      loss_fn(sampled_params[i, ])
    })
  })


  result_df <- sampled_params
  result_df$loss <- unlist(losses)
  result_df <- result_df %>% arrange(loss)

  return(result_df[1:top_k, ])
}
