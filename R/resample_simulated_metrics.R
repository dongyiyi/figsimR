#' Bootstrap Resampling of Simulated Community Metrics
#'
#' This function performs multiple independent runs of a fig wasp community simulation model
#' and computes community metrics for a fixed number of sampled figs from each run.
#' It is used to generate the expected distribution of metrics under a given simulation configuration.
#'
#' @param n_draws Integer. Number of simulation replicates to generate. Default = 500.
#' @param sample_n Integer. Number of figs to sample per replicate for metric calculation. Default = 200.
#' @param num_figs Integer. Number of figs to simulate per call to the simulator. Default = 1000.
#' @param simulator_func A function that simulates fig wasp communities and returns
#'   a list with at least a \code{summary} data.frame including wasp emergence and an \code{is_dropped} column.
#'   Default = \code{\link{simulate_figwasp_community}}.
#' @param calc_func A function to calculate community metrics from a sampled data.frame.
#'   Must accept a data.frame of species abundance (rows = figs, columns = species).
#'   Default = \code{\link{calc_all_metrics}}.
#' @param wasp_cols A character vector indicating which columns in the simulation output
#'   correspond to wasp species.
#' @param ... Additional arguments passed to \code{simulator_func}, such as parameter lists or switches.
#'
#' @return A data.frame where each row represents a replicate and includes all community metrics
#'   returned by \code{calc_func}. A column \code{source = "Simulated"} is added to distinguish simulated outputs.
#'
#' @details
#' In each replicate, the simulation function is called with \code{num_figs} figs. Only figs
#' where \code{is_dropped == 0} are retained, and a sample of size \code{sample_n} is randomly drawn
#' for metric calculation. If a replicate contains fewer than \code{sample_n} valid figs,
#' that iteration is skipped.
#'
#' This function is commonly used to generate the expected distribution of metrics
#' under the Full Intrinsic Model (FIM) or other simulation variants, for comparison with observed data.
#'
#' @examples
#' # Example (requires a working simulator and parameter inputs)
#' # Assuming wasp_cols is defined and simulate_figwasp_community() is available
#' \dontrun{
#' sim_metrics <- resample_simulated_metrics(
#'   n_draws = 100,
#'   sample_n = 100,
#'   num_figs = 500,
#'   simulator_func = simulate_figwasp_community,
#'   calc_func = calc_all_metrics,
#'   wasp_cols = c("Ceratosolen_sp", "Sycophaga_testacea", "Apocrypta_sp"),
#'   fecundity_mean = parameter_list_default$fecundity_mean,
#'   entry_mu = parameter_list_default$entry_mu,
#'   entry_size = parameter_list_default$entry_size,
#'   species_roles = parameter_list_default$species_roles,
#'   max_entry_table = parameter_list_default$max_entry_table
#' )
#' head(sim_metrics)
#' }
#'
#' @importFrom dplyr bind_rows
#' @export
resample_simulated_metrics <- function(
    n_draws = 500,
    sample_n = 200,
    num_figs = 1000,
    simulator_func = simulate_figwasp_community,
    calc_func = calc_all_metrics,
    wasp_cols,
    ...
) {

  if (!is.function(simulator_func)) stop("simulator_func must be a function.")
  if (!is.function(calc_func)) stop("calc_func must be a function.")

  metrics_list <- vector("list", n_draws)

  for (i in 1:n_draws) {

    sim_result <- simulator_func(num_figs = num_figs, ...)
    df <- sim_result$summary
    df_valid <- df[df$is_dropped == 0, ]


    if (nrow(df_valid) < sample_n) next


    sampled_df <- df_valid[sample(nrow(df_valid), sample_n), wasp_cols, drop = FALSE]


    metrics_list[[i]] <- calc_func(sampled_df)
  }


  result_df <- dplyr::bind_rows(metrics_list)
  result_df$source <- "Simulated"
  return(result_df)
}
