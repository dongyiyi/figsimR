#' Resample and Compute Community Metrics from Simulated Data
#'
#' This function performs repeated subsampling of simulated fig wasp data to assess the variability
#' in community-level metrics. It automatically detects columns containing emergence counts
#' (i.e., columns starting with `"emergence_"`) and calculates diversity or network statistics
#' using a user-defined `calc_all_metrics()` function.
#'
#' @param sim_data A data frame returned by the simulator, typically \code{sim_result$summary}.
#' @param n_reps Number of bootstrap replicates to perform. Default is 500.
#' @param sample_size Number of figs to draw per replicate. Default is 200.
#' @param seed Random seed for reproducibility. Default is 42.
#'
#' @return A data frame with one row per replicate, each row containing calculated community metrics.
#'
#' @details
#' Columns with names starting with \code{"emergence_"} are automatically extracted as wasp emergence data.
#' If the metric computation fails for a given replicate, it will be skipped with a warning.
#'
#' @examples
#' \dontrun{
#' sim_result <- simulate_figwasp_community(num_figs = 1000, ...)
#' summary_df <- sim_result$summary
#' metrics_df <- resample_metrics_from_simulation(summary_df)
#' }
#'
#' @export
resample_metrics_from_simulation <- function(sim_data, n_reps = 500, sample_size = 200, seed = 42) {
  set.seed(seed)


  bee_cols <- grep("^emergence_", colnames(sim_data), value = TRUE)

  if (length(bee_cols) == 0) {
    stop("No emergence_ columns found in simulated data.")
  }

  results_list <- vector("list", n_reps)

  for (i in seq_len(n_reps)) {
    sampled_df <- sim_data[sample(1:nrow(sim_data), size = sample_size, replace = TRUE), ]
    sampled_bees <- sampled_df[, bee_cols, drop = FALSE]

    metrics <- calc_all_metrics(sampled_bees)

    if (is.null(metrics)) {
      warning(paste("Metrics computation failed at iteration", i))
      next
    }

    results_list[[i]] <- metrics
  }

  result_df <- do.call(rbind, results_list)
  rownames(result_df) <- NULL
  return(result_df)
}
