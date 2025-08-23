#' Bootstrap Resampling of Observed Community Metrics
#'
#' This function resamples the observed fig wasp community dataset multiple times
#' to estimate the distribution of community metrics, such as diversity and network structure.
#' It is typically used to compare observed patterns against simulated outputs.
#'
#' @param observed_df A data.frame containing fig-level community observations,
#'   with one row per fig and one column per wasp species (typically emergence counts).
#'
#' @param wasp_cols A character vector indicating which columns in \code{observed_df}
#'   correspond to wasp species. These columns will be used to compute metrics.
#'
#' @param n_draws Integer. Number of bootstrap replicates to generate. Default = 500.
#'
#' @param sample_n Integer. Number of figs to sample per bootstrap iteration. Default = 200.
#'
#' @param seed Optional integer. If provided, sets the random seed for reproducibility.
#'
#' @param calc_func A function to calculate community metrics. Must accept a data.frame
#'   of species abundance (rows = figs, columns = species). Default = \code{\link{calc_all_metrics}}.
#'
#' @return A data.frame containing \code{n_draws} rows, each row representing a bootstrap
#'   replicate. Each row includes all metrics returned by \code{calc_func}, along with
#'   a column \code{source = "Observed"}.
#'
#' @details
#' In each iteration, the function randomly samples \code{sample_n} rows (with replacement)
#' from the observed dataset, calculates community metrics using \code{calc_func},
#' and stores the results.
#'
#' This resampling approach estimates the empirical distribution of community metrics
#' under the observed data, providing a null reference for comparison with simulation outputs.
#'
#' @examples
#' # Example with toy data
#' set.seed(123)
#' df <- data.frame(
#'   fig_id = 1:300,
#'   A = rpois(300, 5),
#'   B = rpois(300, 3),
#'   C = rpois(300, 1)
#' )
#'
#' observed_metrics <- resample_observed_all_metrics(
#'   observed_df = df,
#'   wasp_cols = c("A", "B", "C"),
#'   n_draws = 100,
#'   sample_n = 50,
#'   seed = 42
#' )
#'
#' head(observed_metrics)
#'
#' @importFrom dplyr bind_rows
#' @export
resample_observed_all_metrics <- function(
    observed_df,
    wasp_cols,
    n_draws = 500,
    sample_n = 200,
    seed = NULL,
    calc_func = calc_all_metrics
) {
  if (!is.null(seed)) set.seed(seed)

  df <- observed_df
  if (!all(wasp_cols %in% colnames(df))) {
    stop("Some wasp_cols are not present in the observed data.")
  }

  metrics_list <- vector("list", n_draws)

  for (i in 1:n_draws) {
    sampled_df <- df[sample(nrow(df), sample_n), wasp_cols, drop = FALSE]
    metrics_list[[i]] <- calc_func(sampled_df)
  }

  result_df <- dplyr::bind_rows(metrics_list)
  result_df$source <- "Observed"
  return(result_df)
}
