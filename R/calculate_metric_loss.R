#' Calculate Total Loss Between Simulated and Observed Metrics
#'
#' This function compares the mean values of community metrics from simulated data
#' to a single row of observed values, computing the squared error loss for each metric
#' and returning both the per-metric loss and the total loss.
#'
#' @param simulated_df A data frame of simulated metric results, typically generated
#' from bootstrap resampling (e.g., via \code{resample_metrics_from_simulation()}).
#' @param observed_df A one-row data frame containing observed metric values.
#' Column names should match those in \code{simulated_df}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{metric_table}}{A data frame showing the simulated mean, observed value,
#'   and squared difference for each metric.}
#'   \item{\code{total_loss}}{The total loss, defined as the sum of squared differences
#'   across all shared metrics.}
#' }
#'
#' @examples
#' \dontrun{
#' sim_metrics <- resample_metrics_from_simulation(sim_data)
#' loss_result <- calculate_metric_loss(sim_metrics, observed_result_df)
#' print(loss_result$total_loss)
#' print(loss_result$metric_table)
#' }
#'
#' @export


calculate_metric_loss <- function(simulated_df, observed_df) {

  simulated_means <- colMeans(simulated_df, na.rm = TRUE)


  observed_values <- as.numeric(observed_df[1, ])
  names(observed_values) <- colnames(observed_df)


  common_metrics <- intersect(names(simulated_means), names(observed_values))
  loss_components <- (simulated_means[common_metrics] - observed_values[common_metrics])^2


  result <- data.frame(
    Metric = common_metrics,
    Simulated_Mean = round(simulated_means[common_metrics], 4),
    Observed_Value = round(observed_values[common_metrics], 4),
    Squared_Diff = round(loss_components, 4)
  )

  total_loss <- sum(loss_components, na.rm = TRUE)


  list(
    metric_table = result,
    total_loss = round(total_loss, 4)
  )
}
