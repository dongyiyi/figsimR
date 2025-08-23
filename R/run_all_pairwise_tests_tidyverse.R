#' Perform Pairwise Dunn's Tests for Multiple Metrics (Tidyverse Version)
#'
#' This function applies Dunn's test (with Benjamini-Hochberg adjustment)
#' to all pairwise comparisons of groups for each unique metric.
#' It uses a tidyverse pipeline to nest and map over metrics.
#'
#' @param long_df A long-format data frame containing one row per observation.
#' @param value_col Column name for the numeric values to compare. Default is \code{"Value"}.
#' @param group_col Column name representing the groupings (e.g., "Source" or "Treatment"). Default is \code{"Source"}.
#' @param metric_col Column name indicating the metric type. Default is \code{"Metric"}.
#'
#' @return A tidy data frame with the following columns:
#' \describe{
#'   \item{\code{Metric}}{Metric name.}
#'   \item{\code{Comparison}}{Pairwise group comparison (e.g., "Simulated - Observed").}
#'   \item{\code{Z}}{Z statistic from Dunn's test.}
#'   \item{\code{P.unadj}}{Unadjusted p-value.}
#'   \item{\code{P.adj}}{Benjamini-Hochberg adjusted p-value.}
#' }
#'
#' @details This function is a shortcut for performing multiple group comparisons across several ecological or network metrics, assuming Kruskal-Wallis significance. Dunn’s test is performed using the \code{FSA::dunnTest()} function.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Metric = rep(c("Richness", "Shannon"), each = 30),
#'   Source = rep(c("Observed", "Sim1", "Sim2"), times = 20),
#'   Value = rnorm(60)
#' )
#' result <- run_all_pairwise_tests_tidyverse(df)
#' print(result)
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom FSA dunnTest
#' @export



run_all_pairwise_tests_tidyverse <- function(long_df,
                                             value_col = "Value",
                                             group_col = "Source",
                                             metric_col = "Metric") {


  long_df %>%
    group_by(.data[[metric_col]]) %>%
    nest() %>%
    dplyr::mutate(

      test_results = map(data, ~ FSA::dunnTest(as.formula(paste(value_col, "~", group_col)),
                                               data = .x,
                                               method = "bh")$res)
    ) %>%
    select(.data[[metric_col]], test_results) %>%
    unnest(test_results) %>%
    ungroup()
}
