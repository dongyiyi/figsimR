#' Perform Kruskal-Wallis and Dunn's Post Hoc Test for Multiple Metrics
#'
#' This function applies the Kruskal-Wallis test to each unique metric in the input data
#' and conducts a post hoc Dunn's test with Benjamini-Hochberg adjustment if the
#' Kruskal-Wallis test is significant (p < 0.05).
#'
#' @param df A long-format data frame containing one row per metric observation.
#' @param metric_col Name of the column that identifies the metric (default: "Metric").
#' @param group_col Name of the column that contains grouping information (e.g., model type).
#' @param value_col Name of the column containing the numeric values to compare.
#'
#' @return A data frame with one row per comparison, including:
#' \describe{
#'   \item{\code{Metric}}{The metric name.}
#'   \item{\code{KW_p}}{The Kruskal-Wallis test p-value.}
#'   \item{\code{Comparison}}{Group comparison name (e.g., "A - B").}
#'   \item{\code{Z}}{The z-score from Dunn's test.}
#'   \item{\code{P.unadj}}{Unadjusted p-value from Dunn's test.}
#'   \item{\code{P.adj}}{BH-adjusted p-value.}
#' }
#'
#' If the Kruskal-Wallis test is not significant for a metric, only the KW p-value will be shown, and other fields will be \code{NA}.
#'
#' @examples
#' \dontrun{
#' long_df <- data.frame(
#'   Metric = rep(c("Richness", "Shannon"), each = 30),
#'   Source = rep(c("Observed", "Simulated1", "Simulated2"), 20),
#'   Value = rnorm(60)
#' )
#' result <- run_kruskal_dunn(long_df)
#' print(result)
#' }
#'
#' @import dplyr
#' @importFrom FSA dunnTest
#' @export

run_kruskal_dunn <- function(df, metric_col = "Metric", group_col = "Source", value_col = "Value") {

  metrics <- unique(df[[metric_col]])


  results <- lapply(metrics, function(metric_name) {
    sub_df <- df %>% filter(.data[[metric_col]] == metric_name)


    kw <- kruskal.test(as.formula(paste(value_col, "~", group_col)), data = sub_df)


    if (kw$p.value < 0.05) {
      dunn <- dunnTest(as.formula(paste(value_col, "~", group_col)),
                       data = sub_df, method = "bh")$res

      dunn$Metric <- metric_name
      dunn$KW_p <- round(kw$p.value, 4)
      return(dunn)
    } else {
      return(data.frame(
        Metric = metric_name,
        KW_p = round(kw$p.value, 4),
        Comparison = NA,
        Z = NA,
        P.unadj = NA,
        P.adj = NA
      ))
    }
  })

  dplyr::bind_rows(results)
}
