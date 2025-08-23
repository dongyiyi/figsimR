#' Plot Metric Distributions with Wilcoxon Test Significance
#'
#' This function compares the distributions of community-level metrics between observed
#' and simulated datasets. It visualizes the density plots for each numeric metric side-by-side
#' and annotates each with the Wilcoxon test p-value and corresponding significance level.
#'
#' @param obs_metrics A data frame of observed community metrics (e.g., alpha, beta, network structure).
#'   Must include only numeric columns to be tested and plotted.
#' @param sim_metrics A data frame of simulated metrics from the same set of metrics as \code{obs_metrics}.
#'   The columns should have the same names and types as in \code{obs_metrics}.
#' @param ncol Integer. Number of columns to arrange plots in the output grid. Default is 3.
#'
#' @return A grid of \code{ggplot2} density plots, one for each shared numeric metric.
#'   Each plot includes a Wilcoxon rank-sum test result between the observed and simulated values,
#'   reported as p-value and a significance label (*** < 0.001, ** < 0.01, * < 0.05, ns otherwise).
#'
#' @details
#' This function is useful for model evaluation, allowing users to assess whether the distributions
#' of key diversity or network metrics produced by a simulation (e.g., FIM) differ significantly
#' from those observed in empirical data. It uses nonparametric Wilcoxon rank-sum tests for comparison.
#'
#' Metrics that are not numeric or are not shared between \code{obs_metrics} and \code{sim_metrics}
#' will be ignored.
#'
#' @examples
#' \dontrun{
#' plot_distribution_metrics_grid_pvalue(observed_summary_df, simulated_summary_df, ncol = 4)
#' }
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows mutate
#' @importFrom gridExtra grid.arrange
#' @export
#'
plot_distribution_metrics_grid_pvalue <- function(obs_metrics, sim_metrics, ncol = 3) {




  obs_metrics <- obs_metrics %>% dplyr::mutate(source = "Observed")
  sim_metrics <- sim_metrics %>% dplyr::mutate(source = "Simulated")

  shared_cols <- intersect(
    colnames(obs_metrics)[sapply(obs_metrics, is.numeric)],
    colnames(sim_metrics)[sapply(sim_metrics, is.numeric)]
  )

  all_metrics <- dplyr::bind_rows(obs_metrics, sim_metrics)


  get_significance_label <- function(p) {
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("ns")
  }

  plot_list <- lapply(shared_cols, function(metric) {

    obs_values <- obs_metrics[[metric]]
    sim_values <- sim_metrics[[metric]]


    test <- tryCatch(
      wilcox.test(obs_values, sim_values),
      error = function(e) list(p.value = NA)
    )
    p_value <- test$p.value
    signif_label <- get_significance_label(p_value)


    plot_title <- paste0(metric, ", p = ",
                         format.pval(p_value, digits = 3), ", ", signif_label, ")")

    ggplot2::ggplot(all_metrics, ggplot2::aes(x = !!sym(metric), fill = source)) +
      geom_density(alpha = 0.5) +
      labs(
        title = plot_title,
        x = metric,
        y = "Density"
      ) +
      scale_fill_manual(
        labels = c("Observed Data", "Simulated FIM"),
        values = c("Observed" = "#fb8072", "Simulated" = "steelblue")
      )+
      theme_classic() + ggplot2::theme(legend.position = "none")
  })

  grid.arrange(grobs = plot_list, ncol = ncol)
}
