#' Plot Metric Distributions with KS Test and Shared Legend
#'
#' This function compares the distributions of community-level metrics between observed and simulated datasets
#' using kernel density plots. Each metric is tested with a Kolmogorov-Smirnov (KS) test, and the resulting p-value
#' is displayed as the title of each plot. A shared legend is displayed above all subplots.
#'
#' @param obs_metrics A data frame containing observed metric values. Only numeric columns will be used.
#' @param sim_metrics A data frame containing simulated metric values. Must have the same structure as \code{obs_metrics}.
#' @param ncol Number of columns in the output plot grid. Default is 3.
#'
#' @return A grid of density plots (one per metric) comparing observed and simulated distributions, with
#' a shared legend on top. Each plot shows a KS test p-value between observed and simulated distributions.
#'
#' @details
#' This function helps evaluate model performance by comparing the distribution of each community-level
#' metric (e.g., alpha diversity, richness, network structure) between observed and simulated datasets.
#' A Kolmogorov-Smirnov test is performed for each metric, and its p-value is shown in the plot title.
#'
#' The shared legend distinguishes between observed and simulated data and is extracted from the first plot.
#' All plots are styled consistently for integration into figure panels or reports.
#'
#' @examples
#' \dontrun{
#' plot_distribution_metrics_grid_pvalue_legend(observed_df, simulated_df, ncol = 4)
#' }
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows mutate
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid unit
#' @export


plot_distribution_metrics_grid_pvalue_legend <- function(obs_metrics, sim_metrics, ncol = 3) {





  obs_metrics <- obs_metrics %>% dplyr::mutate(Source = "Observed")
  sim_metrics <- sim_metrics %>% dplyr::mutate(Source = "Simulated")

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
      ks.test(obs_values, sim_values),
      error = function(e) list(p.value = NA)
    )
    p_value <- test$p.value
    signif_label <- get_significance_label(p_value)

    plot_title <- paste0("", "p = ",
                         format.pval(p_value, digits = 3))

    ggplot2::ggplot(all_metrics, ggplot2::aes(x = !!sym(metric), fill = Source)) +
      geom_density(alpha = 0.5) +
      labs(title = plot_title, x = metric, y = "Density") +
      scale_fill_manual(
        labels = c("Observed Data", "Simulated FIM"),
        values = c("Observed" = "#fb8072", "Simulated" = "steelblue")
      ) +
      theme_classic(base_size = 6, base_family = "Arial") +
      ggplot2::theme(

        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5)
      ) +ggplot2::theme(legend.position = "none")
  })


  legend_plot <- ggplot2::ggplot(all_metrics, ggplot2::aes(x = !!sym(shared_cols[1]), fill = Source)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(
      labels = c("Observed Data", "Simulated FIM"),
      values = c("Observed" = "#fb8072", "Simulated" = "steelblue")
    ) +
    ggplot2::theme(legend.position = "top",
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          legend.key.size = unit(3, "mm"))


  get_legend <- function(plot) {
    g <- ggplotGrob(plot)
    legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
    g$grobs[[legend_index]]
  }
  shared_legend <- get_legend(legend_plot)


  combined_plots <- arrangeGrob(grobs = plot_list, ncol = ncol)
  grid.arrange(shared_legend, combined_plots, ncol = 1, heights = c(1, 10))
}
