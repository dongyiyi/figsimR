#' Summarize Simulated Metrics from Ficus Wasp Community Output
#'
#' This function computes summary statistics from the output of the
#' `simulate_figwasp_community()` function, including per-species
#' entry, egg, and emergence means, dropped fruit summary, fig-level richness,
#' and species presence proportions by richness class.
#'
#' A heatmap is plotted to visualize how each species' presence frequency varies
#' across figs with different species richness levels.
#'
#' @param all_metrics Data frame. The `$summary` output from `simulate_figwasp_community()`.
#' @param species_list Character vector of species names (used to extract relevant columns).
#' @param version_label Optional string for labeling plots. Default is "v".
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{stats_table}{Data frame with Entry, Egg, and Emergence means per species}
#'   \item{drop_summary}{Table summarizing dropped vs non-dropped figs}
#'   \item{richness_vector}{Integer vector of species richness per fig}
#'   \item{species_by_richness}{Data frame with species presence proportions by richness}
#'   \item{heatmap}{The richness–species presence ggplot object (printed automatically)}
#' }
#'
#' @examples
#' \dontrun{
#' result <- simulate_figwasp_community(...)
#' summary <- summarize_simulated_metrics(result$summary, species_list)
#' }
#'
#' @export


summarize_simulated_metrics <- function(all_metrics, species_list, version_label = "v") {



  simulated_stats <- data.frame(Species = species_list)

  simulated_stats$Entry_Mean <- sapply(species_list, function(sp) {
    mean(all_metrics[[paste0("entry_", sp)]])
  })

  simulated_stats$Eggs_Mean <- sapply(species_list, function(sp) {
    mean(all_metrics[[paste0("eggs_", sp)]])
  })

  simulated_stats$Emergence_Mean <- sapply(species_list, function(sp) {
    mean(all_metrics[[paste0("emergence_", sp)]])
  })


  drop_table <- table(all_metrics$is_dropped)


  emergence_cols <- grep("^emergence_", colnames(all_metrics), value = TRUE)
  fig_richness <- rowSums(all_metrics[, emergence_cols] > 0)

  hist(fig_richness,
       breaks = seq(-0.5, max(fig_richness) + 0.5, 1),
       main = paste("Number of Emerged Species per Fig (", version_label, ")", sep = ""),
       xlab = "Species Richness per Fig",
       col = "steelblue")


  df_long <- all_metrics[, c(emergence_cols)]
  df_long$FigID <- seq_len(nrow(df_long))
  df_long$richness <- fig_richness

  df_long_tidy <- df_long %>%
    pivot_longer(
      cols = all_of(emergence_cols),
      names_to = "Species",
      values_to = "Count"
    ) %>%
    dplyr::mutate(
      Species = gsub("emergence_", "", Species),
      Present = Count > 0
    )


  richness_species_prop <- df_long_tidy %>%
    group_by(richness, Species) %>%
    summarise(
      N_present = sum(Present),
      Total_figs = n_distinct(FigID),
      Proportion = N_present / Total_figs,
      .groups = "drop"
    )


  heatmap_plot <- ggplot2::ggplot(richness_species_prop, ggplot2::aes(x = factor(richness), y = Species, fill = Proportion)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(x = "Richness", y = "Species",
         title = paste("Species Presence Proportion by Richness Level (", version_label, ")", sep = "")) +
    theme_minimal()

  print(heatmap_plot)

  return(list(
    stats_table = simulated_stats,
    drop_summary = drop_table,
    richness_vector = fig_richness,
    species_by_richness = richness_species_prop,
    heatmap = heatmap_plot
  ))
}
