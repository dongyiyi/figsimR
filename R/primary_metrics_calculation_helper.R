#' Calculate a mean with empty-vector protection
#'
#' Returns the mean of a numeric vector while ignoring missing values. If the
#' input vector has length zero, the function returns \code{NA_real_}.
#'
#' @param x A numeric vector.
#'
#' @return A numeric value.
#'
#' @export

safe_mean <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

#' Calculate a correlation with safety checks
#'
#' Returns the Pearson correlation between two numeric vectors. If either vector
#' has fewer than two values or zero variance, the function returns \code{NA_real_}.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return A numeric value.
#'
#' @export
safe_cor <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  if (stats::sd(x, na.rm = TRUE) == 0 || stats::sd(y, na.rm = TRUE) == 0) {
    return(NA_real_)
  }
  stats::cor(x, y, use = "complete.obs")
}


#' Calculate primary diagnostic metrics from a community matrix
#'
#' Calculates process-proximal diagnostic metrics from a fruit-level community
#' matrix, including total abundance, species-specific abundance, occurrence
#' frequency, guild-level contribution, host-parasitoid conditional relationships,
#' mean-variance relationships, and extreme-value figs.
#'
#' @param comm_mat A matrix or data frame with figs as rows and wasp species as columns.
#' @param dataset_name Character. Label for the dataset, such as `"Observed"` or `"Simulated"`.
#' @param guild_lookup Optional data frame with columns `species` and `guild`.
#' @param host_pairs Optional data frame with columns `parasitoid` and `host`.
#' @param low_total_cutoff Optional numeric cutoff defining very low total abundance.
#' @param high_total_cutoff Optional numeric cutoff defining very high total abundance.
#' @param dominance_threshold Numeric. Relative-abundance threshold used to define single-species dominance.
#'
#' @return A list of data frames containing primary diagnostic summaries.
#' @export
calc_primary_metrics <- function(comm_mat,
                                 dataset_name,
                                 guild_lookup = NULL,
                                 host_pairs = NULL,
                                 low_total_cutoff = NULL,
                                 high_total_cutoff = NULL,
                                 dominance_threshold = 0.80) {

  comm_mat <- as.data.frame(comm_mat)
  comm_mat[] <- lapply(comm_mat, as.numeric)

  species_cols <- colnames(comm_mat)
  n_figs <- nrow(comm_mat)

  # ------------------------------------------------------------
  # 1) Total abundance per fig
  # ------------------------------------------------------------
  total_per_fig <- tibble(
    dataset = dataset_name,
    fig_id = seq_len(n_figs),
    total_abundance = rowSums(comm_mat, na.rm = TRUE)
  )

  total_summary <- total_per_fig %>%
    summarise(
      dataset = dataset_name,
      n_figs = n(),
      mean_total_abundance = mean(total_abundance, na.rm = TRUE),
      median_total_abundance = median(total_abundance, na.rm = TRUE),
      var_total_abundance = var(total_abundance, na.rm = TRUE),
      sd_total_abundance = sd(total_abundance, na.rm = TRUE),
      q05_total_abundance = quantile(total_abundance, 0.05, na.rm = TRUE),
      q25_total_abundance = quantile(total_abundance, 0.25, na.rm = TRUE),
      q75_total_abundance = quantile(total_abundance, 0.75, na.rm = TRUE),
      q95_total_abundance = quantile(total_abundance, 0.95, na.rm = TRUE)
    )

  # ------------------------------------------------------------
  # 2) Species-specific abundance per fig
  # ------------------------------------------------------------
  species_abundance_long <- comm_mat %>%
    mutate(fig_id = seq_len(n_figs)) %>%
    pivot_longer(
      cols = all_of(species_cols),
      names_to = "species",
      values_to = "abundance"
    ) %>%
    mutate(dataset = dataset_name) %>%
    select(dataset, fig_id, species, abundance)

  species_abundance_summary <- species_abundance_long %>%
    group_by(dataset, species) %>%
    summarise(
      mean_abundance = mean(abundance, na.rm = TRUE),
      median_abundance = median(abundance, na.rm = TRUE),
      var_abundance = var(abundance, na.rm = TRUE),
      sd_abundance = sd(abundance, na.rm = TRUE),
      q05_abundance = quantile(abundance, 0.05, na.rm = TRUE),
      q25_abundance = quantile(abundance, 0.25, na.rm = TRUE),
      q75_abundance = quantile(abundance, 0.75, na.rm = TRUE),
      q95_abundance = quantile(abundance, 0.95, na.rm = TRUE),
      zero_proportion = mean(abundance == 0, na.rm = TRUE),
      .groups = "drop"
    )

  # ------------------------------------------------------------
  # 3) Frequency of occurrence per species
  # ------------------------------------------------------------
  occurrence_frequency <- species_abundance_long %>%
    group_by(dataset, species) %>%
    summarise(
      occurrence_frequency = mean(abundance > 0, na.rm = TRUE),
      zero_proportion = mean(abundance == 0, na.rm = TRUE),
      n_present_figs = sum(abundance > 0, na.rm = TRUE),
      n_total_figs = n(),
      .groups = "drop"
    )

  # ------------------------------------------------------------
  # 4) Relative contribution of functional guilds
  # ------------------------------------------------------------
  guild_per_fig <- NULL
  guild_summary <- NULL

  if (!is.null(guild_lookup)) {
    guild_per_fig <- species_abundance_long %>%
      left_join(guild_lookup, by = "species") %>%
      group_by(dataset, fig_id, guild) %>%
      summarise(
        guild_abundance = sum(abundance, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      group_by(dataset, fig_id) %>%
      mutate(
        total_abundance = sum(guild_abundance, na.rm = TRUE),
        guild_relative_contribution = ifelse(
          total_abundance > 0,
          guild_abundance / total_abundance,
          NA_real_
        )
      ) %>%
      ungroup()

    guild_summary <- guild_per_fig %>%
      group_by(dataset, guild) %>%
      summarise(
        mean_guild_abundance = mean(guild_abundance, na.rm = TRUE),
        median_guild_abundance = median(guild_abundance, na.rm = TRUE),
        var_guild_abundance = var(guild_abundance, na.rm = TRUE),
        sd_guild_abundance = sd(guild_abundance, na.rm = TRUE),
        mean_relative_contribution = mean(guild_relative_contribution, na.rm = TRUE),
        median_relative_contribution = median(guild_relative_contribution, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # ------------------------------------------------------------
  # 5) Conditional relationships between species
  # ------------------------------------------------------------
  conditional_per_fig <- NULL
  conditional_summary <- NULL


  if (!is.null(host_pairs) && nrow(host_pairs) > 0) {

    # Make sure host_pairs columns are plain character vectors
    host_pairs <- host_pairs %>%
      mutate(
        parasitoid = as.character(parasitoid),
        host = as.character(host)
      )

    conditional_per_fig <- purrr::map_dfr(seq_len(nrow(host_pairs)), function(i) {

      parasitoid <- host_pairs$parasitoid[[i]]
      host <- host_pairs$host[[i]]

      if (!parasitoid %in% species_cols || !host %in% species_cols) {
        stop(
          "Host or parasitoid name not found in community matrix: ",
          parasitoid, " - ", host
        )
      }

      host_abundance <- comm_mat[, host, drop = TRUE]
      parasitoid_abundance <- comm_mat[, parasitoid, drop = TRUE]

      tibble::tibble(
        dataset = dataset_name,
        fig_id = seq_len(n_figs),
        parasitoid = parasitoid,
        host = host,
        host_abundance = as.numeric(host_abundance),
        parasitoid_abundance = as.numeric(parasitoid_abundance),
        host_present = host_abundance > 0,
        parasitoid_present = parasitoid_abundance > 0
      )
    })

    conditional_summary <- conditional_per_fig %>%
      group_by(dataset, parasitoid, host) %>%
      summarise(
        n_figs = n(),
        n_host_present = sum(host_present, na.rm = TRUE),
        n_parasitoid_present = sum(parasitoid_present, na.rm = TRUE),
        p_host_present = mean(host_present, na.rm = TRUE),
        p_parasitoid_present = mean(parasitoid_present, na.rm = TRUE),
        p_parasitoid_given_host_present =
          safe_mean(parasitoid_present[host_present]),
        p_parasitoid_given_host_absent =
          safe_mean(parasitoid_present[!host_present]),
        mean_parasitoid_abundance_when_host_present =
          safe_mean(parasitoid_abundance[host_present]),
        mean_parasitoid_abundance_when_host_absent =
          safe_mean(parasitoid_abundance[!host_present]),
        host_parasitoid_abundance_cor =
          safe_cor(host_abundance, parasitoid_abundance),
        .groups = "drop"
      )
  }

  # ------------------------------------------------------------
  # 6) Mean-variance relationship
  # ------------------------------------------------------------
  mean_variance_species <- species_abundance_summary %>%
    transmute(
      dataset,
      level = "species",
      group = species,
      mean_abundance = mean_abundance,
      variance_abundance = var_abundance
    )

  mean_variance_total <- total_per_fig %>%
    summarise(
      dataset = dataset_name,
      level = "total",
      group = "all_species",
      mean_abundance = mean(total_abundance, na.rm = TRUE),
      variance_abundance = var(total_abundance, na.rm = TRUE)
    )

  mean_variance_guild <- NULL

  if (!is.null(guild_per_fig)) {
    mean_variance_guild <- guild_per_fig %>%
      group_by(dataset, guild) %>%
      summarise(
        level = "guild",
        group = first(guild),
        mean_abundance = mean(guild_abundance, na.rm = TRUE),
        variance_abundance = var(guild_abundance, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      select(dataset, level, group, mean_abundance, variance_abundance)
  }

  mean_variance_relationship <- bind_rows(
    mean_variance_total,
    mean_variance_species,
    mean_variance_guild
  )

  # ------------------------------------------------------------
  # 7) Frequency of figs with extreme values
  # ------------------------------------------------------------
  if (is.null(low_total_cutoff)) {
    low_total_cutoff <- quantile(total_per_fig$total_abundance, 0.05, na.rm = TRUE)
  }

  if (is.null(high_total_cutoff)) {
    high_total_cutoff <- quantile(total_per_fig$total_abundance, 0.95, na.rm = TRUE)
  }

  max_species_abundance <- apply(comm_mat, 1, max, na.rm = TRUE)
  dominant_species <- species_cols[apply(comm_mat, 1, which.max)]
  dominant_share <- ifelse(
    total_per_fig$total_abundance > 0,
    max_species_abundance / total_per_fig$total_abundance,
    NA_real_
  )

  extreme_per_fig <- total_per_fig %>%
    mutate(
      low_total_cutoff = as.numeric(low_total_cutoff),
      high_total_cutoff = as.numeric(high_total_cutoff),
      no_wasps = total_abundance == 0,
      very_low_abundance = total_abundance <= low_total_cutoff,
      very_high_abundance = total_abundance >= high_total_cutoff,
      dominant_species = dominant_species,
      dominant_share = dominant_share,
      single_species_dominated = dominant_share >= dominance_threshold
    )

  if (!is.null(guild_lookup)) {
    pollinator_species <- guild_lookup %>%
      filter(guild == "pollinator") %>%
      pull(species)

    if (length(pollinator_species) > 0) {
      pollinator_abundance <- rowSums(
        comm_mat[, intersect(pollinator_species, species_cols), drop = FALSE],
        na.rm = TRUE
      )

      extreme_per_fig <- extreme_per_fig %>%
        mutate(
          pollinator_abundance = pollinator_abundance,
          no_pollinator = pollinator_abundance == 0
        )
    }
  }

  extreme_summary <- extreme_per_fig %>%
    summarise(
      dataset = dataset_name,
      n_figs = n(),
      low_total_cutoff = first(low_total_cutoff),
      high_total_cutoff = first(high_total_cutoff),
      prop_no_wasps = mean(no_wasps, na.rm = TRUE),
      prop_very_low_abundance = mean(very_low_abundance, na.rm = TRUE),
      prop_very_high_abundance = mean(very_high_abundance, na.rm = TRUE),
      prop_single_species_dominated = mean(single_species_dominated, na.rm = TRUE),
      mean_dominant_share = mean(dominant_share, na.rm = TRUE),
      median_dominant_share = median(dominant_share, na.rm = TRUE),
      prop_no_pollinator = if ("no_pollinator" %in% names(extreme_per_fig)) {
        mean(extreme_per_fig$no_pollinator, na.rm = TRUE)
      } else {
        NA_real_
      }
    )

  return(list(
    total_per_fig = total_per_fig,
    total_summary = total_summary,
    species_abundance_long = species_abundance_long,
    species_abundance_summary = species_abundance_summary,
    occurrence_frequency = occurrence_frequency,
    guild_per_fig = guild_per_fig,
    guild_summary = guild_summary,
    conditional_per_fig = conditional_per_fig,
    conditional_summary = conditional_summary,
    mean_variance_relationship = mean_variance_relationship,
    extreme_per_fig = extreme_per_fig,
    extreme_summary = extreme_summary
  ))
}

#' Calculate species richness per fig
#'
#' Calculates the number of species present in each fig from a community matrix.
#'
#' @param comm_mat A matrix or data frame with figs as rows and wasp species as columns.
#' @param dataset_name Character. Label for the dataset, such as `"Observed"` or `"Simulated"`.
#'
#' @return A data frame with columns `dataset`, `fig_id`, and `richness`.
#' @export
calc_richness_per_fig <- function(comm_mat, dataset_name) {

  comm_mat <- as.data.frame(comm_mat)
  comm_mat[is.na(comm_mat)] <- 0

  data.frame(
    dataset = dataset_name,
    fig_id = seq_len(nrow(comm_mat)),
    richness = rowSums(comm_mat > 0)
  )
}
