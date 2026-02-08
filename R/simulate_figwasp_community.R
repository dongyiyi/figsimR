#' Simulate Fig Wasp Community Assembly (Observed on Ficus racemosa)
#'
#' This function simulates fig wasp community structure within a set of fig fruits,
#' incorporating priority effects, host–parasitoid interactions, resource constraints,
#' pollination, oviposition success, and fig abortion mechanisms.
#'
#' Each simulation models community assembly in a number of fig fruits, where:
#' \itemize{
#'   \item Flower number is determined by fig diameter and stochastic heterogeneity.
#'   \item Wasp species arrive in \code{entry_priority} phases (priority effects).
#'   \item Pollinators and gallers lay eggs into limited ovules.
#'   \item Parasitoids require host presence and attack host eggs.
#'   \item Layer preference can restrict oviposition to \code{core}, \code{mid}, or \code{outer} layers.
#'   \item Figs may abort (drop) if flower usage exceeds \code{host_sanction}.
#' }
#'
#' Species roles, hosts, and parasitoid links must be defined in \code{species_roles}.
#' Default species and parameter values are provided in \code{\link{species_list}} and \code{\link{parameter_list_default}}.
#'
#' @param num_figs Integer. Number of fig fruits to simulate. Each fig will be treated as a discrete community.
#' @param fig_diameter_mean Numeric. Mean fig diameter (in cm).
#' @param fig_diameter_sd Numeric. SD of fig diameter.
#' @param k Numeric. Scaling constant for estimating flower number from fig diameter.
#' @param alpha Numeric. Exponent in the diameter–flower power-law.
#' @param max_entry_table Named numeric vector. Per-species entry caps (scaled by fig size).
#' @param fecundity_mean, fecundity_dispersion Named numeric vectors for NB fecundity.
#' @param entry_mu, entry_size Named numeric vectors controlling entry distribution.
#' @param entry_priority Named list of phases and species.
#' @param species_roles List with $guild, $hosts, $parasitoid.
#' @param entry_distribution "nb" or "lognormal".
#' @param interaction_matrix Optional numeric matrix (pairwise interactions).
#' @param interaction_weight Numeric. Strength of interaction effects.
#' @param seed RNG seed (optional).
#' @param egg_success_prob Global oviposition success (if not using phase-specific).
#' @param egg_success_prob_by_phase List of per-phase success probs by species.
#' @param parasitism_prob Named numeric vector for supplemental parasitism.
#' @param layer_preference Named list of layer probabilities per species.
#' @param p_pollination_per_ovule Numeric [0–1]. Background pollination to seeds.
#' @param p_no_entry Numeric [0–1]. Probability a fig receives no entries.
#' @param host_sanction Numeric [0–1]. Legacy overuse threshold for drop logic.
#' @param fig_diameter_min, fig_diameter_max Numeric. Truncation bounds (if used).
#' @param record_individual Logical. If TRUE, per-individual oviposition recorded.
#' @param use_layering,use_layer_preference Logical switches for layering.
#' @param enable_drop Logical. Enable legacy drop (host_sanction) mechanism.
#' @param drop_cancels_emergence Logical. If FALSE, dropped figs yield zero output.
#' @param use_supplemental_parasitism Logical.
#' @param use_flower_limit Logical.
#' @param use_egg_success_by_phase Logical.
#'
#' @param use_sink_strength Logical. If TRUE, compute sink-strength metrics but
#'   do NOT alter legacy drop decision. Adds columns: sink_strength, sink_prop,
#'   sink_below_min, sink_above_max.
#' @param sink_w_gall Numeric in [0, +). Per-ovule sink weight for galled ovules (default 1.0).
#' @param sink_w_seed Numeric in [0, +). Per-ovule sink weight for seeds (default 1.5).
#' @param sink_linear_coef Numeric. Global multiplier on sink strength (default 1.0).
#' @param sink_min_prop Numeric in [0,1]. Reference lower bound of sink proportion (default 0.20).
#' @param sink_max_prop Numeric in [0,1]. Reference upper bound of sink proportion (default 0.95).
#'
#' @return list(summary = data.frame(...), individual_eggs = list() if requested)
#' @examples
#' data(species_list)
#' data(parameter_list_default)
#' set.seed(123)
#' sim_result <- simulate_figwasp_community(
#'   num_figs = 10,
#'   fecundity_mean = parameter_list_default$fecundity_mean,
#'   fecundity_dispersion = parameter_list_default$fecundity_dispersion,
#'   entry_mu = parameter_list_default$entry_mu,
#'   entry_size = parameter_list_default$entry_size,
#'   entry_priority = parameter_list_default$entry_priority,
#'   species_roles = parameter_list_default$species_roles,
#'   max_entry_table = parameter_list_default$max_entry_table,
#'   interaction_matrix = parameter_list_default$interaction_matrix,
#'   interaction_weight = parameter_list_default$interaction_weight,
#'   egg_success_prob = parameter_list_default$egg_success_prob,
#'   egg_success_prob_by_phase = parameter_list_default$egg_success_prob_by_phase,
#'   parasitism_prob = parameter_list_default$parasitism_prob,
#'   layer_preference = parameter_list_default$layer_preference,
#'   record_individual = TRUE,
#'   use_sink_strength = TRUE,
#'   sink_w_gall = 1.0,
#'   sink_w_seed = 1.5,
#'   sink_linear_coef = 1.0,
#'   sink_min_prop = 0.20,
#'   sink_max_prop = 0.95
#' )
#' head(sim_result$summary)
#' @export
simulate_figwasp_community <- function(
    num_figs = 1000,
    fig_diameter_mean = 2.5,
    fig_diameter_sd = 1.2,
    k = 300,
    alpha = 1.3,
    max_entry_table = NULL,
    fecundity_mean,
    fecundity_dispersion,
    entry_mu,
    entry_size,
    entry_priority,
    species_roles,
    entry_distribution = "lognormal",
    interaction_matrix = NULL,
    interaction_weight = 0,
    seed = NULL,
    egg_success_prob = NULL,
    egg_success_prob_by_phase = NULL,
    parasitism_prob = NULL,
    layer_preference = NULL,
    p_pollination_per_ovule = 0.98,
    p_no_entry = 0.002,
    host_sanction = 0.8,
    fig_diameter_min = 1.3,
    fig_diameter_max = 4.0,
    record_individual = FALSE,

    # module control
    use_layering = TRUE,
    use_layer_preference = TRUE,
    enable_drop = TRUE,
    drop_cancels_emergence = FALSE,
    use_supplemental_parasitism = FALSE,
    use_flower_limit = TRUE,
    use_egg_success_by_phase = TRUE,

    # sink metrics (diagnostic-only)
    use_sink_strength = FALSE,
    sink_w_gall = 1.0,
    sink_w_seed = 1.5,
    sink_linear_coef = 1.0,
    sink_min_prop = 0.20,
    sink_max_prop = 0.95
) {
  if (!is.null(seed)) set.seed(seed)

  species_names <- names(species_roles$guild)

  if (is.null(max_entry_table)) {
    max_entry_table <- c(
      Sycophaga_testacea = 20,
      Apocrypta_sp = 20,
      Sycophaga_mayri = 20,
      Ceratosolen_sp = 20,
      Sycophaga_agraensis = 20,
      Apocrypta_westwoodi = 20
    )
  }

  df <- data.frame(fig_id = 1:num_figs)
  diameters <- rnorm(num_figs, mean = fig_diameter_mean, sd = fig_diameter_sd)
  if (use_flower_limit) {
    diameters <- pmin(pmax(diameters, fig_diameter_min), fig_diameter_max)
  }
  df$fig_diameter <- diameters
  heterogeneity_noise <- rgamma(num_figs, shape = 20, scale = 0.1)
  df$flower_count <- round(k * df$fig_diameter^alpha * heterogeneity_noise)
  df$resource_use <- 0L
  no_entry_vec <- runif(num_figs) < p_no_entry
  df$richness_skipped <- as.integer(no_entry_vec)

  for (sp in species_names) {
    df[[paste0("entry_", sp)]] <- 0L
    df[[paste0("eggs_", sp)]] <- 0L
  }

  entry_mat <- matrix(0L, nrow = num_figs, ncol = length(species_names))
  colnames(entry_mat) <- species_names
  if (record_individual) individual_records <- list()

  if (use_layering && use_layer_preference && is.null(layer_preference)) {
    layer_preference <- list(
      Ceratosolen_sp = c(core = 0.6, mid = 0.3, outer = 0.1),
      Sycophaga_mayri = c(core = 0.1, mid = 0.4, outer = 0.5),
      Sycophaga_testacea = c(core = 0.1, mid = 0.4, outer = 0.5),
      Apocrypta_sp = c(core = 0.1, mid = 0.4, outer = 0.5),
      Apocrypta_westwoodi = c(core = 0.1, mid = 0.4, outer = 0.5),
      Sycophaga_agraensis = c(core = 0.4, mid = 0.4, outer = 0.2)
    )
  }

  for (phase_name in names(entry_priority)) {
    phase <- entry_priority[[phase_name]]
    for (sp in phase) {
      guild <- species_roles$guild[sp]
      hosts <- species_roles$hosts[[sp]]
      fec_mean <- fecundity_mean[sp]
      fec_disp <- fecundity_dispersion[sp]
      prob_active <- 1
      if (use_egg_success_by_phase) {
        if (!is.null(egg_success_prob_by_phase[[sp]][[phase_name]])) {
          prob_active <- egg_success_prob_by_phase[[sp]][[phase_name]]
        }
      } else {
        if (!is.null(egg_success_prob[sp])) {
          prob_active <- egg_success_prob[sp]
        }
      }

      host_present <- if (length(hosts) > 0) rowSums(entry_mat[, hosts, drop = FALSE]) > 0 else rep(TRUE, num_figs)

      interaction_factor <- rep(1, num_figs)
      if (!is.null(interaction_matrix)) {
        for (j in 1:num_figs) {
          present_vec <- entry_mat[j, , drop = FALSE]
          adjustment <- sum(interaction_matrix[sp, ] * (present_vec > 0))
          interaction_factor[j] <- exp(interaction_weight * adjustment)
        }
      }

      n_vec <- rep(0L, num_figs)
      idx <- which(host_present & !no_entry_vec)

      if (entry_distribution == "nb") {
        mu_vec <- entry_mu[sp] * interaction_factor[idx]
        n_vec[idx] <- rnbinom(length(idx), size = entry_size[sp], mu = mu_vec)
      } else {
        meanlog <- log(entry_mu[sp]) - 0.5
        sdlog <- 1 / sqrt(entry_size[sp])
        raw <- rlnorm(length(idx), meanlog = meanlog, sdlog = sdlog) * interaction_factor[idx]
        n_vec[idx] <- rpois(length(idx), lambda = raw)
      }

      max_entry_vec <- round(max_entry_table[sp] * df$fig_diameter / fig_diameter_mean)
      n_vec <- pmin(n_vec, max_entry_vec)

      entry_mat[, sp] <- entry_mat[, sp] + n_vec
      df[[paste0("entry_", sp)]] <- df[[paste0("entry_", sp)]] + n_vec

      eggs_vec <- rep(0L, num_figs)
      for (j in which(n_vec > 0)) {
        n <- n_vec[j]
        active <- rbinom(n, 1, prob_active)
        individual_eggs <- rnbinom(n, size = fec_disp, mu = fec_mean) * active
        if (record_individual) {
          individual_records[[paste0("fig", j, "_", sp)]] <- individual_eggs
        }
        if (guild %in% c("pollinator", "galler")) {
          available <- df$flower_count[j] - df$resource_use[j]
          if (use_layering && !is.null(layer_preference[[sp]])) {
            layer_prob <- layer_preference[[sp]]
            assigned_layer <- sample(names(layer_prob), size = length(individual_eggs), replace = TRUE, prob = layer_prob)
            individual_eggs <- ifelse(assigned_layer == "core", individual_eggs, 0)
          }
          eggs_cumsum <- cumsum(individual_eggs)
          if (any(eggs_cumsum > available)) {
            cutoff <- which(eggs_cumsum > available)[1] - 1
            individual_eggs <- if (cutoff > 0) individual_eggs[1:cutoff] else integer(0)
          }
          total_eggs <- sum(individual_eggs)
          df$resource_use[j] <- df$resource_use[j] + total_eggs
          eggs_vec[j] <- total_eggs
        } else if (guild == "parasitoid" && length(hosts) > 0) {
          host_eggs <- sum(sapply(hosts, function(h) df[[paste0("eggs_", h)]][j]))
          total_eggs <- min(sum(individual_eggs), host_eggs)
          eggs_vec[j] <- total_eggs
        }
      }
      df[[paste0("eggs_", sp)]] <- df[[paste0("eggs_", sp)]] + eggs_vec
    }
  }

  if (use_supplemental_parasitism) {
    for (para in species_names[species_roles$guild == "parasitoid"]) {
      hosts <- species_roles$hosts[[para]]
      if (length(hosts) == 0) next
      p_para <- parasitism_prob[[para]] %||% 0.9
      for (j in 1:num_figs) {
        if (sum(sapply(hosts, function(h) df[[paste0("eggs_", h)]][j])) > 0 && runif(1) < p_para) {
          n <- rpois(1, lambda = entry_mu[para])
          fec <- rnbinom(n, size = fecundity_dispersion[para], mu = fecundity_mean[para])
          host_eggs <- sum(sapply(hosts, function(h) df[[paste0("eggs_", h)]][j]))
          total_eggs <- min(sum(fec), host_eggs)
          df[[paste0("eggs_", para)]][j] <- df[[paste0("eggs_", para)]][j] + total_eggs
          df[[paste0("entry_", para)]][j] <- df[[paste0("entry_", para)]][j] + 1
        }
      }
    }
  }

  # Emergence after parasitism
  for (sp in species_names) {
    parasitoids <- species_roles$parasitoid[[sp]]
    parasite_total <- if (length(parasitoids) == 0) rep(0, nrow(df)) else {
      egg_cols <- paste0("eggs_", parasitoids)
      egg_cols <- egg_cols[egg_cols %in% names(df)]
      if (length(egg_cols) == 0) rep(0, nrow(df)) else rowSums(df[, egg_cols, drop = FALSE], na.rm = TRUE)
    }
    df[[paste0("emergence_", sp)]] <- pmax(df[[paste0("eggs_", sp)]] - parasite_total, 0)
  }

  # Seeds & failures from unoccupied flowers
  df$unoccupied_flowers <- pmax(df$flower_count - df$resource_use, 0)
  df$seeds <- rbinom(n = num_figs, size = df$unoccupied_flowers, prob = p_pollination_per_ovule)
  df$failed_ovules <- df$unoccupied_flowers - df$seeds
  df$used_flowers <- df$resource_use + df$seeds + df$failed_ovules
  df$resource_ratio <- ifelse(df$flower_count > 0, df$resource_use / df$flower_count, 0)

  # --- sink-strength metrics ---
  if (use_sink_strength) {
    pollinators <- names(species_roles$guild)[species_roles$guild == "pollinator"]
    gallers     <- names(species_roles$guild)[species_roles$guild == "galler"]

    galled_total <- 0L
    if (length(pollinators)) galled_total <- galled_total + rowSums(df[, paste0("eggs_", pollinators), drop = FALSE])
    if (length(gallers))     galled_total <- galled_total + rowSums(df[, paste0("eggs_", gallers),     drop = FALSE])

    seed_count <- df$seeds

    sink_strength <- sink_linear_coef * (sink_w_gall * galled_total + sink_w_seed * seed_count)
    w_max <- max(sink_w_seed, sink_w_gall)
    sink_prop <- ifelse(df$flower_count > 0, sink_strength / (w_max * df$flower_count), 0)

    df$sink_strength  <- as.numeric(sink_strength)
    df$sink_prop      <- as.numeric(sink_prop)
    df$sink_below_min <- sink_prop < sink_min_prop
    df$sink_above_max <- sink_prop > sink_max_prop
  } else {
    df$sink_strength  <- NA_real_
    df$sink_prop      <- NA_real_
    df$sink_below_min <- NA
    df$sink_above_max <- NA
  }


  if (enable_drop) {
    drop_slope <- 10
    df$drop_prob <- plogis((df$resource_ratio - host_sanction) * drop_slope)
    df$is_dropped <- rbinom(nrow(df), 1, df$drop_prob)
  } else {
    df$drop_prob <- rep(0, nrow(df))
    df$is_dropped <- rep(0L, nrow(df))
  }


  if (enable_drop && !drop_cancels_emergence) {
    dropped_idx <- which(df$is_dropped == 1)
    if (length(dropped_idx) > 0) {
      emergence_cols <- grep("^emergence_", names(df), value = TRUE)
      df[dropped_idx, emergence_cols] <- 0
      df$seeds[dropped_idx] <- 0
      df$failed_ovules[dropped_idx] <- 0
    }
  }

  df$entry_matrix <- split(entry_mat, row(entry_mat))

  if (record_individual) {
    return(list(summary = df, individual_eggs = individual_records))
  } else {
    return(list(summary = df))
  }
}
