#' Simulate Fig Wasp Community Assembly (Observed on \emph{Ficus racemosa})
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
#'   \item The ovary-layer preference filter can be bypassed with \code{"none"}; Or layer preference can restrict oviposition to \code{core}, \code{mid}, or \code{outer} layers.
#'   \item Figs may abort (drop) if flower usage exceeds \code{host_sanction} or \code{sink_strength}.
#' }
#'
#' Species roles, hosts, and parasitoid links must be defined in \code{species_roles}.
#' Default species and parameter values are provided in \code{\link{species_list}} and \code{\link{parameter_list_default}}.
#'
#' @param num_figs Integer. Number of fig fruits to simulate (i.e., number of independent communities, default = 1000). Each fig will be treated as a discrete community.
#' @param fig_diameter_mean Numeric. Mean fig diameter (in cm, default = 2.5). It is used to determine ovule number. Affects resource availability.
#' @param fig_diameter_sd Numeric. Standard deviation of fig diameter under the normal distribution used to generate fig diameters, default = 1.2).
#' @param k Numeric. Scaling constant for estimating flower number from fig diameter (default = 300). Controls ovule-based competition.
#' @param alpha Numeric. Exponent in the diameter–flower power-law used to estimate flower number from fig diameter; controls shape of seed-yield curve (default = 1.3).
#' @param max_entry_table Named numeric vector giving the baseline maximum number of
#'   arrivals or oviposition-attempting wasps per species. Names must match
#'   \code{species_roles$guild}. If \code{NULL}, the function uses the default
#'   \emph{Ficus racemosa} values when the supplied species set matches the worked
#'   example; otherwise, a generic value is assigned to all species, with a warning recommending system-specific values.
#'
#' @param fecundity_mean The mean potential number of eggs per individual wasp.
#' @param fecundity_dispersion Named numeric vectors for fecundity distribution. The dispersion parameter for the negative binomial distribution modeling individual fecundity.
#' @param entry_mu Named numeric vectors controlling the expected number of wasp arrival
#'   or oviposition-attempt events per fig. For pollinators, this corresponds to physical
#'   entry into the fig cavity; for non-pollinating fig wasps, it represents arrival at the
#'   fig and oviposition attempts from outside the fig wall. \code{entry_mu} should not be interpreted as physical
#'   entry into the fig cavity for all species.
#' @param entry_size Named numeric vectors controlling entry distribution (the dispersion of the
#'   arrival or oviposition-attempt distribution). See \code{entry_distribution}.
#'   For \code{entry_distribution = "nb"}: used as NB size parameter. For \code{lognormal}: used to derive \code{sdlog = 1 / sqrt(entry_size)}.
#'   In the negative binomial option, this corresponds to the "size" parameter; smaller values produce stronger aggregation and greater among-fig variation.
#'   In the lognormal option, it controls the log-scale variance of the arrival process.
#'   This parameter is not a sample size and does not represent the number of eggs.
#' @param entry_priority Named list of phases and species.
#' @param species_roles List with $guild, $hosts, $parasitoid.
#' @param entry_distribution "nb" or "lognormal". See \code{entry_size}.
#' @param interaction_matrix Optional numeric matrix (pairwise interactions).
#' @param interaction_weight Numeric. Strength of interaction effects.
#' @param seed RNG seed (optional).
#' @param egg_success_prob Global probability that a single oviposition attempt results in a successful egg. Applies when \code{use_egg_success_by_phase = FALSE}.
#' @param use_egg_success_by_phase List of per-phase success probs by species. Default \code{TRUE}.
#' @param egg_success_prob_by_phase Named list of species-specific, phase-specific
#'   egg-success probabilities. Each element is a named numeric vector whose names
#'   correspond to phases in \code{entry_priority}. When
#'   \code{use_egg_success_by_phase = TRUE}, a matching species-phase value
#'   overrides the baseline value in \code{egg_success_prob}. If no phase-specific
#'   value is provided, the baseline \code{egg_success_prob} is used; if that is
#'   also missing, the default success probability is 1.
#'
#' @param parasitism_prob Named numeric vector for supplemental parasitism.
#' @param layer_preference Named list specifying species-specific ovary-layer preference or accessibility. Each species can be assigned either a numeric
#'   vector giving relative probabilities for \code{core}, \code{mid}, or \code{outer} ovary
#'   layers; Or \code{none} to indicate no species-specific ovary-layer preference.
#'   If \code{NULL}, figsimR uses the default \emph{Ficus racemosa}-specific
#'   layer preference when the supplied species set matches the worked example;
#'   otherwise, species are assigned \code{none} and users are encouraged to provide system-specific values.
#' @param p_pollination_per_ovule Numeric [0–1]. Probability that each unoccupied ovule becomes a seed.
#' @param p_no_entry Numeric [0–1]. Probability a fig receives no entries.
#' @param host_sanction Numeric [0–1]. Legacy overuse threshold for drop logic.
#' @param fig_diameter_min Numeric. Minimum fig diameter used to cap values (optional). Truncation bounds (if used). Truncate simulated diameters (and thus flower counts) to a biologically reasonable range.
#' @param fig_diameter_max Numeric. Maximum fig diameter used to cap values (optional). Truncation bounds (if used). Truncate simulated diameters (and thus flower counts) to a biologically reasonable range.
#' @param record_individual Logical. If \code{TRUE}, per-individual oviposition recorded.
#' @param use_layering Logical switches for layering. If \code{TRUE}, fig flowers are divided into spatial layers (core, mid, outer).
#' @param use_layer_preference Logical. If \code{TRUE}, apply species-specific
#'   ovary-layer preference or accessibility values supplied through
#'   \code{layer_preference}. If \code{FALSE}, ovary-layer resources are still used when
#'   \code{use_layering = TRUE}, but no species-specific layer preference is applied.
#' @param enable_drop Logical. Enable legacy drop (host sanction) mechanism. If flower use proportion exceeds host_sanction, fig may abort.
#'   In the output "drop_prob" column: Ranges from 0 (never drop) to 1 (always drop); "is_dropped" column: no drop (0) and drop (1).
#'   When sink strength is used (\code{use_sink_strength = TRUE}) and \code{enable_drop = FALSE} (host sanction is disabled); the columns "is_dropped" and "drop_prob" in the outputs are not meaningful.
#' @param drop_cancels_emergence Logical. If \code{FALSE} (figs are dropped and no wasps emerge), dropped figs yield zero output (no wasps emerge at the end of phase).
#' @param use_supplemental_parasitism Logical. Whether to activate supplemental parasitism. Default FALSE.
#' @param use_flower_limit Logical. If TRUE, clamps fig diameter between \code{fig_diameter_min} and \code{fig_diameter_max} before calculating flower number.
#' @param use_egg_success_by_phase Logical. Use egg success probabilities by phase. Default TRUE.
#'
#' @param use_sink_strength Logical. If \code{TRUE}, compute sink-strength metrics but
#'   do NOT alter legacy drop decision: \code{use_sink_strength = TRUE} and \code{enable_drop = FALSE} and the columns "is_dropped" and "drop_prob" are not meaningful.
#'   The fig is dropped when the value of "sink_prop" falls outside the thresholds defined by \code{sink_min_prop} and \code{sink_max_prop}; The columns "sink_below_min" and "sink_above_max" are logical: TRUE = drop; FALSE = no drop.
#'   Adds columns: \code{sink_strength}, \code{sink_prop},
#'   \code{sink_below_min}, \code{sink_above_max}. Sink strength threshold calculation: \code{sink_strength = sink_linear_coef * (sink_w_gall * galled_total + sink_w_seed * seed_count)}.
#' @param sink_w_gall Numeric in [0, +). Per-ovule sink weight for galled ovules (default 1.0).
#' @param sink_w_seed Numeric in [0, +). Per-ovule sink weight for seeds (default 1.5).
#' @param sink_linear_coef Numeric. Global multiplier on sink strength (default 1.0).
#' @param sink_min_prop Numeric in [0,1]. Reference lower bound of sink proportion (default 0.20).
#' @param sink_max_prop Numeric in [0,1]. Reference upper bound of sink proportion (default 0.95).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{summary}: a compact data frame containing fig-level variables,
#'   total entry, egg, and emergence counts for each species, seed and ovule
#'   outcomes, sink-strength metrics, drop status, and the entry matrix.
#'   \item \code{original_summary}: the complete simulation output, including
#'   layer-specific resource-use and egg-count columns and other diagnostic
#'   variables.
#'   \item \code{individual_eggs}: a list of per-individual egg counts, returned
#'   only when \code{record_individual = TRUE}.
#' }
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

    # sink metrics (diagnostic-only; do not alter drop logic)
    use_sink_strength = FALSE,
    sink_w_gall = 1.0,
    sink_w_seed = 1.5,
    sink_linear_coef = 1.0,
    sink_min_prop = 0.20,
    sink_max_prop = 0.95
) {
  if (!is.null(seed)) set.seed(seed)

  species_names <- names(species_roles$guild)

  check_named_vector <- function(x, object_name, species_names, positive = TRUE) {
    if (is.null(x)) {
      stop(object_name, " cannot be NULL.")
    }
    if (is.null(names(x))) {
      stop(object_name, " must be a named vector.")
    }
    missing_species <- setdiff(species_names, names(x))
    if (length(missing_species) > 0) {
      stop(
        object_name, " is missing values for: ",
        paste(missing_species, collapse = ", ")
      )
    }
    x <- x[species_names]
    if (any(is.na(x))) {
      stop(object_name, " contains NA values after matching species names.")
    }
    if (positive && any(x <= 0)) {
      stop(object_name, " must contain positive values.")
    }
    x
  }

  entry_mu <- check_named_vector(entry_mu, "entry_mu", species_names)
  entry_size <- check_named_vector(entry_size, "entry_size", species_names)
  fecundity_mean <- check_named_vector(fecundity_mean, "fecundity_mean", species_names)
  fecundity_dispersion <- check_named_vector(fecundity_dispersion, "fecundity_dispersion", species_names)
  entry_species <- unique(unlist(entry_priority))
  unknown_entry_species <- setdiff(entry_species, species_names)
  if (length(unknown_entry_species) > 0) {
    stop(
      "entry_priority contains species not found in species_roles$guild: ",
      paste(unknown_entry_species, collapse = ", ")
    )
  }

  if (!is.null(interaction_matrix)) {
    if (!is.matrix(interaction_matrix)) {
      stop("interaction_matrix must be a matrix.")
    }
    if (!all(species_names %in% rownames(interaction_matrix)) ||
        !all(species_names %in% colnames(interaction_matrix))) {
      stop("interaction_matrix rownames and colnames must include all species.")
    }
    interaction_matrix <- interaction_matrix[species_names, species_names, drop = FALSE]
  }


  is_no_layer_preference <- function(x) {
    is.character(x) && length(x) == 1 && tolower(x) == "none"
  }

  normalize_layer_preference <- function(x, sp) {

    if (is_no_layer_preference(x)) {
      return("none")
    }

    if (!is.numeric(x)) {
      stop("layer_preference for ", sp, " must be either 'none' or a numeric vector.")
    }

    required_layers <- c("core", "mid", "outer")

    if (is.null(names(x)) || !all(required_layers %in% names(x))) {
      stop(
        "layer_preference for ", sp,
        " must be a named numeric vector containing core, mid, and outer."
      )
    }

    x <- x[required_layers]

    if (any(is.na(x)) || any(x < 0)) {
      stop("layer_preference for ", sp, " cannot contain NA or negative values.")
    }

    if (sum(x) <= 0) {
      stop("layer_preference for ", sp, " must sum to a positive value.")
    }

    x / sum(x)
  }

  normalize_prob_vector <- function(x, required_names, object_name) {
    if (is.null(names(x)) || !all(required_names %in% names(x))) {
      stop(
        object_name,
        " must be a named numeric vector containing: ",
        paste(required_names, collapse = ", ")
      )
    }

    x <- x[required_names]

    if (any(is.na(x)) || any(x < 0)) {
      stop(object_name, " cannot contain NA or negative values.")
    }

    if (sum(x) <= 0) {
      stop(object_name, " must sum to a positive value.")
    }

    x / sum(x)
  }

  cap_eggs_by_available <- function(eggs, available) {
    if (length(eggs) == 0 || available <= 0) {
      return(integer(0))
    }

    eggs_cumsum <- cumsum(eggs)

    if (any(eggs_cumsum > available)) {
      cutoff <- which(eggs_cumsum > available)[1] - 1
      eggs <- if (cutoff > 0) eggs[1:cutoff] else integer(0)
    }

    eggs
  }

  if (is.null(max_entry_table)) {

    racemosa_species <- c(
      "Sycophaga_testacea",
      "Apocrypta_sp",
      "Sycophaga_mayri",
      "Ceratosolen_sp",
      "Sycophaga_agraensis",
      "Apocrypta_westwoodi"
    )

    is_racemosa_default <- setequal(species_names, racemosa_species)

    if (is_racemosa_default) {

      message(
        "No max_entry_table supplied. Using the default Ficus racemosa-specific max_entry_table."
      )

      max_entry_table <- c(
        Sycophaga_testacea = 20,
        Apocrypta_sp = 20,
        Sycophaga_mayri = 20,
        Ceratosolen_sp = 20,
        Sycophaga_agraensis = 20,
        Apocrypta_westwoodi = 20
      )

      max_entry_table <- max_entry_table[species_names]

    } else {

      warning(
        "No max_entry_table supplied for a non-Ficus racemosa community. ",
        "Using a generic default value of 20 for all species. ",
        "For biologically realistic simulations, provide a named max_entry_table."
      )

      max_entry_table <- setNames(rep(20, length(species_names)), species_names)
    }
  }
  max_entry_table <- check_named_vector(max_entry_table, "max_entry_table", species_names)

  df <- data.frame(fig_id = 1:num_figs)
  diameters <- rnorm(num_figs, mean = fig_diameter_mean, sd = fig_diameter_sd)
  if (use_flower_limit) {
    diameters <- pmin(pmax(diameters, fig_diameter_min), fig_diameter_max)
  }
  df$fig_diameter <- diameters
  heterogeneity_noise <- rgamma(num_figs, shape = 20, scale = 0.1)
  df$flower_count <- round(k * df$fig_diameter^alpha * heterogeneity_noise)
  df$resource_use <- 0L
  layer_names <- c("core", "mid", "outer")

  if (use_layering) {
    # Default ovary-layer capacity proportions.
    ovary_layer_prop <- c(core = 1/3, mid = 1/3, outer = 1/3)
    ovary_layer_prop <- normalize_prob_vector(
      ovary_layer_prop,
      required_names = layer_names,
      object_name = "ovary_layer_prop"
    )

    df$flower_core <- floor(df$flower_count * ovary_layer_prop["core"])
    df$flower_mid <- floor(df$flower_count * ovary_layer_prop["mid"])
    df$flower_outer <- df$flower_count - df$flower_core - df$flower_mid

    df$resource_use_core <- 0L
    df$resource_use_mid <- 0L
    df$resource_use_outer <- 0L
  }
  no_entry_vec <- runif(num_figs) < p_no_entry
  df$richness_skipped <- as.integer(no_entry_vec)

  for (sp in species_names) {
    df[[paste0("entry_", sp)]] <- 0L
    df[[paste0("eggs_", sp)]] <- 0L

    if (use_layering) {
      for (lay in layer_names) {
        df[[paste0("eggs_", sp, "_", lay)]] <- 0L
        df[[paste0("parasitized_", sp, "_", lay)]] <- 0L
      }
    }
  }

  if (use_layering) {
    for (para in species_names) {

      hosts_para <- species_roles$hosts[[para]]
      guild_para <- species_roles$guild[[para]]

      if (length(hosts_para) > 0 && guild_para == "parasitoid") {
        for (h in hosts_para) {
          for (lay in layer_names) {
            df[[paste0("eggs_", para, "_on_", h, "_", lay)]] <- 0L
          }
        }
      }
    }
  }

  entry_mat <- matrix(0L, nrow = num_figs, ncol = length(species_names))
  colnames(entry_mat) <- species_names
  if (record_individual) individual_records <- list()

  if (use_layering && use_layer_preference && is.null(layer_preference)) {

    racemosa_species <- c(
      "Sycophaga_testacea",
      "Apocrypta_sp",
      "Sycophaga_mayri",
      "Ceratosolen_sp",
      "Sycophaga_agraensis",
      "Apocrypta_westwoodi"
    )

    is_racemosa_default <- setequal(species_names, racemosa_species)

    if (is_racemosa_default) {

      message(
        "No layer_preference supplied. Using the default Ficus racemosa-specific ovary-layer preference."
      )

      layer_preference <- list(
        Ceratosolen_sp = c(core = 0.6, mid = 0.3, outer = 0.1),
        Sycophaga_mayri = c(core = 0.1, mid = 0.4, outer = 0.5),
        Sycophaga_testacea = c(core = 0.1, mid = 0.4, outer = 0.5),
        Apocrypta_sp = c(core = 0.1, mid = 0.4, outer = 0.5),
        Apocrypta_westwoodi = c(core = 0.1, mid = 0.4, outer = 0.5),
        Sycophaga_agraensis = c(core = 0.4, mid = 0.4, outer = 0.2)
      )

      layer_preference <- layer_preference[species_names]

    } else {

      warning(
        "No layer_preference supplied for a non-Ficus racemosa community. ",
        "Using 'none' for all species. ",
        "For biologically realistic simulations, provide a named layer_preference list."
      )

      layer_preference <- setNames(
        rep(list("none"), length(species_names)),
        species_names
      )
    }
  }

  if (use_layering && use_layer_preference && !is.null(layer_preference)) {

    # Add "none" for species not specified by the user
    missing_layer_species <- setdiff(species_names, names(layer_preference))

    if (length(missing_layer_species) > 0) {
      warning(
        "layer_preference is missing for: ",
        paste(missing_layer_species, collapse = ", "),
        ". Using 'none' for these species."
      )

      layer_preference[missing_layer_species] <- rep(
        list("none"),
        length(missing_layer_species)
      )
    }

    # Remove extra species if any
    extra_layer_species <- setdiff(names(layer_preference), species_names)

    if (length(extra_layer_species) > 0) {
      warning(
        "layer_preference contains species not found in species_roles$guild: ",
        paste(extra_layer_species, collapse = ", "),
        ". These entries will be ignored."
      )

      layer_preference <- layer_preference[setdiff(names(layer_preference), extra_layer_species)]
    }

    # Reorder and normalize
    layer_preference <- layer_preference[species_names]

    for (sp in species_names) {
      layer_preference[[sp]] <- normalize_layer_preference(layer_preference[[sp]], sp)
    }
  }

  get_named_value <- function(x, nm, default = NA_real_) {
    if (is.null(x) || is.null(names(x)) || !(nm %in% names(x))) {
      return(default)
    }
    unname(x[[nm]])
  }

  for (phase_name in names(entry_priority)) {
    phase <- entry_priority[[phase_name]]
    for (sp in phase) {
      guild <- species_roles$guild[sp]
      hosts <- species_roles$hosts[[sp]]
      fec_mean <- fecundity_mean[sp]
      fec_disp <- fecundity_dispersion[sp]
      prob_active <- get_named_value(egg_success_prob, sp, default = 1)
      if (
        use_egg_success_by_phase &&
        !is.null(egg_success_prob_by_phase) &&
        !is.null(egg_success_prob_by_phase[[sp]]) &&
        phase_name %in% names(egg_success_prob_by_phase[[sp]])
      ) {
        prob_active <- unname(egg_success_prob_by_phase[[sp]][[phase_name]])
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

      entry_distribution <- match.arg(entry_distribution, c("lognormal", "nb"))
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

          if (use_layering) {

            remaining_by_layer <- c(
              core  = df$flower_core[j]  - df$resource_use_core[j],
              mid   = df$flower_mid[j]   - df$resource_use_mid[j],
              outer = df$flower_outer[j] - df$resource_use_outer[j]
            )

            remaining_by_layer <- pmax(remaining_by_layer, 0)

            if (sum(remaining_by_layer) <= 0 || length(individual_eggs) == 0) {
              total_eggs <- 0L
            } else {

              if (
                use_layer_preference &&
                !is.null(layer_preference[[sp]]) &&
                !is_no_layer_preference(layer_preference[[sp]])
              ) {
                # Species-specific ovary-layer preference.
                assign_prob <- layer_preference[[sp]][layer_names]
                assign_prob <- assign_prob / sum(assign_prob)
              } else {
                # No species-specific ovary-layer preference.
                # Allocate attempts according to remaining ovary-layer availability.
                assign_prob <- remaining_by_layer / sum(remaining_by_layer)
              }

              assigned_layer <- sample(
                layer_names,
                size = length(individual_eggs),
                replace = TRUE,
                prob = assign_prob
              )

              total_eggs <- 0L

              for (lay in layer_names) {
                eggs_layer <- individual_eggs[assigned_layer == lay]
                available_layer <- remaining_by_layer[lay]

                eggs_layer <- cap_eggs_by_available(
                  eggs = eggs_layer,
                  available = available_layer
                )

                layer_total <- sum(eggs_layer)

                df[[paste0("resource_use_", lay)]][j] <-
                  df[[paste0("resource_use_", lay)]][j] + layer_total

                df[[paste0("eggs_", sp, "_", lay)]][j] <-
                  df[[paste0("eggs_", sp, "_", lay)]][j] + layer_total

                total_eggs <- total_eggs + layer_total
              }
            }

            df$resource_use[j] <- df$resource_use[j] + total_eggs
            eggs_vec[j] <- total_eggs

          } else {

            # Original pooled-resource logic when ovary layering is disabled.
            available <- df$flower_count[j] - df$resource_use[j]

            individual_eggs <- cap_eggs_by_available(
              eggs = individual_eggs,
              available = available
            )

            total_eggs <- sum(individual_eggs)
            df$resource_use[j] <- df$resource_use[j] + total_eggs
            eggs_vec[j] <- total_eggs
          }

        } else if (guild == "parasitoid" && length(hosts) > 0) {

          if (use_layering) {

            total_eggs <- 0L

            # Build host-by-layer availability matrix
            available_host_layer <- matrix(
              0,
              nrow = length(hosts),
              ncol = length(layer_names),
              dimnames = list(hosts, layer_names)
            )

            for (h in hosts) {
              for (lay in layer_names) {
                host_egg_col <- paste0("eggs_", h, "_", lay)
                parasitized_col <- paste0("parasitized_", h, "_", lay)

                host_eggs_layer <- if (host_egg_col %in% names(df)) {
                  df[[host_egg_col]][j]
                } else {
                  0
                }

                already_parasitized_layer <- if (parasitized_col %in% names(df)) {
                  df[[parasitized_col]][j]
                } else {
                  0
                }

                available_host_layer[h, lay] <- max(
                  host_eggs_layer - already_parasitized_layer,
                  0
                )
              }
            }

            if (sum(available_host_layer) <= 0 || length(individual_eggs) == 0) {
              total_eggs <- 0L
            } else {

              available_by_layer <- colSums(available_host_layer)

              if (
                use_layer_preference &&
                !is.null(layer_preference[[sp]]) &&
                !is_no_layer_preference(layer_preference[[sp]])
              ) {
                assign_prob_layer <- layer_preference[[sp]][layer_names]

                # Avoid assigning attempts to layers with no available host eggs.
                assign_prob_layer[available_by_layer <= 0] <- 0

                if (sum(assign_prob_layer) > 0) {
                  assign_prob_layer <- assign_prob_layer / sum(assign_prob_layer)
                } else {
                  assign_prob_layer <- available_by_layer / sum(available_by_layer)
                }

              } else {
                # No species-specific layer preference:
                # assign according to host-egg availability across layers.
                assign_prob_layer <- available_by_layer / sum(available_by_layer)
              }

              assigned_layer <- sample(
                layer_names,
                size = length(individual_eggs),
                replace = TRUE,
                prob = assign_prob_layer
              )

              for (lay in layer_names) {
                eggs_this_layer <- individual_eggs[assigned_layer == lay]

                if (length(eggs_this_layer) == 0) next

                available_hosts_this_layer <- available_host_layer[, lay, drop = TRUE]
                available_hosts_this_layer <- pmax(available_hosts_this_layer, 0)

                # pmax() may drop names, so restore host names explicitly
                names(available_hosts_this_layer) <- rownames(available_host_layer)

                valid_hosts <- names(available_hosts_this_layer)[available_hosts_this_layer > 0]

                if (length(valid_hosts) == 0) next

                if (length(valid_hosts) == 1) {
                  assigned_host <- rep(valid_hosts, length(eggs_this_layer))
                } else {
                  assigned_host <- sample(
                    valid_hosts,
                    size = length(eggs_this_layer),
                    replace = TRUE,
                    prob = available_hosts_this_layer[valid_hosts] /
                      sum(available_hosts_this_layer[valid_hosts])
                  )
                }

                for (h in hosts) {
                  eggs_host_layer <- eggs_this_layer[assigned_host == h]

                  if (length(eggs_host_layer) == 0) next

                  available_h_lay <- available_host_layer[h, lay]

                  eggs_host_layer <- cap_eggs_by_available(
                    eggs = eggs_host_layer,
                    available = available_h_lay
                  )

                  layer_host_total <- sum(eggs_host_layer)

                  if (layer_host_total <= 0) next

                  # Update parasitoid eggs in this layer
                  df[[paste0("eggs_", sp, "_", lay)]][j] <-
                    df[[paste0("eggs_", sp, "_", lay)]][j] + layer_host_total

                  # Update parasitoid-on-host-by-layer diagnostic column
                  para_host_layer_col <- paste0("eggs_", sp, "_on_", h, "_", lay)
                  if (para_host_layer_col %in% names(df)) {
                    df[[para_host_layer_col]][j] <-
                      df[[para_host_layer_col]][j] + layer_host_total
                  }

                  # Update host-layer parasitized amount
                  df[[paste0("parasitized_", h, "_", lay)]][j] <-
                    df[[paste0("parasitized_", h, "_", lay)]][j] + layer_host_total

                  # Reduce local availability so later allocations in the same loop
                  # cannot exceed the same host-layer pool.
                  available_host_layer[h, lay] <-
                    max(available_host_layer[h, lay] - layer_host_total, 0)

                  total_eggs <- total_eggs + layer_host_total
                }
              }
            }

            eggs_vec[j] <- total_eggs

          } else {

            host_eggs <- sum(sapply(hosts, function(h) {
              df[[paste0("eggs_", h)]][j]
            }))

            total_eggs <- min(sum(individual_eggs), host_eggs)
            eggs_vec[j] <- total_eggs
          }
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

    if (use_layering) {

      parasitized_cols <- paste0("parasitized_", sp, "_", layer_names)
      parasitized_cols <- parasitized_cols[parasitized_cols %in% names(df)]

      if (length(parasitized_cols) > 0) {
        parasite_total <- rowSums(df[, parasitized_cols, drop = FALSE], na.rm = TRUE)
      } else {
        parasite_total <- rep(0, nrow(df))
      }

    } else {

      parasitoids <- species_roles$parasitoid[[sp]]

      parasite_total <- if (length(parasitoids) == 0) {
        rep(0, nrow(df))
      } else {
        egg_cols <- paste0("eggs_", parasitoids)
        egg_cols <- egg_cols[egg_cols %in% names(df)]

        if (length(egg_cols) == 0) {
          rep(0, nrow(df))
        } else {
          rowSums(df[, egg_cols, drop = FALSE], na.rm = TRUE)
        }
      }
    }

    df[[paste0("emergence_", sp)]] <-
      pmax(df[[paste0("eggs_", sp)]] - parasite_total, 0)
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

  # Keep the complete output, including layer-specific columns and other
  # diagnostic variables.
  original_summary <- df

  # Build a compact user-facing summary.
  base_cols <- c(
    "fig_id",
    "fig_diameter",
    "flower_count",
    "resource_use",
    "richness_skipped"
  )

  entry_egg_cols <- unlist(lapply(species_names, function(sp) {
    c(
      paste0("entry_", sp),
      paste0("eggs_", sp)
    )
  }))

  emergence_cols <- paste0("emergence_", species_names)

  final_cols <- c(
    "unoccupied_flowers",
    "seeds",
    "failed_ovules",
    "used_flowers",
    "resource_ratio",
    "sink_strength",
    "sink_prop",
    "sink_below_min",
    "sink_above_max",
    "drop_prob",
    "is_dropped",
    "entry_matrix"
  )

  summary_cols <- c(
    base_cols,
    entry_egg_cols,
    emergence_cols,
    final_cols
  )

  # Keep only columns that exist. This makes the output robust if some optional
  # modules are disabled or future columns are added/removed.
  summary_cols <- summary_cols[summary_cols %in% names(original_summary)]

  summary <- original_summary[, summary_cols, drop = FALSE]

  if (record_individual) {
    return(list(
      summary = summary,
      original_summary = original_summary,
      individual_eggs = individual_records
    ))
  } else {
    return(list(
      summary = summary,
      original_summary = original_summary
    ))
  }
}
