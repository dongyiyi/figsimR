#' Convert a Parameter List into a Param Range List for Sampling
#'
#' Converts a parameter list into a list of value ranges suitable for
#' Latin Hypercube Sampling or Bayesian optimization.
#'
#' @param param_list A named list of species parameters used in the fig wasp simulator.
#'   Must include at least: \code{entry_mu}, \code{entry_size}, \code{fecundity_mean},
#'   \code{fecundity_dispersion}, and \code{egg_success_prob}. Optional entries include
#'   \code{layer_preference}, \code{parasitism_prob}, and \code{egg_success_prob_by_phase}.
#'
#' @param prop A numeric scalar (e.g., 0.3) indicating the percentage range around each
#'   parameter's base value to use for lower and upper bounds.
#'
#' @param floor_list A named list setting minimum values for each parameter type.
#'   Example values:
#'   \code{default = 0.01}
#'   \code{egg_success_prob = 0.01}
#'   \code{layer = 0.01}
#'
#' @param ceiling_list A named list setting maximum values for each parameter type.
#'   Example values:
#'   \code{default = Inf};
#'   \code{egg_success_prob = 1};
#'   \code{layer = 1};
#'   \code{parasitism_prob = 1}.
#'
#' @return A named list of numeric vectors of length 2 (lower, upper bounds), including:
#'   \code{entry_mu_*}, \code{entry_size_*}, \code{fecundity_mean_*}, \code{egg_success_prob_*}
#'   If applicable: \code{layer_core_raw_*}, \code{parasitism_prob_*}, and \code{interaction_weight}
#'
#' @details
#' This function prevents zero or NA values by assigning default values.
#' Bounds are adjusted based on proportional range \code{prop}, and clipped
#' to valid floor and ceiling values per parameter type.
#'
#' @examples
#' data(parameter_list_default)
#' param_input <- parameter_list_default
#' param_input[c("interaction_matrix", "interaction_weight")] <- NULL
#' ranges <- convert_parameter_list_to_param_ranges(param_input)
#' str(ranges[1:3])
#'
#' @export
convert_parameter_list_to_param_ranges <- function(param_list,
                                                   prop = 0.3,
                                                   floor_list = list(default = 0.01),
                                                   ceiling_list = list(
                                                     default = Inf,
                                                     egg_success_prob = 1,
                                                     layer = 1,
                                                     parasitism_prob = 1
                                                   )) {
  param_ranges <- list()
  species_list <- names(param_list$entry_mu)

  # A fallback operator: if 'a' is not NULL, return 'a'; otherwise return 'b'
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # Internal helper to compute (lower, upper) bounds around a numeric value
  get_range <- function(val, type = "default") {
    # Handle NULL, NA, or non-numeric input with a default fallback
    if (is.null(val) || all(is.na(val)) || !is.numeric(val)) {
      val <- 0.5
    }

    # Use first value if vector > 1
    if (length(val) > 1) {
      warning("get_range() received vector of length > 1. Using first element only.")
      val <- val[1]
    }

    # Prevent zero as base value (to avoid generating zero-width ranges)
    if (val == 0) val <- 0.01

    lower <- max(floor_list[[type]] %||% floor_list$default, val * (1 - prop))
    upper <- min(ceiling_list[[type]] %||% ceiling_list$default, val * (1 + prop))

    # Fix cases where upper == lower (zero-width range)
    sorted <- sort(c(lower, upper))
    if (sorted[1] == sorted[2]) {
      sorted <- c(sorted[1] * 0.9, sorted[2] * 1.1)
    }

    return(sorted)
  }

  # Main loop: generate range for each species-specific parameter
  for (sp in species_list) {
    param_ranges[[paste0("entry_mu_", sp)]] <- get_range(param_list$entry_mu[[sp]], "default")
    param_ranges[[paste0("entry_size_", sp)]] <- get_range(param_list$entry_size[[sp]], "default")
    param_ranges[[paste0("fecundity_mean_", sp)]] <- get_range(param_list$fecundity_mean[[sp]], "default")
    param_ranges[[paste0("fecundity_dispersion_", sp)]] <- get_range(param_list$fecundity_dispersion[[sp]], "default")
    param_ranges[[paste0("egg_success_prob_", sp)]] <- get_range(param_list$egg_success_prob[[sp]], "egg_success_prob")

    # Optional: if layer preference info is available (core/mid/outer)
    if (!is.null(param_list$layer_preference[[sp]])) {
      lp <- param_list$layer_preference[[sp]]
      param_ranges[[paste0("layer_core_raw_", sp)]] <- get_range(lp["core"], "layer")
      param_ranges[[paste0("layer_mid_raw_", sp)]] <- get_range(lp["mid"], "layer")
      param_ranges[[paste0("layer_outer_raw_", sp)]] <- get_range(lp["outer"], "layer")
    }

    # Optional: parasitism probability (used only for parasitoids)
    if (!is.null(param_list$parasitism_prob) && sp %in% names(param_list$parasitism_prob)) {
      val <- param_list$parasitism_prob[sp]
      param_ranges[[paste0("parasitism_prob_", sp)]] <- get_range(val, "parasitism_prob")
    }
  }

  # Handle phase-specific oviposition success probability
  if (!is.null(param_list$egg_success_prob_by_phase)) {
    for (sp in names(param_list$egg_success_prob_by_phase)) {
      for (ph in names(param_list$egg_success_prob_by_phase[[sp]])) {
        val <- param_list$egg_success_prob_by_phase[[sp]][[ph]]
        param_ranges[[paste0("egg_success_prob_by_phase_", sp, "_", ph)]] <- get_range(val, "egg_success_prob")
      }
    }
  }

  # Handle global parameter
  if (!is.null(param_list$interaction_weight)) {
    param_ranges[["interaction_weight"]] <- get_range(param_list$interaction_weight, "default")
  }

  return(param_ranges)
}
