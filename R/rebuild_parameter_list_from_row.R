#' Rebuild Parameter List from a Row of Flat Parameter Values
#'
#' This function reconstructs a nested \code{parameter_list} object by inserting flat parameter values
#' (e.g., from optimization or Latin Hypercube Sampling results) into the appropriate fields of a template list.
#'
#' @param row A named vector or single-row data.frame containing flat parameter values.
#'   Names must follow specific naming patterns (see Details).
#' @param parameter_list_template A complete simulation parameter list (e.g., \code{parameter_list_default})
#'   to be used as the base for reconstruction.
#'
#' @return A fully structured \code{parameter_list} list, where values in \code{row} are inserted
#'   into their corresponding nested positions.
#'
#' @details
#' This function supports multiple naming conventions for parameter fields:
#'
#' \itemize{
#'   \item \strong{\code{egg_success_prob_by_phase_Species_phaseX}}: assigns phase-specific egg success probability
#'     (e.g., \code{egg_success_prob_by_phase_Apocrypta_sp_phase2 = 0.6}).
#'   \item \strong{\code{layer_core_raw_Species}}, \code{layer_mid_raw_Species}, \code{layer_outer_raw_Species}:
#'     assigns spatial layer preferences for species (e.g., \code{layer_mid_raw_PollA = 0.3}).
#'   \item \strong{\code{egg_success_prob_Species}}: assigns non-phase-specific egg success probability.
#'   \item \strong{\code{prefix_species}}: general pattern for parameters such as \code{entry_mu}, \code{entry_size},
#'     \code{fecundity_mean}, etc.
#'   \item \strong{\code{prefix_subprefix_species}}: supports two-part prefixes for compatibility with list structures
#'     such as \code{entry_priority$submodel$species} (less common).
#' }
#'
#' The function parses each parameter name, extracts the relevant field and species identity,
#' and updates the \code{parameter_list_template} accordingly. Unrecognized fields are ignored.
#'
#' @examples
#' \dontrun{
#' row <- list(
#'   entry_mu_PollA = 0.5,
#'   fecundity_mean_PollA = 12,
#'   egg_success_prob_by_phase_Apocrypta_sp_phase2 = 0.6,
#'   layer_mid_raw_PollA = 0.3
#' )
#' param_list <- rebuild_parameter_list_from_row(row, parameter_list_default)
#' }
#'
#' @export

rebuild_parameter_list_from_row <- function(row, parameter_list_template) {
  param_list <- parameter_list_template

  for (param_name in names(row)) {
    val <- row[[param_name]]

    # 1. egg_success_prob_by_phase_Species_phaseX
    if (grepl("^egg_success_prob_by_phase_", param_name)) {
      m <- regmatches(param_name, regexec("^egg_success_prob_by_phase_(.+)_phase(\\d+)", param_name))[[1]]
      species <- m[2]
      phase <- paste0("phase", m[3])
      param_list$egg_success_prob_by_phase[[species]][[phase]] <- val

      # 2. layer_core_raw_Species or layer_mid_raw_Species
    } else if (grepl("^layer_(core|mid|outer)_raw_", param_name)) {
      m <- regmatches(param_name, regexec("^layer_(core|mid|outer)_raw_(.+)", param_name))[[1]]
      layer <- m[2]
      species <- m[3]
      param_list$layer_preference[[species]][[layer]] <- val

      #  3. egg_success_prob_Species
    } else if (grepl("^egg_success_prob_", param_name)) {
      species <- sub("^egg_success_prob_", "", param_name)
      param_list$egg_success_prob[[species]] <- val

      # 4. General pattern: prefix_species
    } else if (grepl("^[^_]+_[^_]+_.+", param_name)) {
      parts <- strsplit(param_name, "_")[[1]]
      prefix <- paste(parts[1:2], collapse = "_")
      species <- paste(parts[-(1:2)], collapse = "_")
      if (prefix %in% names(param_list)) {
        param_list[[prefix]][[species]] <- val
      }

      # 5. Fallback for prefix_species (single part prefix)
    } else if (grepl("^[^_]+_.+", param_name)) {
      parts <- strsplit(param_name, "_")[[1]]
      prefix <- parts[1]
      species <- paste(parts[-1], collapse = "_")
      if (prefix %in% names(param_list)) {
        param_list[[prefix]][[species]] <- val
      }
    }
  }

  return(param_list)
}
