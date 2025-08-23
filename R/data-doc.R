#' List of All Fig Wasp Species Observed on \emph{Ficus racemosa}
#'
#' This character vector contains the standardized names of all fig wasp species
#' observed on \emph{Ficus racemosa}, used in the fig wasp community simulator.
#' Each species name is formatted using an underscore to join the genus and species,
#' e.g., `Genus_species`. This vector defines the species order used in all parameter
#' lists (e.g., `entry_mu`, `fecundity_mean`, etc.), and it serves as the default
#' species pool for simulations.
#'
#' @format A character vector of length 6.
#' @usage data(species_list)
#'
#' @examples
#' data(species_list)
#' print(species_list)
#' # [1] "Sycophaga_testacea"   "Apocrypta_sp"
#' #     "Sycophaga_mayri"     "Ceratosolen_sp"
#' #     "Sycophaga_agraensis" "Apocrypta_westwoodi"
#'
#' @seealso \code{\link{parameter_list_default}} for parameters corresponding to each species.
"species_list"

#' Default Parameter List for Simulating Fig Wasp Communities
#'
#' This named list contains all required input parameters for running
#' \code{\link{simulate_figwasp_community}} using the six fig wasp species observed on
#' \emph{Ficus racemosa}. It includes species-specific trait values, community entry priority,
#' host–parasitoid relationships, and ecological interaction settings.
#'
#' All vectors and lists are named using standardized species names
#' from \code{\link{species_list}}.
#'
#' @format A named list with the following components:
#' \describe{
#'   \item{\code{entry_mu}}{Named numeric vector of mean number of individuals of each species entering a fig.}
#'   \item{\code{entry_size}}{Named numeric vector controlling entry variation; used as size in NB or sdlog in lognormal.}
#'   \item{\code{fecundity_mean}}{Named numeric vector of the mean number of eggs per wasp (i.e., per female individual).}
#'   \item{\code{fecundity_dispersion}}{Named numeric vector of fecundity dispersion (size parameter in negative binomial).}
#'   \item{\code{egg_success_prob}}{Named numeric vector of baseline oviposition success probabilities (0–1).}
#'   \item{\code{egg_success_prob_by_phase}}{Nested list of success probabilities for species-by-phase combinations. Overrides \code{egg_success_prob} if enabled.}
#'   \item{\code{layer_preference}}{Named list giving relative probabilities of oviposition in fig layers: \code{core}, \code{mid}, \code{outer}. Sum should be 1. Used when \code{use_layering = TRUE}.}
#'   \item{\code{max_entry_table}}{Named numeric vector of maximum number of individuals of each species allowed to enter a fig, adjusted by fig size.}
#'   \item{\code{parasitism_prob}}{Named numeric vector of supplemental parasitism probabilities (0–1), used only if \code{use_supplemental_parasitism = TRUE}.}
#'   \item{\code{entry_priority}}{Named list of phases. Each element is a character vector of species names representing entry order. Entry is processed by phase.}
#'   \item{\code{interaction_matrix}}{Square matrix of interaction coefficients among species. Used to simulate facilitation or inhibition. Default is all zero.}
#'   \item{\code{interaction_weight}}{Scalar controlling strength of interaction effects. Higher values amplify the influence of \code{interaction_matrix}.}
#'   \item{\code{species_roles}}{A list defining ecological roles and trophic relationships. Contains:
#'     \describe{
#'       \item{\code{guild}}{Named character vector assigning each species to a guild: \code{"pollinator"}, \code{"galler"}, or \code{"parasitoid"}.}
#'       \item{\code{hosts}}{Named list of host species for each parasitoid. Empty character vector if no host.}
#'       \item{\code{parasitoid}}{Named list of parasitoids attacking each host species. Empty if species is not parasitized.}
#'     }}
#' }
#'
#' @usage data(parameter_list_default)
#'
#' @examples
#' data(parameter_list_default)
#' names(parameter_list_default)
#' parameter_list_default$entry_priority
#' parameter_list_default$species_roles$guild
#'
#' @seealso \code{\link{simulate_figwasp_community}}, \code{\link{species_list}}
"parameter_list_default"

