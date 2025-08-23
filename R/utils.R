#' @importFrom magrittr %>%
NULL

#' Null coalescing operator
#' @name null_coalesce
#' @keywords internal
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b
utils::globalVariables(c(
  ".data", "Species", "Count", "Source", "Present", "richness",
  "Proportion", "FigID", "N_present", "Total_figs", "loss", "test_results",
  "lhs_best_param_list", "observed_interaction_matrix"
))
