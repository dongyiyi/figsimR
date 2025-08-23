#' Compute Alpha Diversity at the Fig Level
#'
#' Calculates sample-wise (e.g., per-fig) alpha diversity metrics, including species richness,
#' Shannon diversity, Simpson diversity, and Pielou's evenness.
#'
#' @param data A data frame where rows represent individual samples (e.g., figs)
#'   and columns represent species abundances.
#' @param sample_id_col Optional. Name of the column containing sample (fig) IDs.
#'   If provided, it is preserved in the output. If \code{NULL}, row names are used.
#'
#' @return A data frame with alpha diversity metrics for each sample (row), including:
#'   \itemize{
#'     \item \code{Sample_ID}: The sample identifier.
#'     \item \code{Richness}: Number of species (non-zero abundances).
#'     \item \code{Shannon}: Shannon diversity index.
#'     \item \code{Simpson}: Simpson diversity index.
#'     \item \code{Evenness}: Pielou's evenness, computed as \code{Shannon / log(Richness)}.
#'   }
#' All values are rounded to 3 decimal places.
#'
#' @details
#' This function is useful for fig-level community analysis, allowing quantification
#' of diversity across replicates or treatments. All input species columns must be numeric.
#' Non-numeric columns (except \code{sample_id_col}) are ignored.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   FigID = paste0("F", 1:3),
#'   SpeciesA = c(3, 0, 2),
#'   SpeciesB = c(1, 4, 0),
#'   SpeciesC = c(0, 2, 1)
#' )
#' compute_alpha_diversity(data, sample_id_col = "FigID")
#' }
#'
#' @importFrom vegan specnumber diversity
#' @export

compute_alpha_diversity <- function(data, sample_id_col = NULL) {

  if (!is.null(sample_id_col)) {
    sample_ids <- data[[sample_id_col]]
    df <- data[, sapply(data, is.numeric)]
  } else {
    sample_ids <- rownames(data)
    df <- data
  }


  richness  <- vegan::specnumber(df)
  shannon   <- vegan::diversity(df, index = "shannon")
  simpson   <- vegan::diversity(df, index = "simpson")
  evenness  <- ifelse(richness >= 2, shannon / log(richness), NA_real_)


  data.frame(
    Sample_ID = sample_ids,
    Richness  = round(richness, 3),
    Shannon   = round(shannon, 3),
    Simpson   = round(simpson, 3),
    Evenness  = round(evenness, 3)
  )
}
