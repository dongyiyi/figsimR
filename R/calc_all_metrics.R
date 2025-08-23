#' Calculate Alpha, Beta, and Network Metrics from Community Data
#'
#' This function computes key community-level summary metrics from a fig wasp emergence matrix.
#' It includes alpha diversity indices (richness, Shannon, Simpson, evenness),
#' beta diversity (Bray-Curtis, Jaccard), and bipartite network structure metrics (nestedness,
#' connectance, linkage density, modularity).
#'
#' Input is typically a matrix of species emergence abundances across multiple figs,
#' with optional "emergence_" prefixes in column names.
#'
#' @param df A data.frame where each row represents a fig (or sampling unit)
#'   and each column corresponds to a wasp species. Values are species abundances (counts).
#'   Column names may include or omit the prefix "emergence_".
#'
#' @return A single-row data.frame containing the following metrics:
#'   \itemize{
#'     \item \code{mean_richness}: Mean number of species per fig
#'     \item \code{mean_shannon}: Mean Shannon diversity
#'     \item \code{mean_simpson}: Mean Simpson diversity
#'     \item \code{mean_evenness}: Mean Shannon evenness (Shannon/log(richness))
#'     \item \code{mean_bray_curtis}: Mean Bray-Curtis dissimilarity (beta diversity)
#'     \item \code{mean_jaccard}: Mean Jaccard dissimilarity (beta diversity, presence-absence)
#'     \item \code{nestedness}: Nestedness based on NODF metric (from \pkg{bipartite})
#'     \item \code{connectance}: Proportion of realized links in the species-fig network
#'     \item \code{links_per_species}: Average number of links per species (linkage density)
#'     \item \code{modularity}: Modularity score from community detection (\pkg{igraph})
#'   }
#'
#' @details
#' If the input matrix contains zero-only rows, they will be excluded from
#' beta diversity and network structure calculations. Alpha diversity is computed
#' using \pkg{vegan}, beta diversity metrics are derived from Bray-Curtis and Jaccard distances,
#' and network metrics are calculated using \pkg{bipartite} and \pkg{igraph}.
#'
#' Column names with the "emergence_" prefix will be automatically stripped.
#' If fewer than two rows remain after removing zero-only rows, beta and network metrics return NA.
#'
#' @examples
#' # Simulated emergence data
#' set.seed(1)
#' df <- data.frame(
#'   emergence_A = sample(0:10, 100, replace = TRUE),
#'   emergence_B = sample(0:5, 100, replace = TRUE),
#'   emergence_C = sample(0:3, 100, replace = TRUE)
#' )
#' metrics <- calc_all_metrics(df)
#' print(metrics)
#'
#' @importFrom vegan specnumber diversity vegdist
#' @importFrom bipartite nested networklevel
#' @importFrom igraph graph_from_biadjacency_matrix cluster_fast_greedy modularity
#' @export
calc_all_metrics <- function(df) {
  if (!is.data.frame(df) || ncol(df) < 2 || nrow(df) < 2) {
    warning("Invalid input data: not a proper community matrix.")
    return(NULL)
  }

  colnames(df) <- gsub("^emergence_", "", colnames(df))

  # alpha diversity
  richness <- vegan::specnumber(df)
  shannon <- vegan::diversity(df, index = "shannon")
  simpson <- vegan::diversity(df, index = "simpson")
  evenness <- ifelse(richness >= 2, shannon / log(richness), NA_real_)

  alpha_summary <- data.frame(
    mean_richness = mean(richness, na.rm = TRUE),
    mean_shannon = mean(shannon, na.rm = TRUE),
    mean_simpson = mean(simpson, na.rm = TRUE),
    mean_evenness = mean(evenness, na.rm = TRUE)
  )


  df_nonzero <- df[rowSums(df) > 0, , drop = FALSE]
  if (nrow(df_nonzero) > 1) {
    # Bray-Curtis
    bc_dist <- vegan::vegdist(df_nonzero, method = "bray")
    mean_bc <- mean(bc_dist, na.rm = TRUE)

    # Jaccard (presence-absence)
    pa_mat <- (df_nonzero > 0) * 1
    jaccard_dist <- vegan::vegdist(pa_mat, method = "jaccard")
    mean_jaccard <- mean(jaccard_dist, na.rm = TRUE)

  } else {
    mean_bc <- NA
    mean_jaccard <- NA
    pa_mat <- (df > 0) * 1
  }

  beta_summary <- data.frame(
    mean_bray_curtis = mean_bc,
    mean_jaccard = mean_jaccard
  )

  # network
  bin_mat <- as.matrix(t(pa_mat))

  nested_val <- tryCatch({
    res <- bipartite::nested(bin_mat, method = "NODF")
    as.numeric(res["NODF"])
  }, error = function(e) {
    warning("Nestedness failed: ", conditionMessage(e))
    NA_real_
  })

  connectance_val <- tryCatch(
    as.numeric(bipartite::networklevel(bin_mat, index = "connectance")),
    error = function(e) NA_real_
  )

  links_per_species_val <- tryCatch(
    as.numeric(bipartite::networklevel(bin_mat, index = "linkage density")),
    error = function(e) NA_real_
  )

  g <- igraph::graph_from_biadjacency_matrix(bin_mat, mode = "all")
  igraph::V(g)$type <- rep(c(FALSE, TRUE), times = c(nrow(bin_mat), ncol(bin_mat)))
  clust <- igraph::cluster_fast_greedy(g)
  modularity_score <- igraph::modularity(clust)

  network_summary <- data.frame(
    nestedness = nested_val,
    connectance = connectance_val,
    links_per_species = links_per_species_val,
    modularity = modularity_score
  )

  return(cbind(alpha_summary, beta_summary, network_summary))
}
