#' Scale a numeric vector with zero-variance protection
#'
#' Standardizes a numeric vector. If the vector has zero variance, the function
#' returns a vector of zeros with the same length.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector.
#'
#' @export
scale_safe <- function(x) {
  if (stats::sd(x, na.rm = TRUE) == 0) {
    rep(0, length(x))
  } else {
    as.numeric(scale(x))
  }
}



#' Calculate emergent secondary community metrics
#'
#' Calculates four emergent secondary metrics from a fig-level community
#' matrix: mean species richness per fig, mean Pielou's evenness per fig,
#' mean pairwise Bray-Curtis dissimilarity, and multivariate dispersion.
#'
#' @param comm_mat A matrix or data frame with figs as rows and wasp species
#'   as columns. Values should be species abundances.
#'
#' @return A one-row data frame with columns \code{mean_richness},
#'   \code{mean_pielou_evenness}, \code{mean_bray_curtis}, and
#'   \code{multivariate_dispersion}.
#'
#' @export

calc_secondary_metrics <- function(comm_mat) {

  comm_mat <- as.data.frame(comm_mat)
  comm_mat[is.na(comm_mat)] <- 0

  row_totals <- rowSums(comm_mat)
  comm_nonzero <- comm_mat[row_totals > 0, , drop = FALSE]

  # 1. Species richness per fig
  richness_per_fig <- rowSums(comm_mat > 0)

  # 2. Pielou's evenness per fig
  shannon_per_fig <- apply(comm_mat, 1, function(x) {
    total <- sum(x)
    if (total <= 0) return(NA_real_)
    p <- x[x > 0] / total
    -sum(p * log(p))
  })

  pielou_per_fig <- ifelse(
    richness_per_fig > 1,
    shannon_per_fig / log(richness_per_fig),
    NA_real_
  )

  # 3. Mean Bray-Curtis dissimilarity
  if (nrow(comm_nonzero) >= 2) {
    bray_dist <- vegan::vegdist(comm_nonzero, method = "bray")
    mean_bray_curtis <- mean(as.numeric(bray_dist), na.rm = TRUE)
  } else {
    mean_bray_curtis <- NA_real_
  }

  # 4. Multivariate dispersion
  if (nrow(comm_nonzero) >= 3) {
    bray_dist <- vegan::vegdist(comm_nonzero, method = "bray")
    group <- factor(rep("community", nrow(comm_nonzero)))
    bd <- vegan::betadisper(bray_dist, group = group)
    multivariate_dispersion <- mean(bd$distances, na.rm = TRUE)
  } else {
    multivariate_dispersion <- NA_real_
  }

  data.frame(
    mean_richness = mean(richness_per_fig, na.rm = TRUE),
    mean_pielou_evenness = mean(pielou_per_fig, na.rm = TRUE),
    mean_bray_curtis = mean_bray_curtis,
    multivariate_dispersion = multivariate_dispersion
  )
}


#' Create interval summaries for secondary metrics
#'
#' Creates simulated predictive intervals and standardized observed deviations
#' for secondary community metrics.
#'
#' @param sim_df A data frame of simulated secondary metrics, usually generated
#'   from repeated simulations or bootstrap resampling.
#' @param obs_df A one-row data frame of observed secondary metrics.
#'
#' @return A data frame containing simulated means, standard deviations,
#'   80% and 95% intervals, observed values, standardized observed deviations,
#'   and interval-status labels.
#'
#' @export

make_secondary_metric_interval_df <- function(sim_df, obs_df) {

  metric_cols <- c(
    "mean_richness",
    "mean_pielou_evenness",
    "mean_bray_curtis",
    "multivariate_dispersion"
  )

  metric_labels <- c(
    mean_richness = "Species richness per fig",
    mean_pielou_evenness = "Pielou's evenness per fig",
    mean_bray_curtis = "Mean Bray-Curtis",
    multivariate_dispersion = "Multivariate dispersion"
  )

  out <- lapply(metric_cols, function(m) {

    sim_values <- sim_df[[m]]
    obs_value <- obs_df[[m]][1]

    sim_mean <- mean(sim_values, na.rm = TRUE)
    sim_sd <- sd(sim_values, na.rm = TRUE)

    data.frame(
      metric = metric_labels[[m]],
      metric_raw = m,
      sim_mean = sim_mean,
      sim_sd = sim_sd,
      lower_80 = quantile(sim_values, 0.10, na.rm = TRUE),
      upper_80 = quantile(sim_values, 0.90, na.rm = TRUE),
      lower_95 = quantile(sim_values, 0.025, na.rm = TRUE),
      upper_95 = quantile(sim_values, 0.975, na.rm = TRUE),
      observed = obs_value,
      observed_z = (obs_value - sim_mean) / sim_sd
    )
  })

  out <- dplyr::bind_rows(out)

  out$lower_80_z <- (out$lower_80 - out$sim_mean) / out$sim_sd
  out$upper_80_z <- (out$upper_80 - out$sim_mean) / out$sim_sd
  out$lower_95_z <- (out$lower_95 - out$sim_mean) / out$sim_sd
  out$upper_95_z <- (out$upper_95 - out$sim_mean) / out$sim_sd

  out$status <- dplyr::case_when(
    out$observed_z >= out$lower_80_z & out$observed_z <= out$upper_80_z ~ "Within 80% PI",
    out$observed_z >= out$lower_95_z & out$observed_z <= out$upper_95_z ~ "Within 95% PI",
    TRUE ~ "Outside 95% PI"
  )

  out
}

#' Calculate fig-level richness and Pielou's evenness
#'
#' Calculates species richness and Pielou's evenness for each fig.
#'
#' @param comm_mat A matrix or data frame with figs as rows and wasp species
#'   as columns.
#' @param dataset Character. Dataset label, such as \code{"Observed"} or
#'   \code{"Simulated"}.
#'
#' @return A data frame with columns \code{dataset}, \code{fig_id},
#'   \code{richness}, and \code{pielou_evenness}.
#'
#' @export

calc_per_fig_richness_evenness <- function(comm_mat, dataset = "Observed") {

  comm_mat <- as.data.frame(comm_mat)
  comm_mat[is.na(comm_mat)] <- 0

  richness <- rowSums(comm_mat > 0)

  shannon <- apply(comm_mat, 1, function(x) {
    total <- sum(x)
    if (total <= 0) return(NA_real_)
    p <- x[x > 0] / total
    -sum(p * log(p))
  })

  pielou <- ifelse(
    richness > 1,
    shannon / log(richness),
    NA_real_
  )

  data.frame(
    dataset = dataset,
    fig_id = seq_len(nrow(comm_mat)),
    richness = richness,
    pielou_evenness = pielou
  )
}


#' Calculate pairwise Bray-Curtis dissimilarity distribution
#'
#' Calculates pairwise Bray-Curtis dissimilarities among figs from a
#' community matrix.
#'
#' @param comm_mat A matrix or data frame with figs as rows and wasp species
#'   as columns.
#' @param dataset Character. Dataset label.
#' @param max_pairs Integer. Maximum number of pairwise distances to retain.
#'   If the number of pairwise distances exceeds this value, distances are
#'   randomly subsampled.
#'
#' @return A data frame with columns \code{dataset} and \code{bray_curtis}.
#'
#' @export

calc_bray_distribution <- function(comm_mat, dataset = "Observed",
                                   max_pairs = 20000) {

  comm_mat <- as.data.frame(comm_mat)
  comm_mat[is.na(comm_mat)] <- 0
  comm_mat <- comm_mat[rowSums(comm_mat) > 0, , drop = FALSE]

  if (nrow(comm_mat) < 2) {
    return(data.frame(dataset = dataset, bray_curtis = numeric(0)))
  }

  bray <- as.numeric(vegan::vegdist(comm_mat, method = "bray"))

  if (length(bray) > max_pairs) {
    bray <- sample(bray, max_pairs)
  }

  data.frame(
    dataset = dataset,
    bray_curtis = bray
  )
}


#' Calculate distance-to-centroid distribution
#'
#' Calculates distances from each fig-level community to the multivariate
#' centroid using Bray-Curtis dissimilarity and \code{vegan::betadisper()}.
#'
#' @param comm_mat A matrix or data frame with figs as rows and wasp species
#'   as columns.
#' @param dataset Character. Dataset label.
#'
#' @return A data frame with columns \code{dataset}, \code{fig_id}, and
#'   \code{distance_to_centroid}.
#'
#' @export

calc_dispersion_distribution <- function(comm_mat, dataset = "Observed") {

  comm_mat <- as.data.frame(comm_mat)
  comm_mat[is.na(comm_mat)] <- 0
  comm_mat <- comm_mat[rowSums(comm_mat) > 0, , drop = FALSE]

  if (nrow(comm_mat) < 3) {
    return(data.frame(dataset = dataset, distance_to_centroid = numeric(0)))
  }

  bray <- vegan::vegdist(comm_mat, method = "bray")
  group <- factor(rep("community", nrow(comm_mat)))
  bd <- vegan::betadisper(bray, group = group, type = "centroid")

  data.frame(
    dataset = dataset,
    fig_id = seq_len(nrow(comm_mat)),
    distance_to_centroid = bd$distances
  )
}


#' Collect simulated distributions of secondary metrics
#'
#' Runs repeated simulations and collects fig-level richness/evenness,
#' pairwise Bray-Curtis dissimilarities, and distance-to-centroid values
#' from sampled simulated figs.
#'
#' @param n_draws Integer. Number of repeated simulation draws.
#' @param sample_n Integer. Number of retained figs sampled from each simulation.
#' @param num_figs Integer. Number of figs simulated in each simulation run.
#' @param wasp_cols Character vector of column names corresponding to wasp
#'   emergence abundances in the simulator output.
#' @param species_names_out Character vector of species names to assign to
#'   the sampled community matrix.
#' @param simulator_func Simulation function. Default is
#'   \code{simulate_figwasp_community}.
#' @param max_bray_pairs_per_draw Integer. Maximum number of Bray-Curtis
#'   pairwise distances retained per draw.
#' @param ... Additional arguments passed to \code{simulator_func}.
#'
#' @return A list with three data frames: \code{per_fig}, \code{bray}, and
#'   \code{dispersion}.
#'
#' @export

collect_simulated_secondary_distributions <- function(
    n_draws = 500,
    sample_n = 200,
    num_figs = 1000,
    wasp_cols,
    species_names_out,
    simulator_func = simulate_figwasp_community,
    max_bray_pairs_per_draw = 5000,
    ...
) {

  per_fig_list <- vector("list", n_draws)
  bray_list <- vector("list", n_draws)
  dispersion_list <- vector("list", n_draws)

  pb <- utils::txtProgressBar(min = 0, max = n_draws, style = 3)
  on.exit(close(pb), add = TRUE)

  for (i in seq_len(n_draws)) {

    sim_result <- simulator_func(
      num_figs = num_figs,
      ...
    )

    df <- sim_result$summary
    df_valid <- df[df$is_dropped == 0, , drop = FALSE]

    if (nrow(df_valid) < sample_n) {
      utils::setTxtProgressBar(pb, i)
      next
    }

    sampled_df <- df_valid[
      sample(seq_len(nrow(df_valid)), sample_n),
      wasp_cols,
      drop = FALSE
    ]

    colnames(sampled_df) <- species_names_out

    per_fig_list[[i]] <- calc_per_fig_richness_evenness(
      comm_mat = sampled_df,
      dataset = "BCM"
    )
    per_fig_list[[i]]$draw <- i

    bray_list[[i]] <- calc_bray_distribution(
      comm_mat = sampled_df,
      dataset = "BCM",
      max_pairs = max_bray_pairs_per_draw
    )
    bray_list[[i]]$draw <- i

    dispersion_list[[i]] <- calc_dispersion_distribution(
      comm_mat = sampled_df,
      dataset = "BCM"
    )
    dispersion_list[[i]]$draw <- i

    utils::setTxtProgressBar(pb, i)
  }

  list(
    per_fig = dplyr::bind_rows(per_fig_list),
    bray = dplyr::bind_rows(bray_list),
    dispersion = dplyr::bind_rows(dispersion_list)
  )
}

#' Compare multivariate dispersion between observed and simulated communities
#'
#' Compares observed and simulated community matrices using multivariate
#' dispersion based on \code{vegan::betadisper()}.
#'
#' @param obs_comm_mat Observed community matrix with figs as rows and species
#'   as columns.
#' @param sim_comm_mat Simulated community matrix with figs as rows and species
#'   as columns.
#' @param transform Character. Transformation applied before distance
#'   calculation. Options are \code{"hellinger"} and \code{"none"}.
#' @param permutations Integer. Number of permutations used in
#'   \code{vegan::permutest()}.
#' @param seed Integer. Random seed for the permutation test.
#'
#' @return A list containing the observed-to-simulated dispersion ratio
#'   \code{rho}, permutation-test \code{p} value, group means, the
#'   \code{betadisper} object, permutation result, distances, and group labels.
#'
#' @export

# obs_comm_mat / sim_comm_mat: rows = figs, cols = species (counts)
dispersion_test <- function(obs_comm_mat, sim_comm_mat,
                            transform = c("hellinger", "none"),
                            permutations = 999, seed = 1) {
  transform <- match.arg(transform)

  # 1) Transform (Hellinger) and bind
  if (transform == "hellinger") {
    Xobs <- vegan::decostand(as.matrix(obs_comm_mat), "hellinger")
    Xsim <- vegan::decostand(as.matrix(sim_comm_mat), "hellinger")
  } else {
    Xobs <- as.matrix(obs_comm_mat)
    Xsim <- as.matrix(sim_comm_mat)
  }
  grp <- factor(c(rep("Observed", nrow(Xobs)),
                  rep("BCM",      nrow(Xsim))),
                levels = c("BCM", "Observed"))
  X <- rbind(Xobs, Xsim)

  # 2) Distances & betadisper
  D  <- dist(X, method = "euclidean")
  bd <- vegan::betadisper(D, group = grp)

  # 3) Permutation test (Pr(>F))
  set.seed(seed)
  prm <- vegan::permutest(bd, permutations = permutations)

  # 4) Means & ratio
  mu  <- tapply(bd$distances, grp, mean)  # mean distance to centroid per group
  rho <- unname(mu["Observed"] / mu["BCM"])
  p   <- prm$tab[1, "Pr(>F)"]

  list(
    rho = rho, p = p,
    means = mu,
    betadisper = bd, perm = prm,
    distances = bd$distances, group = grp   # for plotting Panel B
  )
}

#' Bootstrap a quantile
#'
#' Estimates a quantile by bootstrap resampling a numeric vector.
#'
#' @param v A numeric vector.
#' @param q Numeric. Quantile to estimate, between 0 and 1.
#' @param B Integer. Number of bootstrap replicates.
#'
#' @return A numeric vector of bootstrap quantile estimates.
#'
#' @export
boot_quantile <- function(v, q, B = 1000) {
  replicate(B, {
    stats::quantile(sample(v, replace = TRUE), q, type = 8)
  })
}

#' Build a bootstrap predictive envelope
#'
#' Computes bootstrap summaries for selected quantiles of a numeric vector.
#'
#' @param v Numeric vector of metric values.
#' @param q_levels Numeric vector of quantile levels to evaluate.
#' @param B Integer. Number of bootstrap replicates.
#'
#' @return A data frame with columns \code{q}, \code{mean}, \code{sd},
#'   \code{lo95}, and \code{hi95}.
#'
#' @export

predictive_envelope_full <- function(
    v,
    q_levels = c(0.05, 0.25, 0.50, 0.75, 0.95),
    B = 1000
) {
  purrr::map_dfr(q_levels, function(q) {
    bq <- boot_quantile(v, q, B = B)
    tibble::tibble(
      q = q,
      mean = mean(bq),
      sd = stats::sd(bq),
      lo95 = stats::quantile(bq, 0.025, type = 8),
      hi95 = stats::quantile(bq, 0.975, type = 8)
    )
  })
}



#' Build predictive envelopes for model metrics
#'
#' Computes bootstrap predictive envelopes for selected metrics from one model output.
#'
#' @param df_one_model A data frame containing model metric values.
#' @param metric_cols Character vector of metric columns to include.
#' @param q_levels Numeric vector of quantile levels used to build envelopes.
#' @param B Integer. Number of bootstrap replicates.
#'
#' @return A data frame with predictive-envelope summaries for each metric and quantile.
#'
#' @export

envelope_for_model <- function(
    df_one_model,
    metric_cols,
    q_levels = c(0.05, 0.25, 0.50, 0.75, 0.95),
    B = 1000
) {
  df_one_model |>
    dplyr::select(dplyr::all_of(metric_cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "metric",
      values_to = "val"
    ) |>
    dplyr::group_by(.data$metric) |>
    dplyr::summarise(
      env = list(predictive_envelope_full(.data$val, q_levels = q_levels, B = B)),
      .groups = "drop"
    ) |>
    tidyr::unnest(.data$env)
}


#' Calculate mean normalized exceedance outside a predictive envelope
#'
#' Computes the mean normalized exceedance of observed quantiles relative to
#' simulated predictive-envelope bounds.
#'
#' @param env_df A data frame containing predictive-envelope columns
#'   \code{metric}, \code{q}, \code{lo95}, and \code{hi95}.
#' @param obs_q_df A data frame containing observed quantiles with columns
#'   \code{metric}, \code{q}, and \code{obs_qval}.
#'
#' @return A numeric value giving the mean normalized exceedance. Values of
#'   zero indicate that observed quantiles fall within the envelope.
#'
#' @export

mean_exceed_norm <- function(env_df, obs_q_df){
  env_df |>
    dplyr::inner_join(obs_q_df, by=c("metric","q")) |>
    dplyr::mutate(
      half = (hi95 - lo95)/2,
      ctr  = (hi95 + lo95)/2,

      exceed = dplyr::case_when(
        obs_qval > hi95 ~ (obs_qval - hi95)/half,
        obs_qval < lo95 ~ (lo95 - obs_qval)/half,
        TRUE            ~ 0
      )
    ) |>
    dplyr::summarise(m = mean(exceed), .groups="drop") |>
    dplyr::pull(m)
}


#' Estimate observed-to-simulated dispersion ratio by matched resampling
#'
#' Estimates the ratio of observed to simulated multivariate dispersion by
#' repeatedly sampling simulated figs to match the number of observed figs.
#'
#' @param obs_comm Observed community matrix with figs as rows and species
#'   as columns.
#' @param sim_comm Simulated community matrix with figs as rows and species
#'   as columns.
#' @param reps Integer. Number of matched resampling replicates.
#' @param seed Integer. Random seed.
#'
#' @return A named numeric vector with the mean ratio and 95% interval
#'   bounds: \code{mean}, \code{lo}, and \code{hi}.
#'
#' @export

compute_rho_matched_fast <- function(obs_comm, sim_comm, reps = 200, seed = 42) {
  set.seed(seed)
  # 1) Hellinger
  Xobs <- vegan::decostand(as.matrix(obs_comm), "hellinger")
  Xsim <- vegan::decostand(as.matrix(sim_comm), "hellinger")

  # 2)
  n_obs <- nrow(Xobs)
  c_obs <- colMeans(Xobs)
  m_obs <- mean(sqrt(rowSums((Xobs - c_obs)^2)))

  # 3)
  rhos <- replicate(reps, {
    idx   <- sample.int(nrow(Xsim), n_obs, replace = nrow(Xsim) < n_obs)
    Ymod  <- Xsim[idx, , drop = FALSE]
    c_mod <- colMeans(Ymod)
    m_mod <- mean(sqrt(rowSums((Ymod - c_mod)^2)))
    m_obs / m_mod
  })

  c(mean = mean(rhos),
    lo   = quantile(rhos, 0.025),
    hi   = quantile(rhos, 0.975))
}





