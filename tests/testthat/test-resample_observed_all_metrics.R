test_that("resample_observed_all_metrics returns expected output", {
  set.seed(42)


  observed_df <- data.frame(
    fig_id = paste0("Fig", 1:100),
    species_A = sample(0:5, 100, replace = TRUE),
    species_B = sample(0:5, 100, replace = TRUE),
    species_C = sample(0:5, 100, replace = TRUE),
    species_D = sample(0:5, 100, replace = TRUE),
    species_E = sample(0:5, 100, replace = TRUE),
    species_F = sample(0:5, 100, replace = TRUE)
  )


  wasp_cols <- paste0("species_", LETTERS[1:6])


  metrics_df <- resample_observed_all_metrics(
    observed_df = observed_df,
    wasp_cols = wasp_cols,
    n_draws = 10,
    sample_n = 20,
    seed = 123
  )


  expect_s3_class(metrics_df, "data.frame")
  expect_equal(nrow(metrics_df), 10)
  expect_equal(metrics_df$source, rep("Observed", 10))

  expect_true(all(c(
    "mean_richness", "mean_shannon", "mean_simpson", "mean_evenness",
    "mean_bray_curtis", "mean_jaccard",
    "nestedness", "connectance", "links_per_species", "modularity"
  ) %in% colnames(metrics_df)))
})

test_that("resample_observed_all_metrics handles invalid input gracefully", {
  dummy_df <- data.frame(a = 1:10, b = 2:11)

  expect_error(resample_observed_all_metrics(
    observed_df = dummy_df,
    wasp_cols = c("a", "b", "c")
  ))
})

test_that("resample_observed_all_metrics produces reproducible output with seed", {
  wasp_df <- data.frame(
    A = sample(0:3, 50, replace = TRUE),
    B = sample(0:3, 50, replace = TRUE),
    C = sample(0:3, 50, replace = TRUE)
  )

  df1 <- resample_observed_all_metrics(wasp_df, wasp_cols = c("A", "B", "C"), n_draws = 5, sample_n = 20, seed = 777)
  df2 <- resample_observed_all_metrics(wasp_df, wasp_cols = c("A", "B", "C"), n_draws = 5, sample_n = 20, seed = 777)

  expect_equal(df1, df2)
})
