# test/test-resample_metrics_from_simulation.R

test_that("resample_metrics_from_simulation works correctly with valid input", {
  # Simulate example data
  sim_data <- data.frame(
    emergence_PollA = rpois(1000, lambda = 5),
    emergence_GallerA = rpois(1000, lambda = 3),
    emergence_ParaA = rpois(1000, lambda = 2)
  )

  # Run resampling function
  result <- resample_metrics_from_simulation(sim_data, n_reps = 10, sample_size = 100, seed = 42)

  # Check output structure
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) == 10)
  expect_true(all(c("mean_richness", "mean_shannon", "mean_simpson", "mean_evenness") %in% colnames(result)))
})

test_that("resample_metrics_from_simulation fails gracefully with no emergence columns", {
  bad_data <- data.frame(foo = rnorm(100), bar = rpois(100, 2))

  expect_error(
    resample_metrics_from_simulation(bad_data),
    "No emergence_ columns found in simulated data."
  )
})
