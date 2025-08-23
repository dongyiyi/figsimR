# test/test-calculate_metric_loss.R

test_that("calculate_metric_loss computes correct loss and structure", {
  # Simulated metrics from 3 replicate simulations
  simulated_df <- data.frame(
    Richness = c(10, 11, 9),
    Shannon = c(1.5, 1.6, 1.4),
    Evenness = c(0.5, 0.55, 0.52)
  )

  # Observed metrics (single row)
  observed_df <- data.frame(
    Richness = 10,
    Shannon = 1.5,
    Evenness = 0.51
  )

  # Run loss calculation
  result <- calculate_metric_loss(simulated_df, observed_df)

  # Check returned structure
  expect_type(result, "list")
  expect_named(result, c("metric_table", "total_loss"))
  expect_s3_class(result$metric_table, "data.frame")
  expect_true("total_loss" %in% names(result))

  # Check that metric_table has correct columns
  expect_true(all(c("Metric", "Simulated_Mean", "Observed_Value", "Squared_Diff") %in% names(result$metric_table)))

  # Check loss value is numeric
  expect_type(result$total_loss, "double")

  # Edge case: missing values
  simulated_df_missing <- simulated_df
  simulated_df_missing$Shannon[1] <- NA
  result_missing <- calculate_metric_loss(simulated_df_missing, observed_df)

  expect_type(result_missing$total_loss, "double")
})
