test_that("summarize_simulated_metrics returns expected structure", {
  skip_on_cran()

  # Construct a minimal dummy all_metrics data frame
  set.seed(123)
  test_all_metrics <- data.frame(
    emergence_PollA = rpois(100, lambda = 2),
    emergence_GallerA = rpois(100, lambda = 1),
    eggs_PollA = rpois(100, lambda = 10),
    eggs_GallerA = rpois(100, lambda = 8),
    entry_PollA = rpois(100, lambda = 5),
    entry_GallerA = rpois(100, lambda = 6),
    is_dropped = sample(0:1, 100, replace = TRUE)
  )

  test_species_list <- c("PollA", "GallerA")

  result <- summarize_simulated_metrics(
    all_metrics = test_all_metrics,
    species_list = test_species_list,
    version_label = "unit-test"
  )

  # Check output is a list
  expect_type(result, "list")

  # Check the named components exist
  expect_named(result, c(
    "stats_table", "drop_summary", "richness_vector",
    "species_by_richness", "heatmap"
  ))

  # Check stats_table has correct columns
  expect_true(all(c("Species", "Entry_Mean", "Eggs_Mean", "Emergence_Mean") %in% colnames(result$stats_table)))

  # Check richness vector is numeric and of correct length
  expect_type(result$richness_vector, "double")
  expect_equal(length(result$richness_vector), nrow(test_all_metrics))

  # Check drop_summary is a table
  expect_s3_class(result$drop_summary, "table")

  # Check ggplot object returned
  expect_s3_class(result$heatmap, "gg")

  # Check species_by_richness is a data.frame and has expected columns
  expect_s3_class(result$species_by_richness, "data.frame")
  expect_true(all(c("richness", "Species", "N_present", "Total_figs", "Proportion") %in% colnames(result$species_by_richness)))
})
