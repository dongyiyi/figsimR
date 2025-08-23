test_that("calc_all_metrics returns expected output from valid input", {
  set.seed(123)

  # Construct a fake community matrix: 4 figs × 5 species
  test_df <- data.frame(
    species_A = c(2, 0, 1, 0),
    species_B = c(0, 1, 1, 0),
    species_C = c(1, 1, 0, 0),
    species_D = c(0, 0, 1, 0),
    species_E = c(3, 2, 1, 0)
  )

  # Run the function
  metrics <- calc_all_metrics(test_df)

  # Structure: must return data.frame with 1 row and expected columns
  expect_s3_class(metrics, "data.frame")
  expect_equal(nrow(metrics), 1)

  expect_true(all(c(
    "mean_richness", "mean_shannon", "mean_simpson", "mean_evenness",
    "mean_bray_curtis", "mean_jaccard",
    "nestedness", "connectance", "links_per_species", "modularity"
  ) %in% colnames(metrics)))

  # Value checks
  expect_true(metrics$mean_richness >= 1)
  expect_true(metrics$mean_shannon >= 0)
  expect_true(metrics$modularity <= 1)
})

test_that("calc_all_metrics returns NULL and warning for invalid input", {
  # Invalid input: not enough rows/columns
  expect_warning(out1 <- calc_all_metrics(data.frame(a = 1)))
  expect_null(out1)

  expect_warning(out2 <- calc_all_metrics(data.frame(a = c(1, 2))))
  expect_null(out2)

  expect_warning(out3 <- calc_all_metrics(data.frame()))
  expect_null(out3)
})
