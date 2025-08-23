# test/test-compute_alpha_diversity.R

test_that("compute_alpha_diversity returns correct structure and values", {
  # Create mock community matrix
  test_data <- data.frame(
    SampleID = paste0("S", 1:3),
    SpeciesA = c(5, 0, 3),
    SpeciesB = c(2, 4, 1),
    SpeciesC = c(0, 1, 2)
  )

  # Case 1: With sample_id_col specified
  result_with_id <- compute_alpha_diversity(test_data, sample_id_col = "SampleID")

  # Check structure
  expect_s3_class(result_with_id, "data.frame")
  expect_equal(nrow(result_with_id), 3)
  expect_true(all(c("Sample_ID", "Richness", "Shannon", "Simpson", "Evenness") %in% names(result_with_id)))

  # Check values are numeric
  expect_type(result_with_id$Richness, "double")
  expect_type(result_with_id$Shannon, "double")
  expect_type(result_with_id$Simpson, "double")
  expect_type(result_with_id$Evenness, "double")

  # Case 2: Without sample_id_col (rownames instead)
  test_data_2 <- test_data[, -1]
  rownames(test_data_2) <- paste0("R", 1:3)

  result_no_id <- compute_alpha_diversity(test_data_2)

  expect_equal(result_no_id$Sample_ID, rownames(test_data_2))
})
