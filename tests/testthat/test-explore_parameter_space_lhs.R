test_that("explore_parameter_space_lhs returns expected top-k results", {
  skip_if_not_installed("figsimR")

  # Load default parameter template
  data("parameter_list_default", package = "figsimR")
  data("species_list", package = "figsimR")

  # Use parameter_list_default as the template and base for new_param_ranges
  parameter_list_template <- parameter_list_default

  # Remove interaction-related parameters (not supported in explore version)
  new_param_ranges <- parameter_list_template
  new_param_ranges[c("interaction_matrix", "interaction_weight")] <- NULL

  # Convert all numeric vectors into range lists [min, max] for LHS sampling
  new_param_ranges <- convert_parameter_list_to_param_ranges(new_param_ranges)

  # Prepare dummy observed_summary with correct 10-column structure
  # You should replace this with a real summary in integration testing
  observed_summary <- data.frame(
    richness = 3, shannon = 0.9, simpson = 0.75, evenness = 0.6,
    bray_curtis = 0.5, jaccard = 0.4,
    nestedness = 0.3, connectance = 0.2, links_per_species = 1.5, modularity = 0.6
  )

  # Define wasp columns
  wasp_cols <- species_list

  # Run the function on reduced n_samples for faster testing
  result <- explore_parameter_space_lhs(
    new_param_ranges = new_param_ranges,
    observed_summary = observed_summary,
    parameter_list_template = parameter_list_template,
    num_figs = 5,
    n_draws = 3,
    sample_n = 2,
    wasp_cols = wasp_cols,
    n_samples = 10,
    top_k = 3,
    n_cores = 2
  )

  # --- Assertions ---
  expect_s3_class(result, "data.frame")
  expect_true(all(c("loss") %in% names(result)))
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$loss))
})
