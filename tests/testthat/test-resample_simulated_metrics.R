test_that("resample_simulated_metrics returns expected output", {
  skip_on_cran()
  # Load default parameter template
  data("parameter_list_default", package = "figsimR")
  data("species_list", package = "figsimR")

  test_result <- resample_simulated_metrics(
    n_draws = 5,
    sample_n = 20,
    num_figs = 50,
    simulator_func = simulate_figwasp_community,
    calc_func = calc_all_metrics,
    wasp_cols = c(
      "emergence_Ceratosolen_sp",
      "emergence_Sycophaga_testacea",
      "emergence_Apocrypta_sp",
      "emergence_Sycophaga_mayri",
      "emergence_Sycophaga_agraensis",
      "emergence_Apocrypta_westwoodi"
    ),
    fecundity_mean = parameter_list_default$fecundity_mean,
    fecundity_dispersion = parameter_list_default$fecundity_dispersion,
    entry_mu = parameter_list_default$entry_mu,
    entry_size = parameter_list_default$entry_size,
    entry_priority = parameter_list_default$entry_priority,
    species_roles = parameter_list_default$species_roles,
    max_entry_table = parameter_list_default$max_entry_table,
    enable_drop = FALSE
  )


  expect_s3_class(test_result, "data.frame")
  expect_true(all(c("mean_richness", "mean_shannon", "mean_simpson", "mean_evenness") %in% colnames(test_result)))
  expect_true("source" %in% colnames(test_result))


  expect_lte(nrow(test_result), 5)
  expect_gte(nrow(test_result), 1)


  expect_false(any(apply(test_result[, c("mean_richness", "mean_shannon", "mean_simpson", "mean_evenness")], 1, function(x) all(is.na(x)))))
})
