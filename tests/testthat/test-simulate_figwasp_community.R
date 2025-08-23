test_that("simulate_figwasp_community returns valid summary with default-style parameter list", {
  skip_if_not_installed("figsimR")

  data("parameter_list_default", package = "figsimR")

  result <- simulate_figwasp_community(
    num_figs = 3,
    fecundity_mean = parameter_list_default$fecundity_mean,
    fecundity_dispersion = parameter_list_default$fecundity_dispersion,
    entry_mu = parameter_list_default$entry_mu,
    entry_size = parameter_list_default$entry_size,
    entry_priority = parameter_list_default$entry_priority,
    species_roles = parameter_list_default$species_roles,
    max_entry_table = parameter_list_default$max_entry_table,
    enable_drop = TRUE,
    drop_cancels_emergence = FALSE,
    entry_distribution = "lognormal",
    interaction_matrix = parameter_list_default$interaction_matrix,
    interaction_weight = parameter_list_default$interaction_weight,
    egg_success_prob = parameter_list_default$egg_success_prob,
    egg_success_prob_by_phase = parameter_list_default$egg_success_prob_by_phase,
    layer_preference = parameter_list_default$layer_preference,
    use_layering = TRUE
  )

  expect_type(result, "list")
  expect_true("summary" %in% names(result))

  df <- result$summary

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_true("fig_id" %in% names(df))
  expect_true(any(grepl("^entry_", names(df))))
  expect_true(any(grepl("^eggs_", names(df))))
})
