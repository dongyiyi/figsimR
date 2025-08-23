test_that("convert_parameter_list_to_param_ranges works as expected", {
  # Minimal reproducible example mimicking structure of parameter_list_default
  param_list <- list(
    entry_mu = list(A = 1, B = 2),
    entry_size = list(A = 0.5, B = 0.8),
    fecundity_mean = list(A = 10, B = 12),
    fecundity_dispersion = list(A = 0.3, B = 0.2),
    egg_success_prob = list(A = 0.9, B = 0.85),
    layer_preference = list(
      A = c(core = 0.6, mid = 0.3, outer = 0.1),
      B = c(core = 0.5, mid = 0.25, outer = 0.25)
    ),
    parasitism_prob = list(A = 0.6, B = 0.7),
    egg_success_prob_by_phase = list(
      A = list(phase1 = 0.9, phase2 = 0.8),
      B = list(phase1 = 0.7)
    ),
    interaction_weight = 0.5
  )

  res <- convert_parameter_list_to_param_ranges(param_list)

  # Test output type
  expect_type(res, "list")
  expect_true(length(res) > 0)

  # Check presence of expected keys
  expected_keys <- c(
    "entry_mu_A", "entry_size_B", "fecundity_mean_A", "fecundity_dispersion_B",
    "egg_success_prob_A", "layer_core_raw_B", "parasitism_prob_A",
    "egg_success_prob_by_phase_A_phase1", "interaction_weight"
  )
  expect_true(all(expected_keys %in% names(res)))

  # Each element must be a numeric vector of length 2
  lapply(res, function(x) {
    expect_type(x, "double")
    expect_length(x, 2)
    expect_true(x[1] < x[2])
  })

  # Check known value range logic
  expect_equal(res$entry_mu_A[1], 1 * 0.7, tolerance = 1e-8)
  expect_equal(res$entry_mu_A[2], 1 * 1.3, tolerance = 1e-8)

  # Edge case: value = 0 should fallback to 0.01
  param_list$entry_mu$C <- 0
  res2 <- convert_parameter_list_to_param_ranges(param_list)
  expect_true(res2$entry_mu_C[1] > 0)
})

