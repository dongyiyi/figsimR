# test/test-rebuild_parameter_list_from_row.R

test_that("rebuild_parameter_list_from_row correctly reconstructs nested parameter list", {
  # Minimal template list
  template <- list(
    egg_success_prob_by_phase = list(),
    egg_success_prob = list(),
    layer_preference = list(),
    fecundity_mean = list(),
    entry_mu = list()
  )

  # Example input row
  input_row <- list(
    egg_success_prob_by_phase_Apocrypta_sp_phase1 = 0.2,
    layer_core_raw_Apocrypta_sp = 0.6,
    egg_success_prob_Apocrypta_sp = 0.9,
    fecundity_mean_Apocrypta_sp = 15,
    entry_mu_Apocrypta_sp = 1.2
  )

  # Convert to named list with class "data.frame row"
  df_row <- as.data.frame(input_row)
  row <- df_row[1, ]

  # Run function
  rebuilt <- rebuild_parameter_list_from_row(row, template)

  # Check structure
  expect_type(rebuilt, "list")
  expect_true("egg_success_prob_by_phase" %in% names(rebuilt))
  expect_true("egg_success_prob" %in% names(rebuilt))
  expect_true("layer_preference" %in% names(rebuilt))
  expect_true("fecundity_mean" %in% names(rebuilt))
  expect_true("entry_mu" %in% names(rebuilt))

  # Check nested assignments
  expect_equal(rebuilt$egg_success_prob_by_phase$Apocrypta_sp$phase1, 0.2)
  expect_equal(rebuilt$layer_preference$Apocrypta_sp$core, 0.6)
  expect_equal(rebuilt$egg_success_prob$Apocrypta_sp, 0.9)
  expect_equal(rebuilt$fecundity_mean$Apocrypta_sp, 15)
  expect_equal(rebuilt$entry_mu$Apocrypta_sp, 1.2)
})
