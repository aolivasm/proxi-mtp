test_that("minimum-risk selection remains the default", {
  grid <- compact_grid(
    outer_lambda = c(1, 10),
    inner_lambda = 1,
    outer_bandwidth = c(1, 2),
    inner_bandwidth = 1
  )
  results <- transform(
    grid,
    mean_risk = c(0.10, 0.11, 0.105, 0.115),
    finite_folds = 2L,
    fold1_risk = c(0.08, 0.10, 0.095, 0.105),
    fold2_risk = c(0.12, 0.12, 0.115, 0.125)
  )

  selected <- select_cv_row(results, grid, "test")

  expect_equal(selected$grid_index, 1L)
  expect_equal(selected$minimum_grid_index, 1L)
  expect_equal(selected$mean_risk, 0.10)
  expect_equal(selected$minimum_risk_se, 0.02)
  expect_equal(selected$one_se_threshold, 0.12)
  expect_equal(selected$one_se_candidate_count, 4L)
  expect_identical(selected$selection_rule, "minimum")
})

test_that("one-SE selection prefers regularized and smooth candidates", {
  grid <- compact_grid(
    outer_lambda = c(1, 10),
    inner_lambda = 1,
    outer_bandwidth = c(1, 2),
    inner_bandwidth = 1
  )
  results <- transform(
    grid,
    mean_risk = c(0.10, 0.11, 0.105, 0.115),
    finite_folds = 2L,
    fold1_risk = c(0.08, 0.10, 0.095, 0.105),
    fold2_risk = c(0.12, 0.12, 0.115, 0.125)
  )

  selected <- select_cv_row(
    results, grid, "test", selection_rule = "one_se_regularized"
  )

  expect_equal(selected$grid_index, 4L)
  expect_equal(selected$outer_lambda_scale, 10)
  expect_equal(selected$outer_bandwidth_scale, 2)
  expect_equal(selected$selected_risk_gap, 0.015)
  expect_identical(selected$selection_rule, "one_se_regularized")
})

test_that("one-SE selection never leaves the risk threshold", {
  grid <- compact_grid(
    outer_lambda = c(1, 10),
    inner_lambda = 1,
    outer_bandwidth = c(1, 2),
    inner_bandwidth = 1
  )
  results <- transform(
    grid,
    mean_risk = c(0.10, 0.11, 0.105, 0.13),
    finite_folds = 2L,
    fold1_risk = c(0.08, 0.10, 0.095, 0.12),
    fold2_risk = c(0.12, 0.12, 0.115, 0.14)
  )

  selected <- select_cv_row(
    results, grid, "test", selection_rule = "one_se_regularized"
  )

  expect_equal(selected$grid_index, 2L)
  expect_lte(selected$mean_risk, selected$one_se_threshold)
  expect_equal(selected$outer_lambda_scale, 10)
  expect_equal(selected$outer_bandwidth_scale, 1)
})

test_that("one-SE interior selection avoids boundaries before comparing risk", {
  grid <- compact_grid(
    outer_lambda = c(1, 10, 100),
    inner_lambda = c(1, 10, 100),
    outer_bandwidth = c(1, 2, 4),
    inner_bandwidth = 1
  )
  results <- transform(
    grid,
    mean_risk = 0.101 + seq_len(nrow(grid)) / 1e5,
    finite_folds = 2L,
    fold1_risk = 0.08,
    fold2_risk = 0.12
  )
  results$mean_risk[1] <- 0.10

  selected <- select_cv_row(
    results, grid, "test", selection_rule = "one_se_interior"
  )

  expect_equal(selected$outer_lambda_scale, 10)
  expect_equal(selected$inner_lambda_scale, 10)
  expect_equal(selected$outer_bandwidth_scale, 2)
  expect_lte(selected$mean_risk, selected$one_se_threshold)
  expect_identical(selected$selection_rule, "one_se_interior")
})

test_that("control validates selection rules", {
  expect_identical(pmtp_control()$selection_rule, "minimum")
  expect_identical(
    pmtp_control(selection_rule = "one_se_regularized")$selection_rule,
    "one_se_regularized"
  )
  expect_identical(
    pmtp_control(selection_rule = "one_se_interior")$selection_rule,
    "one_se_interior"
  )
  expect_error(pmtp_control(selection_rule = "unsupported"), "arg")
})
