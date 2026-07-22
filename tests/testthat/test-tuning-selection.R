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

test_that("boundary warnings distinguish near ties from unsupported minima", {
  grid <- compact_grid(
    outer_lambda = c(1, 10, 100),
    inner_lambda = 1,
    outer_bandwidth = 1,
    inner_bandwidth = 1
  )
  near_tied <- transform(
    grid,
    mean_risk = c(0.10, 0.105, 0.11),
    finite_folds = 2L,
    fold1_risk = c(0.08, 0.085, 0.09),
    fold2_risk = c(0.12, 0.125, 0.13)
  )
  selected <- select_cv_row(near_tied, grid, "test")
  diagnostic <- expect_no_warning(warn_grid_boundary(
    selected, near_tied, grid, "Test bridge"
  ))
  expect_true(diagnostic$boundary[["outer_lambda_scale"]])
  expect_true(diagnostic$has_interior_alternative[["outer_lambda_scale"]])
  expect_false(diagnostic$unsupported_boundary[["outer_lambda_scale"]])

  separated <- near_tied
  separated$mean_risk <- c(0.10, 0.20, 0.21)
  separated$fold1_risk <- c(0.099, 0.19, 0.20)
  separated$fold2_risk <- c(0.101, 0.21, 0.22)
  selected <- select_cv_row(separated, grid, "test")
  expect_warning(
    diagnostic <- warn_grid_boundary(
      selected, separated, grid, "Test bridge"
    ),
    "no interior candidate"
  )
  expect_true(diagnostic$unsupported_boundary[["outer_lambda_scale"]])
})

test_that("repeated inner CV retains every fold risk", {
  data <- make_test_data(n = 40L)
  control <- pmtp_control(
    outer_folds = 2L,
    inner_folds = 2L,
    lambda_h = 1e-2,
    lambda_gp = 1,
    lambda_g = 1e-2,
    lambda_hp = 1,
    bandwidth_h = 1,
    bandwidth_gp = 1 / 2,
    bandwidth_g = 1,
    bandwidth_hp = 1 / 2,
    max_norm_h = Inf,
    max_norm_g = Inf,
    seed = 87L,
    keep_cv = TRUE,
    inner_repeats = 2L
  )

  fit <- pmtp(
    data,
    policy = list(identity = function(a) a),
    control = control
  )

  for (outer_fold_index in seq_along(fit$tuning)) {
    outer_fold <- fit$tuning[[outer_fold_index]]
    expected_training_rows <- sum(fit$outer_fold != outer_fold_index)
    components <- list(
      outer_fold$outcome,
      outer_fold$treatment[[1L]]
    )
    for (component in components) {
      expect_equal(dim(component$fold_id), c(expected_training_rows, 2L))
      expect_equal(nrow(component$risk_folds), 4L)
      expect_equal(component$risk_folds$inner_repeat, c(1L, 1L, 2L, 2L))
      expect_equal(component$risk_folds$inner_fold, c(1L, 2L, 1L, 2L))
      expect_true(all(paste0("fold", 1:4, "_risk") %in%
        names(component$results)))
      expect_equal(component$results$finite_folds, 4L)
    }
  }
})

test_that("control validates inner repeat counts", {
  expect_equal(pmtp_control()$inner_repeats, 1L)
  expect_equal(pmtp_control(inner_repeats = 3L)$inner_repeats, 3L)
  expect_error(pmtp_control(inner_repeats = 0L), "positive integer")
  expect_error(pmtp_control(inner_repeats = 1.5), "positive integer")
})
