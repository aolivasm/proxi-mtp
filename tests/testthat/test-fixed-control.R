test_that("fixed control requires scalar hyperparameters", {
  expect_error(
    pmtp_control(tune = FALSE, lambda_h = c(1e-3, 1e-2)),
    "must be scalar"
  )
  control <- pmtp_control_fixed()
  expect_false(control$tune)
  expect_equal(control$inner_folds, 1L)
})

test_that("fixed control skips inner tuning but retains outer cross-fitting", {
  data <- make_test_data(n = 40)
  fit <- pmtp(
    data,
    policy = list(identity = function(a) a),
    control = pmtp_control_fixed(seed = 91)
  )

  expect_length(unique(fit$outer_fold), 2L)
  expect_true(all(vapply(fit$tuning, function(x) x$outcome$fixed, logical(1))))
  expect_true(all(is.finite(coef(fit))))
})
