test_that("legacy and dimension-aware critical-radius rates are correct", {
  n <- 3000
  scale <- 0.25

  expect_equal(
    proximtp:::critical_radius_squared(n, dimension = 3L, rule = "legacy_d1"),
    log(n) / n
  )
  expect_equal(
    proximtp:::actual_outer_lambda(
      scale, n, dimension = 3L, rule = "legacy_d1"
    ),
    scale * sqrt(log(n) / n)
  )
  expect_equal(
    proximtp:::actual_inner_lambda(
      scale, n, dimension = 3L, rule = "legacy_d1"
    ),
    scale * log(n) / n
  )

  dimension_rate <- log(n)^3 / n
  expect_equal(
    proximtp:::critical_radius_squared(
      n, dimension = 3L, rule = "gaussian_dimension"
    ),
    dimension_rate
  )
  expect_equal(
    proximtp:::actual_outer_lambda(
      scale, n, dimension = 3L, rule = "gaussian_dimension"
    ),
    scale * sqrt(dimension_rate)
  )
  expect_equal(
    proximtp:::actual_inner_lambda(
      scale, n, dimension = 3L, rule = "gaussian_dimension"
    ),
    scale * dimension_rate
  )
  expect_equal(
    proximtp:::actual_risk_lambda(
      scale, n, dimension = 3L, rule = "gaussian_dimension"
    ),
    scale * dimension_rate
  )
  expect_equal(
    proximtp:::critical_radius_squared(
      n, dimension = 1L, rule = "gaussian_dimension"
    ),
    log(n) / n
  )
})

test_that("critical-radius inputs and control rules are validated", {
  expect_identical(pmtp_control()$critical_radius_rule, "legacy_d1")
  expect_identical(
    pmtp_control(
      critical_radius_rule = "gaussian_dimension"
    )$critical_radius_rule,
    "gaussian_dimension"
  )
  expect_identical(
    pmtp_control_fixed(
      critical_radius_rule = "gaussian_dimension"
    )$critical_radius_rule,
    "gaussian_dimension"
  )
  expect_identical(
    pmtp_control_fast(
      critical_radius_rule = "gaussian_dimension"
    )$critical_radius_rule,
    "gaussian_dimension"
  )
  expect_error(
    pmtp_control(critical_radius_rule = "unsupported"),
    "arg"
  )
  expect_error(
    proximtp:::critical_radius_squared(1, 3, "gaussian_dimension"),
    "greater than one"
  )
  expect_error(
    proximtp:::critical_radius_squared(100, 0, "gaussian_dimension"),
    "positive integer"
  )
  expect_error(
    proximtp:::critical_radius_squared(100, 1.5, "gaussian_dimension"),
    "positive integer"
  )
})

test_that("dimension-aware fitting uses the bridge-specific dimensions", {
  data <- make_test_data(n = 40L, seed = 181L)
  policy <- list(identity = function(a) a)
  legacy <- pmtp(
    data,
    policy = policy,
    control = minimal_control(
      seed = 182L, critical_radius_rule = "legacy_d1"
    )
  )
  dimension_aware <- pmtp(
    data,
    policy = policy,
    control = minimal_control(
      seed = 182L, critical_radius_rule = "gaussian_dimension"
    )
  )

  for (fold in seq_along(dimension_aware$tuning)) {
    expect_true(is.numeric(
      dimension_aware$tuning[[fold]]$outcome$final$constraint
    ))
    expect_true(is.numeric(
      dimension_aware$tuning[[fold]]$treatment[[1]]$final$constraint
    ))
    expect_equal(
      dimension_aware$tuning[[fold]]$outcome$final$critical_radius_dimension,
      c(outer = 3, inner = 3)
    )
    expect_equal(
      dimension_aware$tuning[[fold]]$treatment[[1]]$final$
        critical_radius_dimension,
      c(outer = 3, inner = 3)
    )
  }
  expect_true(all(is.finite(coef(dimension_aware))))
  expect_false(isTRUE(all.equal(
    legacy$nuisance$g0,
    dimension_aware$nuisance$g0,
    tolerance = 1e-10
  )))
})
