test_that("Sobolev critical-radius rates match the Matérn native-space order", {
  n <- 1000L
  dimension <- 3L
  sobolev_l <- 4
  sobolev_order <- (sobolev_l + dimension) / 2
  expected <- n^(
    -2 * sobolev_order / (2 * sobolev_order + dimension)
  )

  expect_equal(
    proximtp:::critical_radius_squared(
      n,
      dimension = dimension,
      rule = "matern_sobolev",
      sobolev_l = sobolev_l
    ),
    expected
  )
})

test_that("control retains the Sobolev kernel configuration", {
  control <- pmtp_control_fixed(
    critical_radius_rule = "matern_sobolev",
    kernel_family = "matern_sobolev",
    matern_smoothness = 2,
    sobolev_l = 4
  )

  expect_identical(control$kernel_family, "matern_sobolev")
  expect_identical(control$critical_radius_rule, "matern_sobolev")
  expect_equal(control$matern_smoothness, 2)
  expect_equal(control$sobolev_l, 4)
})

test_that("full-rank Sobolev Nyström reproduces exact bridge predictions", {
  set.seed(820)
  n <- 24L
  h <- matrix(stats::rnorm(3L * n), n, 3L)
  g <- matrix(stats::rnorm(3L * n), n, 3L)
  y <- stats::rnorm(n)
  weights <- stats::runif(n, 0.5, 2)
  exact <- pmtp_control_fixed(
    critical_radius_rule = "matern_sobolev",
    kernel_family = "matern_sobolev",
    max_norm_h = Inf,
    max_norm_g = Inf,
    jitter = 1e-10
  )
  nystrom <- pmtp_control_fixed(
    critical_radius_rule = "matern_sobolev",
    kernel_family = "matern_sobolev",
    kernel_approximation = "nystrom",
    nystrom_rank = function(n) n,
    max_norm_h = Inf,
    max_norm_g = Inf,
    jitter = 1e-10
  )

  exact_fit <- proximtp:::fit_outcome_bridge(
    h, g, y, weights, 0.03, 0.1, 1.2, 0.8, Inf, exact
  )
  nystrom_fit <- proximtp:::fit_outcome_bridge(
    h, g, y, weights, 0.03, 0.1, 1.2, 0.8, Inf, nystrom
  )

  expect_equal(
    proximtp:::predict_outcome_bridge(nystrom_fit, h),
    proximtp:::predict_outcome_bridge(exact_fit, h),
    tolerance = 1e-8
  )
})
