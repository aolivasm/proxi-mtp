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

test_that("Sobolev diagonal launcher defines the prespecified 12 cells", {
  skip_if_not(file.exists(test_path("../../simulation/run-sobolev-diagonal.R")),
              "Repository-only simulation integration test")
  environment <- new.env(parent = globalenv())
  sys.source(
    testthat::test_path("../../simulation/run-sobolev-diagonal.R"),
    envir = environment
  )
  config <- environment$parse_sobolev_diagonal_arguments(character())
  grid <- environment$sobolev_diagonal_grid(config$sample_sizes)
  jobs <- lapply(seq_len(nrow(grid)), function(index) {
    environment$sobolev_diagonal_job_arguments(
      grid[index, , drop = FALSE], config
    )
  })

  expect_equal(nrow(grid), 12L)
  expect_true(all(grid$beta_z == -grid$beta_w))
  expect_true(all(vapply(jobs, function(x) {
    any(x == "--kernel-family=matern_sobolev") &&
      any(x == "--critical-radius-rule=matern_sobolev") &&
      any(x == "--tuning=cv_paper") &&
      any(x == "--weighted=false")
  }, logical(1L))))

  expanded <- environment$parse_sobolev_diagonal_arguments(
    "--tuning=cv_expanded"
  )
  expanded_job <- environment$sobolev_diagonal_job_arguments(
    grid[1L, , drop = FALSE], expanded
  )
  expect_true(any(expanded_job == "--tuning=cv_expanded"))

  wide <- environment$parse_sobolev_diagonal_arguments(
    "--tuning=cv_sobolev_wide"
  )
  wide_job <- environment$sobolev_diagonal_job_arguments(
    grid[1L, , drop = FALSE], wide
  )
  expect_true(any(wide_job == "--tuning=cv_sobolev_wide"))
})

test_that("wide Sobolev calibration extends every truncated tuning boundary", {
  skip_if_not(file.exists(test_path("../../simulation/run-paper-dgp-study.R")),
              "Repository-only simulation integration test")
  environment <- new.env(parent = globalenv())
  sys.source(
    testthat::test_path("../../simulation/run-paper-dgp-study.R"),
    envir = environment
  )
  config <- environment$parse_study_arguments(character())
  config$critical_radius_rule <- "matern_sobolev"
  config$kernel_family <- "matern_sobolev"
  control <- environment$cv_sobolev_wide_study_control(
    seed = 91L,
    selection_rule = "minimum",
    inner_repeats = 2L,
    config = config
  )

  expect_equal(control$lambda_h, 10^seq(-5, -1, by = 1))
  expect_equal(control$lambda_gp, 10^seq(-3, 4, by = 1))
  expect_equal(control$bandwidth_h, 2^seq(-2, 8, by = 1))
  expect_equal(control$lambda_g, control$lambda_h)
  expect_equal(control$lambda_hp, control$lambda_gp)
  expect_equal(control$bandwidth_g, control$bandwidth_h)
})
