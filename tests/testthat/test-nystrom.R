test_that("the default Nyström rank grows sublinearly with fold size", {
  rule <- pmtp_nystrom_rank()

  expect_equal(rule(10L), 10L)
  expect_gt(rule(600L), rule(100L))
  expect_lt(rule(600L), 600L)
  expect_equal(rule(600L), 143L)
  expect_error(pmtp_nystrom_rank(exponent = 1), "strictly between")
  expect_error(
    pmtp_control_fixed(
      kernel_approximation = "nystrom", nystrom_rank = 0
    ),
    "positive"
  )
})

test_that("full-rank Nyström features reconstruct the Gaussian kernel", {
  set.seed(120)
  arguments <- matrix(stats::rnorm(72), 24L, 3L)
  weights <- stats::runif(24L, 0.5, 3)
  control <- pmtp_control_fixed(
    kernel_approximation = "nystrom",
    nystrom_rank = function(n) n,
    max_norm_h = Inf,
    max_norm_g = Inf
  )
  feature_map <- proximtp:::fit_nystrom_map(
    arguments, sigma2 = 1.3, weights = weights, control = control
  )

  expect_equal(feature_map$requested_rank, nrow(arguments))
  expect_equal(
    tcrossprod(feature_map$training_features),
    gaussian_kernel(arguments, sigma2 = 1.3),
    tolerance = 1e-10
  )
})

test_that("full-rank Nyström reproduces both bridge fits and risks", {
  set.seed(121)
  n <- 24L
  h <- matrix(stats::rnorm(3L * n), n, 3L)
  g <- matrix(stats::rnorm(3L * n), n, 3L)
  hq <- h
  hq[, 1L] <- 0.85 * hq[, 1L]
  y <- stats::rnorm(n)
  target <- rep(c(1, 0, 1), length.out = n)
  support <- as.numeric(h[, 1L] <= stats::quantile(h[, 1L], 0.8))
  exact_control <- pmtp_control_fixed(
    max_norm_h = Inf, max_norm_g = Inf, jitter = 1e-10
  )
  nystrom_control <- pmtp_control_fixed(
    kernel_approximation = "nystrom",
    nystrom_rank = function(n) n,
    max_norm_h = Inf,
    max_norm_g = Inf,
    jitter = 1e-10
  )

  for (weights in list(rep(1, n), stats::runif(n, 0.5, 4))) {
    outcome_exact <- proximtp:::fit_outcome_bridge(
      h, g, y, weights, 0.03, 0.1, 1.2, 0.8, Inf, exact_control
    )
    outcome_nystrom <- proximtp:::fit_outcome_bridge(
      h, g, y, weights, 0.03, 0.1, 1.2, 0.8, Inf, nystrom_control
    )
    expect_equal(
      proximtp:::predict_outcome_bridge(outcome_nystrom, h),
      proximtp:::predict_outcome_bridge(outcome_exact, h),
      tolerance = 1e-8
    )

    treatment_exact <- proximtp:::fit_treatment_bridge(
      g, h, hq, weights, target, support,
      0.03, 0.1, 0.8, 1.2, Inf, exact_control
    )
    treatment_nystrom <- proximtp:::fit_treatment_bridge(
      g, h, hq, weights, target, support,
      0.03, 0.1, 0.8, 1.2, Inf, nystrom_control
    )
    treatment_exact_value <- proximtp:::predict_treatment_bridge(
      treatment_exact, g
    )
    treatment_nystrom_value <- proximtp:::predict_treatment_bridge(
      treatment_nystrom, g
    )
    expect_equal(
      treatment_nystrom_value,
      treatment_exact_value,
      tolerance = 1e-8
    )

    outcome_exact_risk <- proximtp:::outcome_validation_risk(
      y, g, weights, 0.9, 0.08, exact_control
    )
    outcome_nystrom_risk <- proximtp:::outcome_validation_risk(
      y, g, weights, 0.9, 0.08, nystrom_control
    )
    expect_equal(
      outcome_nystrom_risk, outcome_exact_risk, tolerance = 1e-8
    )

    treatment_exact_risk <- proximtp:::treatment_validation_risk(
      treatment_exact_value, h, hq, weights, target, support,
      0.9, 0.08, exact_control
    )
    treatment_nystrom_risk <- proximtp:::treatment_validation_risk(
      treatment_exact_value, h, hq, weights, target, support,
      0.9, 0.08, nystrom_control
    )
    expect_equal(
      treatment_nystrom_risk, treatment_exact_risk, tolerance = 1e-8
    )
  }
})

test_that("full-rank Nyström preserves RKHS norm constraints", {
  set.seed(123)
  n <- 22L
  h <- matrix(stats::rnorm(3L * n), n, 3L)
  g <- matrix(stats::rnorm(3L * n), n, 3L)
  y <- stats::rnorm(n)
  weights <- stats::runif(n, 0.5, 3)
  exact_control <- pmtp_control_fixed()
  nystrom_control <- pmtp_control_fixed(
    kernel_approximation = "nystrom", nystrom_rank = function(n) n
  )
  bound <- 0.05

  exact <- proximtp:::fit_outcome_bridge(
    h, g, y, weights, 0.001, 0.1, 1, 1, bound, exact_control
  )
  nystrom <- proximtp:::fit_outcome_bridge(
    h, g, y, weights, 0.001, 0.1, 1, 1, bound, nystrom_control
  )

  expect_lte(exact$norm, bound + 1e-7)
  expect_lte(nystrom$norm, bound + 1e-7)
  expect_equal(
    proximtp:::predict_outcome_bridge(nystrom, h),
    proximtp:::predict_outcome_bridge(exact, h),
    tolerance = 1e-7
  )
})

test_that("full-rank Nyström reproduces weighted and unweighted pmtp fits", {
  data <- make_test_data(n = 34L, seed = 122L)
  fixed_arguments <- list(
    outer_folds = 2L,
    lambda_h = 1e-2,
    lambda_gp = 1,
    lambda_g = 1e-2,
    lambda_hp = 1,
    max_norm_h = Inf,
    max_norm_g = Inf,
    seed = 57L
  )
  exact_control <- do.call(pmtp_control_fixed, fixed_arguments)
  nystrom_control <- do.call(
    pmtp_control_fixed,
    c(fixed_arguments, list(
      kernel_approximation = "nystrom",
      nystrom_rank = function(n) n
    ))
  )
  policy <- list(identity = function(a) a)

  for (weighted in c(FALSE, TRUE)) {
    extra_arguments <- if (weighted) {
      list(
        weights = "sampling_weight",
        population_size = sum(data$sampling_weight)
      )
    } else {
      list()
    }
    exact_fit <- do.call(
      pmtp,
      c(list(data = data, policy = policy, control = exact_control),
        extra_arguments)
    )
    nystrom_fit <- do.call(
      pmtp,
      c(list(data = data, policy = policy, control = nystrom_control),
        extra_arguments)
    )

    expect_equal(coef(nystrom_fit), coef(exact_fit), tolerance = 1e-7)
    expect_equal(
      nystrom_fit$nuisance, exact_fit$nuisance, tolerance = 1e-7
    )
    expect_equal(vcov(nystrom_fit), vcov(exact_fit), tolerance = 1e-7)
    expect_identical(
      nystrom_fit$tuning[[1L]]$outcome$final$approximation$method,
      "nystrom"
    )
  }
})

test_that("weighted Nyström fits are invariant to common weight rescaling", {
  data <- make_test_data(n = 36L, seed = 124L)
  data$scaled_sampling_weight <- 11 * data$sampling_weight
  control <- pmtp_control_fixed(
    kernel_approximation = "nystrom",
    nystrom_rank = 12L,
    nystrom_landmarks = "weighted",
    max_norm_h = Inf,
    max_norm_g = Inf,
    seed = 58L
  )
  policy <- list(identity = function(a) a)
  original <- pmtp(
    data,
    policy = policy,
    weights = "sampling_weight",
    population_size = sum(data$sampling_weight),
    control = control
  )
  scaled <- pmtp(
    data,
    policy = policy,
    weights = "scaled_sampling_weight",
    population_size = sum(data$scaled_sampling_weight),
    control = control
  )

  expect_equal(coef(scaled), coef(original), tolerance = 1e-10)
  expect_equal(scaled$nuisance, original$nuisance, tolerance = 1e-10)
  expect_equal(vcov(scaled), vcov(original), tolerance = 1e-10)
})

test_that("low-rank Nyström runs through nested weighted tuning", {
  data <- make_test_data(n = 40L, seed = 125L)
  control <- minimal_control(seed = 59L)
  control$kernel_approximation <- "nystrom"
  control$nystrom_rank <- 7L
  control$nystrom_landmarks <- "uniform"

  fit <- suppressWarnings(pmtp(
    data,
    policy = list(identity = function(a) a),
    weights = "sampling_weight",
    population_size = sum(data$sampling_weight),
    control = control
  ))

  expect_true(all(is.finite(coef(fit))))
  expect_true(all(vapply(fit$tuning, function(fold) {
    all(is.finite(fold$outcome$results$mean_risk))
  }, logical(1))))
  expect_equal(
    fit$tuning[[1L]]$outcome$final$approximation$outer$requested_rank,
    7L
  )
})
