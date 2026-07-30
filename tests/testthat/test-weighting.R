test_that("all-one weights exactly reproduce the unweighted code path", {
  data <- make_test_data()
  data$one <- 1
  policy <- list(identity = function(a) a)
  control <- minimal_control()

  unweighted <- pmtp(
    data, policy = policy, control = control
  )
  weighted <- pmtp(
    data, policy = policy, weights = "one",
    population_size = nrow(data), control = control
  )

  expect_equal(coef(weighted), coef(unweighted), tolerance = 1e-12)
  expect_equal(
    weighted$influence_function,
    unweighted$influence_function,
    tolerance = 1e-12
  )
  expect_equal(vcov(weighted), vcov(unweighted), tolerance = 1e-12)
})

test_that("paper normalization reproduces the displayed outcome matrix formula", {
  data <- make_test_data(n = 20L, seed = 31L)
  h <- as.matrix(data[c("A", "L", "W")])
  gp <- as.matrix(data[c("A", "L", "Z")])
  weights <- data$sampling_weight
  normalizer <- 125
  lambda_h <- 0.02
  lambda_gp <- 0.15
  sigma2_h <- 1.1
  sigma2_gp <- 0.9
  control <- pmtp_control_fixed(
    max_norm_h = Inf, max_norm_g = Inf, jitter = 1e-12
  )

  fit <- proximtp:::fit_outcome_bridge(
    h = h,
    gp = gp,
    y = data$Y,
    weights = weights,
    lambda_h = lambda_h,
    lambda_gp = lambda_gp,
    sigma2_h = sigma2_h,
    sigma2_gp = sigma2_gp,
    max_norm = Inf,
    control = control,
    normalizer = normalizer
  )

  k_h <- proximtp:::controlled_kernel_matrix(
    h, sigma2 = sigma2_h, control = control
  )
  k_gp <- proximtp:::controlled_kernel_matrix(
    gp, sigma2 = sigma2_gp, control = control
  )
  s <- diag(weights)
  g_gp <- solve(
    s %*% k_gp / normalizer +
      lambda_gp * diag(nrow(k_gp))
  ) / 4
  metric <- s %*% k_gp %*% g_gp %*% s
  expected <- solve(
    metric %*% k_h +
      normalizer^2 * lambda_h * diag(nrow(k_h)),
    metric %*% data$Y
  )

  expect_equal(drop(fit$coefficients), drop(expected), tolerance = 1e-7)
})

test_that("both weighted loss normalizations remain available", {
  data <- make_test_data(n = 40L, seed = 33L)
  policy <- list(identity = function(a) a)
  population_size <- 400
  common <- list(
    outer_folds = 2L,
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
    seed = 37L
  )
  hajek_control <- do.call(
    pmtp_control_fixed,
    c(common, list(weighted_loss_normalization = "hajek"))
  )
  paper_control <- do.call(
    pmtp_control_fixed,
    c(common, list(weighted_loss_normalization = "horvitz_thompson"))
  )

  hajek <- pmtp(
    data, policy = policy, weights = "sampling_weight",
    population_size = population_size, control = hajek_control
  )
  paper <- pmtp(
    data, policy = policy, weights = "sampling_weight",
    population_size = population_size, control = paper_control
  )

  expect_identical(hajek$weighted_loss_normalization, "hajek")
  expect_identical(
    paper$weighted_loss_normalization, "horvitz_thompson"
  )
  expect_false(isTRUE(all.equal(
    paper$nuisance, hajek$nuisance, tolerance = 1e-10
  )))
})

test_that("paper and Hajek normalization agree when all weights equal one", {
  data <- make_test_data(n = 36L, seed = 35L)
  data$one <- 1
  policy <- list(identity = function(a) a)
  hajek_control <- minimal_control(
    seed = 41L, weighted_loss_normalization = "hajek"
  )
  paper_control <- minimal_control(
    seed = 41L, weighted_loss_normalization = "horvitz_thompson"
  )

  hajek <- pmtp(
    data, policy = policy, weights = "one",
    population_size = nrow(data), control = hajek_control
  )
  paper <- pmtp(
    data, policy = policy, weights = "one",
    population_size = nrow(data), control = paper_control
  )

  expect_equal(coef(paper), coef(hajek), tolerance = 1e-12)
  expect_equal(paper$nuisance, hajek$nuisance, tolerance = 1e-12)
})

test_that("weighted target denominator and two-phase variance are coherent", {
  data <- make_test_data()
  fit <- pmtp(
    data,
    policy = list(identity = function(a) a),
    weights = "sampling_weight",
    target = "target",
    population_size = sum(data$sampling_weight),
    control = minimal_control()
  )

  denominator <- sum(data$sampling_weight * data$target)
  expected_estimate <- sum(data$sampling_weight * fit$contributions[, 1]) /
    denominator
  expected_tau2 <- sum(
    (data$sampling_weight * fit$influence_function[, 1])^2
  ) / fit$population_size
  expected_se <- sqrt(expected_tau2 / fit$population_size)

  expect_equal(coef(fit)[[1]], expected_estimate)
  expect_equal(fit$target_probability, denominator / fit$population_size)
  expect_equal(fit$asymptotic_variance[[1]], expected_tau2)
  expect_equal(fit$estimates$std_error[[1]], expected_se)
})

test_that("common rescaling of sampling weights leaves the fit unchanged", {
  data <- make_test_data()
  data$scaled_sampling_weight <- 11 * data$sampling_weight
  policy <- list(identity = function(a) a)
  control <- minimal_control(seed = 29L)

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
