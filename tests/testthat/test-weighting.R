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
