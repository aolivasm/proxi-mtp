test_that("two-phase paper-DGP sampling is reproducible and correctly weighted", {
  data <- simulate_pmtp_dgp(500, seed = 121)
  first <- sample_pmtp_two_phase(data, seed = 122)
  second <- sample_pmtp_two_phase(data, seed = 122)

  expect_equal(first$included, second$included)
  expect_equal(first$phase_two, second$phase_two)
  expect_equal(
    first$phase_two$ipw,
    1 / first$phase_two$inclusion_probability
  )
  expect_true(all(first$inclusion_probability > 0 & first$inclusion_probability < 1))
  expect_identical(attr(first$phase_two, "spec"), attr(data, "spec"))
})

test_that("two-phase intercept calibration targets the final analysis size", {
  data <- simulate_pmtp_dgp(1500, seed = 126)
  sampled <- sample_pmtp_two_phase(
    data,
    seed = 127,
    target_sample_size = 100
  )

  expect_equal(sampled$expected_phase_two_n, 100, tolerance = 1e-8)
  expect_lt(abs(nrow(sampled$phase_two) - 100), 35)
  expect_true(is.finite(sampled$intercept))
})

test_that("all-one weights reproduce unweighted parametric estimates", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(300, spec, seed = 131)
  parameters <- pmtp_oracle_bridge_parameters(spec)
  data$one <- 1
  unweighted <- pmtp_parametric(
    data, spec, start_h = parameters$phi, start_g = parameters$eta
  )
  weighted <- pmtp_parametric(
    data, spec, weights = "one",
    start_h = parameters$phi, start_g = parameters$eta
  )

  expect_equal(coef(weighted), coef(unweighted), tolerance = 1e-12)
  expect_equal(weighted$coefficients, unweighted$coefficients, tolerance = 1e-12)
  expect_equal(
    weighted$standard_error,
    unweighted$standard_error,
    tolerance = 1e-10
  )
  expect_equal(vcov(weighted), vcov(unweighted), tolerance = 1e-10)
})

test_that("parametric two-phase covariance uses squared inverse weights", {
  spec <- pmtp_dgp_spec()
  phase_one <- simulate_pmtp_dgp(1200, spec, seed = 133)
  sampled <- sample_pmtp_two_phase(
    phase_one, seed = 134, target_sample_size = 250
  )
  data <- sampled$phase_two
  parameters <- pmtp_analytic_bridge_parameters(spec)

  fit <- pmtp_parametric(
    data, spec,
    weights = "ipw",
    population_size = nrow(phase_one),
    start_h = parameters$phi,
    start_g = parameters$eta
  )
  expected <- crossprod(data$ipw * fit$influence_function) /
    nrow(phase_one)^2

  expect_equal(vcov(fit), expected, tolerance = 1e-12)
  expect_equal(fit$population_size, nrow(phase_one))
  expect_equal(
    fit$target_probability,
    sum(data$ipw * fit$target) / nrow(phase_one)
  )
  expect_true(all(is.finite(fit$standard_error)))
})

test_that("common weight rescaling leaves parametric inference unchanged", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(400, spec, seed = 135)
  parameters <- pmtp_analytic_bridge_parameters(spec)
  data$weight <- seq(0.75, 2.25, length.out = nrow(data))
  data$scaled_weight <- 11 * data$weight

  original <- pmtp_parametric(
    data, spec,
    weights = "weight",
    population_size = sum(data$weight),
    start_h = parameters$phi,
    start_g = parameters$eta
  )
  scaled <- pmtp_parametric(
    data, spec,
    weights = "scaled_weight",
    population_size = sum(data$scaled_weight),
    start_h = parameters$phi,
    start_g = parameters$eta
  )

  expect_equal(coef(scaled), coef(original), tolerance = 1e-10)
  expect_equal(
    scaled$standard_error,
    original$standard_error,
    tolerance = 1e-8
  )
  expect_equal(vcov(scaled), vcov(original), tolerance = 1e-8)
})

test_that("parametric suite supports two-phase sandwich inference", {
  spec <- pmtp_dgp_spec()
  phase_one <- simulate_pmtp_dgp(1000, spec, seed = 136)
  sampled <- sample_pmtp_two_phase(
    phase_one, seed = 137, target_sample_size = 220
  )
  data <- sampled$phase_two
  suite <- pmtp_parametric_suite(
    data, spec,
    weights = "ipw",
    population_size = nrow(phase_one)
  )
  expected <- crossprod(data$ipw * suite$influence_function) /
    nrow(phase_one)^2

  expect_true(all(suite$converged))
  expect_equal(vcov(suite), expected, tolerance = 1e-12)
  expect_equal(suite$population_size, nrow(phase_one))
  expect_true(all(is.finite(suite$standard_error)))
  expect_true(all(suite$standard_error > 0))
})
