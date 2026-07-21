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
})
