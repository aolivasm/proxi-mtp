test_that("paper policy agrees with equation C.1 at its breakpoints", {
  spec <- pmtp_dgp_spec()
  a <- c(-2, -1.6, 0.6, 1, 2)
  expected <- c(-1.6, -1.2, 1, 1 + 0.4 / 1.4, 2)

  expect_equal(pmtp_taper_policy(a, spec), expected, tolerance = 1e-14)
  expect_equal(
    proximtp:::pmtp_policy_support(a, spec),
    c(0, 1, 1, 1, 1)
  )
})

test_that("paper DGP is reproducible and respects truncated supports", {
  spec <- pmtp_dgp_spec()
  first <- simulate_pmtp_dgp(100, spec, seed = 81)
  second <- simulate_pmtp_dgp(100, spec, seed = 81)

  expect_equal(first, second)
  expect_true(all(first$A >= spec$c & first$A <= spec$d))
  expect_true(all(abs(first$L) <= 3))
  expect_true(all(abs(first$U - spec$beta["beta1"] * first$L) <= 3))
  expect_true(all(first$qA >= spec$c + spec$delta & first$qA <= spec$d))
  expect_equal(attr(first, "sample_truth"), mean(first$p_policy))
})

test_that("analytic bridge coefficients reproduce the supplement values", {
  parameters <- pmtp_analytic_bridge_parameters(pmtp_dgp_spec())

  expect_equal(
    parameters$phi,
    c(-1.217796672, -1.5, 0.25, 0.5, -0.75),
    tolerance = 1e-8
  )
  expect_equal(parameters$eta, c(0.4, -0.08, -0.2, -0.1), tolerance = 1e-12)
  expect_true(parameters$treatment_exact)
  expect_false(parameters$outcome_exact)
  expect_identical(parameters$outcome_reference, "appendix_f1_approximation")
  expect_identical(parameters$treatment_reference, "exact_beta8_zero")
  expect_equal(
    pmtp_oracle_bridge_parameters(pmtp_dgp_spec()),
    parameters
  )
})

test_that("outcome bridge is exact in the algebraic sanity-check DGP", {
  spec <- pmtp_dgp_spec(beta10 = 0, beta12 = 0)
  data <- simulate_pmtp_dgp(200, spec, seed = 92)
  bridge <- pmtp_analytic_bridges(data$A, data$L, data$Z, data$W, spec)
  parameters <- pmtp_analytic_bridge_parameters(spec)

  expect_true(parameters$outcome_exact)
  expect_true(parameters$treatment_exact)
  expect_identical(parameters$outcome_reference, "exact_sanity_check")
  expect_equal(bridge$h, data$p_y, tolerance = 1e-13)
  expect_equal(
    pmtp_oracle_bridges(data$A, data$L, data$Z, data$W, spec),
    bridge
  )
})

test_that("exact treatment bridge recovers the target in a large diagnostic sample", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(30000, spec, seed = 103)
  analytic <- pmtp_analytic_estimates(data, spec)

  expect_lt(
    abs(analytic$conditional_mean[["DQW"]] - analytic$sample_truth), 0.01
  )
  expect_lt(
    abs(analytic$conditional_mean[["DR"]] - analytic$sample_truth), 0.01
  )
  expect_equal(pmtp_oracle_estimates(data, spec), analytic)
})

test_that("parametric estimating equations solve on the paper DGP", {
  spec <- pmtp_dgp_spec()
  parameters <- pmtp_analytic_bridge_parameters(spec)
  data <- simulate_pmtp_dgp(600, spec, seed = 20260721)
  fit <- pmtp_parametric(
    data,
    spec,
    start_h = parameters$phi,
    start_g = parameters$eta
  )

  expect_true(all(fit$converged))
  expect_lt(max(fit$residual_norm), 1e-6)
  expect_true(all(is.finite(coef(fit))))
})
