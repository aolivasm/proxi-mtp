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

test_that("truncated-normal sampler is stable for extreme means", {
  means <- c(-15, -10, -4, 0, 4, 10, 15)
  draw <- withr::with_seed(
    103,
    proximtp:::rtruncnorm_vector(means, -2, 2)
  )

  expect_true(all(is.finite(draw)))
  expect_true(all(draw >= -2 & draw <= 2))
  expect_gt(draw[1L], -2)
  expect_lt(draw[length(draw)], 2)
})

test_that("supplementary DGP specifications are represented directly", {
  no_confounding <- pmtp_dgp_spec(beta7 = 0)
  direct_effects <- pmtp_dgp_spec(beta8 = 0.3, beta12 = -0.3)
  restricted <- pmtp_dgp_spec(epsilon = 0, r = 1)

  expect_equal(unname(no_confounding$beta["beta7"]), 0)
  expect_equal(unname(direct_effects$beta[c("beta8", "beta12")]), c(0.3, -0.3))
  expect_equal(
    pmtp_taper_policy(c(-2, 1.6, 1.7), restricted),
    c(-1.6, 2, NA_real_)
  )
  expect_equal(
    proximtp:::pmtp_policy_v(c(-1.7, -1.6, 0, 2), restricted),
    c(0, 1, 1, 1)
  )

  data <- simulate_pmtp_dgp(200, restricted, seed = 104)
  expect_equal(data$target, as.numeric(data$A <= 1.6))
  expect_true(all(is.finite(data$p_policy[data$target == 1])))
})

test_that("paper scenario helper centralizes policies, support, and truths", {
  c7 <- pmtp_paper_scenario("c7", beta_z = 1, beta_w = -1)
  expect_identical(c7$target, "target")
  expect_equal(c7$truth, 0.2728130435)
  expect_equal(c7$policy(c(-2, 1.6, 1.8)), c(-1.6, 2, 1.8))
  expect_equal(c7$policy_support(c(-2, -1.6, 2)), c(0, 1, 1))

  c8 <- pmtp_paper_scenario("c8", beta_z = 0.75, beta_w = -0.75)
  expect_equal(c8$truth, 0.2445412003)
  expect_equal(unname(c8$spec$beta[c("beta8", "beta12")]), c(0.3, -0.3))
  expect_equal(
    pmtp_paper_scenario("c8", beta_z = 2, beta_w = -2)$truth,
    0.2184408879
  )
  expect_equal(
    pmtp_paper_scenario("c8", beta_z = 0.5, beta_w = -0.5)$truth,
    0.2487874399
  )
  expect_error(
    pmtp_paper_scenario("c8", beta_z = 0.6, beta_w = -0.6),
    "supply `truth`"
  )

  c9 <- pmtp_paper_scenario("c9")
  data <- data.frame(A = c(-1.5, -1.5, 1.5), L = c(-1, 1, -1))
  expect_equal(c9$policy(data, "A"), c(-0.9, -1.1, 1.6875))
  expect_equal(c9$policy_support(data, "A"), c(0, 1, 1))
  expect_equal(c9$truth, 0.2406921500)
})

test_that("C.7 restricted policy is finite for core estimation", {
  scenario <- pmtp_paper_scenario("c7", beta_z = 2, beta_w = -2)
  data <- simulate_pmtp_dgp(80, scenario$spec, seed = 106)
  fit <- pmtp(
    data,
    policy = list(restricted_shift = scenario$policy),
    policy_support = list(scenario$policy_support),
    target = scenario$target,
    control = pmtp_control_fixed(outer_folds = 2L, seed = 107)
  )

  expect_true(is.finite(fit$estimates$estimate[1L]))
  expect_gt(fit$target_probability, 0)
  expect_lt(fit$target_probability, 1)
})

test_that("covariate-dependent taper implements Supplement C.9", {
  spec <- pmtp_dgp_spec()
  observed <- pmtp_covariate_taper_policy(
    a = c(0, 0, 1.5, 1.5),
    l = c(-1, 1, -1, 1),
    spec = spec
  )
  expected <- c(
    0.6,
    0.4,
    1.5 + 0.6 / 1.6 * 0.5,
    1.5 + 0.4 / 1.4 * 0.5
  )
  expect_equal(observed, expected, tolerance = 1e-14)
})

test_that("core analysis preparation accepts data-dependent policies", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(80, spec, seed = 105)
  policy <- function(data, treatment) {
    pmtp_covariate_taper_policy(data[[treatment]], data$L, spec)
  }
  policy_support <- function(data, treatment) {
    delta <- ifelse(data$L < 0, 0.6, 0.4)
    as.numeric(data[[treatment]] >= spec$c + delta & data[[treatment]] <= spec$d)
  }
  policies <- proximtp:::normalize_policies(list(covariate_taper = policy))

  expect_error(
    proximtp:::build_analysis_data(
      data, "A", "Y", "L", "Z", "W",
      NULL, NULL, policies, NULL
    ),
    "Supply its `policy_support`"
  )
  prepared <- proximtp:::build_analysis_data(
    data, "A", "Y", "L", "Z", "W",
    NULL, NULL, policies, policy_support
  )
  expect_equal(
    prepared$q$covariate_taper,
    pmtp_covariate_taper_policy(data$A, data$L, spec)
  )
  expect_equal(
    prepared$policy_support$covariate_taper,
    policy_support(data, "A")
  )

  counterfactual <- pmtp_dgp_counterfactual(data, policy)
  expect_equal(counterfactual$qA, prepared$q$covariate_taper)
  expect_true(all(is.finite(counterfactual$p_policy)))
  expect_equal(counterfactual$target, rep(1, nrow(data)))
})

test_that("DGP counterfactual helper reproduces the default policy values", {
  spec <- pmtp_dgp_spec(epsilon = 0, r = 1)
  data <- simulate_pmtp_dgp(100, spec, seed = 106)
  counterfactual <- pmtp_dgp_counterfactual(
    data, function(a) pmtp_taper_policy(a, spec)
  )

  expect_equal(counterfactual$qA, data$qA)
  expect_equal(counterfactual$p_policy, data$p_policy)
  expect_equal(counterfactual$target, data$target)
  expect_equal(counterfactual$sample_truth, attr(data, "sample_truth"))
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
