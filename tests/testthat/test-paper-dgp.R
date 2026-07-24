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
  expect_equal(
    pmtp_paper_scenario("c8", beta_z = 2, beta_w = -0.5)$truth,
    0.2564580611
  )
  expect_equal(
    pmtp_paper_scenario("c8", beta_z = 0.5, beta_w = -2)$truth,
    0.2135172919
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
  expect_identical(fit$n_sample, nrow(data))
  expect_identical(length(fit$nuisance$h0), nrow(data))
  expect_identical(length(fit$nuisance$g0[, 1L]), nrow(data))
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
  expect_true(all(is.finite(fit$standard_error)))
  expect_true(all(fit$standard_error > 0))
  expect_equal(
    unname(diag(vcov(fit))),
    unname(fit$standard_error^2),
    tolerance = 1e-12
  )
  expect_equal(
    vcov(fit),
    crossprod(fit$weights * fit$influence_function) /
      fit$population_size^2,
    tolerance = 1e-12
  )
  expect_equal(fit$inference$jacobian_rank, 12L)
  expect_lt(
    max(abs(colSums(fit$weights * fit$influence_function) /
      fit$population_size)),
    1e-5
  )

  interval <- summary(fit, conf_level = 0.90)
  expect_equal(interval$estimate, unname(coef(fit)))
  expect_equal(interval$std_error, unname(fit$standard_error))
  expect_true(all(interval$conf_low < interval$conf_high))
})

test_that("parametric outcome-bridge derivatives include the multiplier term", {
  spec <- pmtp_dgp_spec()
  parameters <- pmtp_analytic_bridge_parameters(spec)
  data <- simulate_pmtp_dgp(40, spec, seed = 20260724)
  phi <- parameters$phi
  truncated_variance <- parameters$truncated_variance

  linear_predictor <- phi[1] + phi[2] * data$A + phi[3] * data$L +
    phi[4] * data$W + phi[5] * data$A^2
  probability <- stats::plogis(linear_predictor)
  multiplier <- 1 + truncated_variance * phi[4]^2 / 2
  expected <- multiplier * probability * (1 - probability) *
    cbind(1, data$A, data$L, data$W, data$A^2)
  expected[, 4] <- expected[, 4] +
    truncated_variance * phi[4] * probability

  numerical <- finite_difference_jacobian(
    function(candidate) {
      pmtp_parametric_h_value(
        data$A, data$L, data$W, candidate, truncated_variance
      )
    },
    phi
  )

  expect_equal(unname(numerical), unname(expected), tolerance = 1e-7)
})

test_that("parametric treatment-bridge derivatives include normalization", {
  spec <- pmtp_dgp_spec()
  parameters <- pmtp_analytic_bridge_parameters(spec)
  data <- simulate_pmtp_dgp(60, spec, seed = 20260725)
  eta <- parameters$eta
  policy_v <- pmtp_policy_v(data$A, spec)
  bridge <- pmtp_parametric_g_value(
    data$A, data$L, data$Z, eta, spec
  )
  denominator <- stats::pnorm(3 - eta[3] * policy_v) -
    stats::pnorm(-3 - eta[3] * policy_v)
  normalizer_derivative <- policy_v * (
    stats::dnorm(3 - eta[3] * policy_v) -
      stats::dnorm(-3 - eta[3] * policy_v)
  ) / denominator
  expected <- bridge * cbind(
    data$A * policy_v,
    data$L * policy_v,
    data$Z * policy_v + normalizer_derivative,
    policy_v^2
  )

  numerical <- finite_difference_jacobian(
    function(candidate) {
      pmtp_parametric_g_value(
        data$A, data$L, data$Z, candidate, spec
      )
    },
    eta
  )

  expect_equal(unname(numerical), unname(expected), tolerance = 1e-6)
})

test_that("parametric solver labels a finite minimum-distance solution", {
  solution <- solve_bridge_moments(
    function(theta) theta^2 + 1,
    starts = list(0.5, -0.5),
    max_iterations = 100L,
    jacobian_function = function(theta) matrix(2 * theta, 1L, 1L)
  )

  expect_true(solution$converged)
  expect_identical(solution$solution_type, "minimum_distance")
  expect_equal(solution$coefficients, 0, tolerance = 1e-5)
  expect_equal(solution$residual_norm, 1, tolerance = 1e-8)
  expect_lt(solution$stationarity_norm, 1e-6)
})

test_that("case-2 parametric fits do not reuse case-1 calibrated starts", {
  spec <- pmtp_dgp_spec(
    beta_z = 2, beta_w = -2,
    beta8 = 0.3, beta12 = -0.3
  )
  fallback_h <- seq_len(5)
  fallback_g <- seq_len(4)

  expect_equal(
    paper_parametric_start_h(spec, fallback_h),
    fallback_h
  )
  expect_equal(
    paper_parametric_start_g_misspecified(spec, fallback_g),
    fallback_g
  )

  data <- simulate_pmtp_dgp(750, spec, seed = 20300724)
  fit <- pmtp_parametric(data, spec)
  expect_true(all(fit$converged))
  expect_true(all(is.finite(fit$standard_error)))
  expect_true(all(fit$standard_error > 0))
})

test_that("parametric suite shares the four supplement bridge fits", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(500, spec, seed = 20260726)
  standalone <- pmtp_parametric(data, spec)
  suite <- pmtp_parametric_suite(data, spec)

  expect_identical(
    names(coef(suite)),
    c(
      "OR_h_correct", "OR_h_misspecified",
      "DQW_g_correct", "DQW_g_misspecified",
      "DR_h_correct_g_correct", "DR_h_correct_g_misspecified",
      "DR_h_misspecified_g_correct",
      "DR_h_misspecified_g_misspecified"
    )
  )
  expect_true(all(suite$converged))
  expect_lt(max(suite$residual_norm), 1e-6)
  expect_length(suite$bridge_solutions$h$correct$coefficients, 5L)
  expect_length(suite$bridge_solutions$h$misspecified$coefficients, 4L)
  expect_length(suite$bridge_solutions$g$correct$coefficients, 4L)
  expect_length(suite$bridge_solutions$g$misspecified$coefficients, 4L)
  expect_equal(
    unname(standalone$estimates),
    unname(coef(suite)[c(
      "OR_h_correct", "DQW_g_correct", "DR_h_correct_g_correct"
    )]),
    tolerance = 1e-9
  )
  expect_equal(
    unname(standalone$standard_error),
    unname(suite$standard_error[c(
      "OR_h_correct", "DQW_g_correct", "DR_h_correct_g_correct"
    )]),
    tolerance = 1e-9
  )
  expect_true(all(is.finite(suite$standard_error)))
  expect_true(all(suite$standard_error > 0))
  expect_equal(
    vcov(suite),
    crossprod(suite$weights * suite$influence_function) /
      suite$population_size^2,
    tolerance = 1e-12
  )
  expect_equal(
    unname(diag(vcov(suite))),
    unname(suite$standard_error^2),
    tolerance = 1e-12
  )
})

test_that("individual misspecified parametric models match suite components", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(400, spec, seed = 20260727)
  suite <- pmtp_parametric_suite(data, spec)
  h_misspecified <- pmtp_parametric(
    data, spec,
    outcome_model = "omit_quadratic",
    treatment_model = "correct"
  )
  g_misspecified <- pmtp_parametric(
    data, spec,
    outcome_model = "correct",
    treatment_model = "constant_v"
  )

  expect_equal(
    h_misspecified$estimates[c("OR", "DR")],
    setNames(
      coef(suite)[c(
        "OR_h_misspecified", "DR_h_misspecified_g_correct"
      )],
      c("OR", "DR")
    ),
    tolerance = 1e-9
  )
  expect_equal(
    g_misspecified$estimates[c("DQW", "DR")],
    setNames(
      coef(suite)[c(
        "DQW_g_misspecified", "DR_h_correct_g_misspecified"
      )],
      c("DQW", "DR")
    ),
    tolerance = 1e-9
  )
  expect_identical(h_misspecified$models$outcome, "omit_quadratic")
  expect_identical(g_misspecified$models$treatment, "constant_v")
})
