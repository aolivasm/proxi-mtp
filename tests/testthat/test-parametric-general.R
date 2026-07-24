make_general_parametric_policy <- function(spec) {
  force(spec)
  function(a) pmtp_taper_policy(a, spec)
}

make_general_parametric_support <- function(spec) {
  force(spec)
  function(a) pmtp_policy_support(a, spec)
}

test_that("general parametric estimator solves the paper moments", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(300, spec, seed = 9101)
  fit <- pmtp_parametric_general(
    data = data,
    treatment = "A",
    outcome = "Y",
    covariates = "L",
    negative_control_treatment = "Z",
    negative_control_outcome = "W",
    policy = list(delta = make_general_parametric_policy(spec)),
    policy_support = list(make_general_parametric_support(spec)),
    outcome_link = "identity"
  )

  expect_s3_class(fit, "pmtp_parametric_general_fit")
  expect_true(fit$valid)
  expect_equal(fit$estimates$estimator, c("OR", "DQW", "DR"))
  expect_true(all(is.finite(fit$estimates$estimate)))
  expect_true(all(is.finite(fit$estimates$std_error)))
  expect_true(all(fit$estimates$std_error > 0))
  expect_true(all(fit$diagnostics$converged))
  expect_lt(max(fit$diagnostics$stationarity_norm), 1e-6)
  expect_equal(
    unname(diag(vcov(fit))),
    unname(fit$standard_error^2),
    tolerance = 1e-12
  )
  expect_equal(coef(fit), fit$estimate_vector)
  expect_equal(nrow(summary(fit)), 3L)
})

test_that("general parametric estimator supports multiple covariates and policies", {
  spec_one <- pmtp_dgp_spec(delta = 0.5)
  spec_two <- pmtp_dgp_spec(delta = 0.8)
  data <- simulate_pmtp_dgp(350, spec_one, seed = 9102)
  data$L_extra <- as.numeric(data$L > 0)
  weights <- exp(0.15 * data$L)
  population_size <- sum(weights)
  fit <- pmtp_parametric_general(
    data = data,
    treatment = "A",
    outcome = "Y",
    covariates = c("L", "L_extra"),
    negative_control_treatment = "Z",
    negative_control_outcome = "W",
    policy = list(
      delta_05 = make_general_parametric_policy(spec_one),
      delta_08 = make_general_parametric_policy(spec_two)
    ),
    policy_support = list(
      make_general_parametric_support(spec_one),
      make_general_parametric_support(spec_two)
    ),
    weights = weights,
    population_size = population_size,
    outcome_link = "identity"
  )

  expect_true(fit$weighted)
  expect_equal(nrow(fit$estimates), 6L)
  expect_equal(
    unique(fit$estimates$policy),
    c("delta_05", "delta_08")
  )
  expect_equal(dim(vcov(fit)), c(6L, 6L))
  expect_equal(
    vcov(fit),
    crossprod(fit$weights * fit$influence_function) /
      population_size^2,
    tolerance = 1e-12
  )
  expect_lt(
    max(abs(colSums(
      fit$weights * fit$influence_function
    ) / population_size)),
    1e-6
  )
  expect_equal(
    fit$inference$jacobian_rank,
    nrow(fit$inference$jacobian)
  )
  expect_true(all(fit$diagnostics$bound_fraction[
    fit$diagnostics$component == "treatment"
  ] == 0))
})

test_that("general parametric weighted results are invariant to weight units", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(260, spec, seed = 9103)
  weights <- exp(0.1 * data$L)
  arguments <- list(
    data = data,
    treatment = "A",
    outcome = "Y",
    covariates = "L",
    negative_control_treatment = "Z",
    negative_control_outcome = "W",
    policy = list(delta = make_general_parametric_policy(spec)),
    policy_support = list(make_general_parametric_support(spec)),
    outcome_link = "identity"
  )
  original <- do.call(
    pmtp_parametric_general,
    c(arguments, list(
      weights = weights,
      population_size = sum(weights)
    ))
  )
  scaled <- do.call(
    pmtp_parametric_general,
    c(arguments, list(
      weights = 7 * weights,
      population_size = 7 * sum(weights)
    ))
  )

  expect_equal(coef(original), coef(scaled), tolerance = 1e-9)
  expect_equal(vcov(original), vcov(scaled), tolerance = 1e-9)
  expect_equal(
    original$standard_error,
    scaled$standard_error,
    tolerance = 1e-9
  )
})

test_that("general parametric estimator supports overidentified moments", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(300, spec, seed = 9401)
  fit <- pmtp_parametric_general(
    data = data,
    treatment = "A",
    outcome = "Y",
    covariates = "L",
    negative_control_treatment = "Z",
    negative_control_outcome = "W",
    policy = list(delta = make_general_parametric_policy(spec)),
    policy_support = list(make_general_parametric_support(spec)),
    outcome_bridge = ~ A + L + W,
    outcome_instruments = ~ A + I(A^2) + L + Z,
    outcome_link = "identity",
    standardize = FALSE
  )

  expect_identical(
    fit$bridge_solutions$outcome$solution_type,
    "minimum_distance"
  )
  expect_gt(fit$bridge_solutions$outcome$residual_norm, 1e-6)
  expect_lt(fit$bridge_solutions$outcome$stationarity_norm, 1e-6)
  expect_true(all(is.finite(fit$estimates$std_error)))
  expect_equal(
    fit$inference$jacobian_rank,
    nrow(fit$inference$jacobian)
  )
})

test_that("general minimum-distance scores are contamination derivatives", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(120, spec, seed = 9402)
  prepared <- prepare_general_parametric_inputs(
    data = data,
    treatment = "A",
    outcome = "Y",
    covariates = "L",
    negative_control_treatment = "Z",
    negative_control_outcome = "W",
    policy = list(delta = make_general_parametric_policy(spec)),
    weights = NULL,
    target = NULL,
    policy_support = list(make_general_parametric_support(spec)),
    population_size = nrow(data),
    outcome_bridge = ~ A + L + W,
    outcome_instruments = ~ A + I(A^2) + L + Z,
    treatment_bridge = NULL,
    treatment_instruments = NULL,
    standardize = FALSE,
    log_g_bounds = c(-30, 30)
  )
  outcome_coefficients <- c(0.1, -0.2, 0.05, 0.15)
  treatment_coefficients <- c(0.1, -0.05, 0.08, -0.1)
  outcome <- general_outcome_moment_components(
    outcome_coefficients, prepared, "identity"
  )
  treatment <- general_treatment_moment_components(
    treatment_coefficients, prepared, 1L
  )
  row <- 17L
  epsilon <- 1e-6
  perturbed <- prepared
  perturbed$weights <- rep(1 - epsilon, nrow(data))
  perturbed$weights[row] <- perturbed$weights[row] +
    epsilon * nrow(data)
  perturbed_outcome <- general_outcome_moment_components(
    outcome_coefficients, perturbed, "identity"
  )
  perturbed_treatment <- general_treatment_moment_components(
    treatment_coefficients, perturbed, 1L
  )

  expect_equal(
    (perturbed_outcome$stationarity - outcome$stationarity) / epsilon,
    outcome$scores[row, ],
    tolerance = 1e-5
  )
  expect_equal(
    (perturbed_treatment$stationarity -
       treatment$stationarity) / epsilon,
    treatment$scores[row, ],
    tolerance = 1e-5
  )
})

test_that("general parametric formulas enforce identification dimensions", {
  spec <- pmtp_dgp_spec()
  data <- simulate_pmtp_dgp(80, spec, seed = 9104)

  expect_error(
    pmtp_parametric_general(
      data = data,
      treatment = "A",
      outcome = "Y",
      covariates = "L",
      negative_control_treatment = "Z",
      negative_control_outcome = "W",
      policy = make_general_parametric_policy(spec),
      policy_support = make_general_parametric_support(spec),
      outcome_bridge = ~ A + I(A^2) + L + W,
      outcome_instruments = ~ 1
    ),
    "more coefficients than outcome instruments"
  )
  expect_error(
    pmtp_parametric_general(
      data = data,
      treatment = "A",
      outcome = "Y",
      covariates = "L",
      negative_control_treatment = "Z",
      negative_control_outcome = "W",
      policy = make_general_parametric_policy(spec),
      policy_support = make_general_parametric_support(spec),
      outcome_bridge = ~ A + L + W,
      outcome_instruments = ~ A + I(A^2) + L + Z
    ),
    "Overidentified.*standardize = FALSE"
  )
  expect_error(
    pmtp_parametric_general(
      data = data,
      treatment = "A",
      outcome = "Y",
      covariates = "L",
      negative_control_treatment = "Z",
      negative_control_outcome = "W",
      policy = make_general_parametric_policy(spec),
      policy_support = make_general_parametric_support(spec),
      log_g_bounds = c(2, -2)
    ),
    "increasing finite"
  )
})
