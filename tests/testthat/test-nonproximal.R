test_that("paper learner libraries preserve reproducibility and weight support", {
  original <- pmtp_paper_learners()
  weighted <- pmtp_paper_learners(weighted = TRUE)
  expanded <- pmtp_paper_learners(weighted = TRUE, xgboost = TRUE)

  expect_length(original, 10L)
  expect_true(all(c("SL.step", "SL.step.interaction") %in% original))
  expect_false(any(c("SL.step", "SL.step.interaction") %in% weighted))
  expect_setequal(setdiff(original, c("SL.step", "SL.step.interaction")), weighted)
  expect_identical(tail(expanded, 1L), "SL.xgboost")
})

test_that("weighted point estimator is invariant to common weight scaling", {
  skip_if_not_installed("SuperLearner")
  spec <- pmtp_dgp_spec(beta7 = 0)
  data <- simulate_pmtp_dgp(140, spec, seed = 401)
  policy <- function(a) pmtp_taper_policy(a, spec)
  learners <- c("SL.glm", "SL.mean")

  first <- pmtp_nonproximal(
    data,
    policy = policy,
    weights = rep(1, nrow(data)),
    engine = "weighted_point",
    folds = 2,
    learner_folds = 2,
    learners_outcome = learners,
    estimators = c("sdr", "tmle"),
    seed = 402
  )
  second <- pmtp_nonproximal(
    data,
    policy = policy,
    weights = rep(7, nrow(data)),
    engine = "weighted_point",
    folds = 2,
    learner_folds = 2,
    learners_outcome = learners,
    estimators = c("sdr", "tmle"),
    seed = 402
  )

  expect_identical(first$weighting, "nuisance_targeting_empirical_and_variance")
  expect_equal(first$estimates$estimate, second$estimates$estimate, tolerance = 1e-10)
  expect_equal(first$estimates$std_error, second$estimates$std_error,
               tolerance = 1e-10)
  expect_true(all(is.finite(first$predictions$density_ratio)))
  expect_length(first$learner_diagnostics, 2L)
  expect_setequal(
    first$learner_diagnostics[[1L]]$outcome$learner,
    paste0(learners, "_All")
  )
  expect_true(all(is.finite(
    first$learner_diagnostics[[1L]]$treatment$coefficient
  )))
  tmle_score <- sum(
    first$weights * first$predictions$density_ratio *
      (data$Y - first$predictions$outcome_natural_tmle)
  )
  expect_lt(abs(tmle_score), 1e-7)
})

test_that("unweighted lmtp comparator runs through the package wrapper", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")
  spec <- pmtp_dgp_spec(beta7 = 0)
  data <- simulate_pmtp_dgp(120, spec, seed = 411)
  policy <- function(a) pmtp_taper_policy(a, spec)

  fit <- pmtp_nonproximal(
    data,
    policy = policy,
    estimators = "sdr",
    engine = "lmtp",
    folds = 2,
    learner_folds = 2,
    learners_outcome = c("SL.glm", "SL.mean"),
    seed = 412
  )

  expect_identical(fit$engine, "lmtp")
  expect_identical(fit$weighting, "none")
  expect_true(all(is.finite(unlist(fit$estimates[c("estimate", "std_error")]))) )
})

test_that("weighted comparator supports restricted and data-dependent policies", {
  skip_if_not_installed("SuperLearner")
  learners <- c("SL.glm", "SL.mean")

  restricted_spec <- pmtp_dgp_spec(epsilon = 0, r = 1)
  restricted_phase_one <- simulate_pmtp_dgp(
    360, restricted_spec, seed = 421
  )
  restricted_sample <- sample_pmtp_two_phase(
    restricted_phase_one, seed = 422, target_sample_size = 180
  )$phase_two
  restricted_fit <- pmtp_nonproximal(
    restricted_sample,
    policy = function(a) pmtp_taper_policy(a, restricted_spec),
    weights = "ipw",
    target = "target",
    population_size = sum(restricted_phase_one$target),
    engine = "weighted_point",
    folds = 2,
    learner_folds = 2,
    learners_outcome = learners,
    estimators = "sdr",
    seed = 423
  )

  expect_identical(
    restricted_fit$n_sample,
    sum(restricted_sample$target == 1)
  )
  expect_equal(
    restricted_fit$population_size,
    sum(restricted_phase_one$target)
  )
  expect_true(all(is.finite(unlist(
    restricted_fit$estimates[c("estimate", "std_error")]
  ))))

  covariate_spec <- pmtp_dgp_spec()
  covariate_data <- simulate_pmtp_dgp(160, covariate_spec, seed = 424)
  covariate_policy <- function(data, treatment) {
    pmtp_covariate_taper_policy(
      data[[treatment]], data$L, covariate_spec
    )
  }
  covariate_fit <- pmtp_nonproximal(
    covariate_data,
    policy = covariate_policy,
    weights = rep(1, nrow(covariate_data)),
    engine = "weighted_point",
    folds = 2,
    learner_folds = 2,
    learners_outcome = learners,
    estimators = "sdr",
    seed = 425
  )

  expect_true(all(is.finite(covariate_fit$predictions$density_ratio)))
  expect_true(all(is.finite(covariate_fit$estimates$estimate)))
})

test_that("non-proximal comparator runs under direct proxy effects", {
  skip_if_not_installed("SuperLearner")
  spec <- pmtp_dgp_spec(beta8 = 0.3, beta12 = -0.3)
  data <- simulate_pmtp_dgp(160, spec, seed = 431)

  fit <- pmtp_nonproximal(
    data,
    policy = function(a) pmtp_taper_policy(a, spec),
    weights = rep(1, nrow(data)),
    engine = "weighted_point",
    folds = 2,
    learner_folds = 2,
    learners_outcome = c("SL.glm", "SL.mean"),
    estimators = c("sdr", "tmle"),
    seed = 432
  )

  expect_true(all(is.finite(fit$estimates$estimate)))
  expect_true(all(is.finite(fit$estimates$std_error)))
})
