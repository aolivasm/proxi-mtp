#' Specification for the paper's simulation data-generating process
#'
#' Creates the parameter specification used by the simulation study in the
#' proximal modified-treatment-policy paper. Individual parameters can be
#' changed to construct algebraic sanity checks or additional scenarios.
#'
#' @param beta_z Coefficient of `U` in the negative-control treatment `Z`.
#' @param beta_w Coefficient of `U` in the negative-control outcome `W`.
#' @param beta8 Direct effect of `Z` on treatment.
#' @param beta10 Direct effect of `U` on the outcome logit.
#' @param beta12 Direct effect of `W` on the outcome logit.
#' @param delta,epsilon,r Policy parameters from equation (C.1) of the
#'   supplement.
#' @param c,d Lower and upper bounds of the treatment distribution.
#'
#' @return A validated list of DGP and policy parameters.
#' @export
pmtp_dgp_spec <- function(beta_z = 2, beta_w = -2, beta8 = 0,
                          beta10 = -1, beta12 = 0,
                          delta = 0.4, epsilon = 1, r = 0L,
                          c = -2, d = 2) {
  values <- c(
    beta_z = beta_z, beta_w = beta_w, beta8 = beta8,
    beta10 = beta10, beta12 = beta12, delta = delta,
    epsilon = epsilon, c = c, d = d
  )
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values))) {
    stop("All DGP parameters must be finite numeric values.", call. = FALSE)
  }
  if (delta <= 0 || epsilon < 0 || c + delta >= d - epsilon) {
    stop("The policy requires delta > 0, epsilon >= 0, and c + delta < d - epsilon.",
         call. = FALSE)
  }
  if (epsilon == 0) {
    stop("The current taper-policy implementation requires epsilon > 0.",
         call. = FALSE)
  }
  if (length(r) != 1L || is.na(r) || !r %in% c(0, 1)) {
    stop("`r` must be either 0 or 1.", call. = FALSE)
  }

  structure(list(
    beta = c(
      beta1 = 0.5, beta2 = 0.2, beta3 = beta_z,
      beta4 = 0.5, beta5 = beta_w, beta6 = 0.3,
      beta7 = 1, beta8 = beta8, beta9 = 0.5,
      beta10 = beta10, beta11 = -1.5, beta12 = beta12
    ),
    mu = -1,
    gamma = -0.75,
    c = c,
    d = d,
    delta = delta,
    epsilon = epsilon,
    r = as.integer(r)
  ), class = "pmtp_dgp_spec")
}

validate_dgp_spec <- function(spec) {
  if (!inherits(spec, "pmtp_dgp_spec")) {
    stop("`spec` must be created by `pmtp_dgp_spec()`.", call. = FALSE)
  }
  spec
}

#' Tapered shift policy used in the paper
#'
#' Evaluates equation (C.1) of the supplement. Values outside the policy's
#' source set are returned as `NA`.
#'
#' @param a Numeric treatment values.
#' @param spec A specification created by [pmtp_dgp_spec()].
#'
#' @return Numeric policy-shifted treatment values.
#' @export
pmtp_taper_policy <- function(a, spec = pmtp_dgp_spec()) {
  spec <- validate_dgp_spec(spec)
  if (!is.numeric(a) || anyNA(a) || any(!is.finite(a))) {
    stop("`a` must contain finite numeric values.", call. = FALSE)
  }
  lower_end <- spec$d - spec$epsilon - spec$delta
  source_end <- spec$d - spec$r * (spec$epsilon + spec$delta)
  first <- a >= spec$c & a <= lower_end
  second <- a > lower_end & a <= source_end
  out <- rep(NA_real_, length(a))
  out[first] <- a[first] + spec$delta
  out[second] <- a[second] +
    spec$delta / (spec$delta + spec$epsilon) * (spec$d - a[second])
  out
}

pmtp_policy_target <- function(a, spec) {
  as.numeric(
    a >= spec$c &
      a <= spec$d - spec$r * (spec$epsilon + spec$delta)
  )
}

pmtp_policy_support <- function(a, spec) {
  as.numeric(a >= spec$c + spec$delta & a <= spec$d - spec$r * spec$epsilon)
}

pmtp_policy_v <- function(a, spec) {
  middle <- a >= spec$c + spec$delta & a <= spec$d - spec$epsilon
  upper <- a > spec$d - spec$epsilon & a <= spec$d - spec$r * spec$epsilon
  as.numeric(middle) + (spec$d - a) / spec$epsilon * as.numeric(upper)
}

truncated_standard_variance <- function(bound = 3) {
  1 - 2 * bound * stats::dnorm(bound) /
    (stats::pnorm(bound) - stats::pnorm(-bound))
}

rtruncnorm_vector <- function(mean, lower, upper) {
  n <- length(mean)
  if (length(lower) == 1L) lower <- rep(lower, n)
  if (length(upper) == 1L) upper <- rep(upper, n)
  lower_probability <- stats::pnorm(lower, mean = mean)
  upper_probability <- stats::pnorm(upper, mean = mean)
  probability <- stats::runif(n, lower_probability, upper_probability)
  stats::qnorm(probability, mean = mean)
}

#' Simulate the paper's proximal MTP data-generating process
#'
#' In addition to the observed variables, the returned data contain the latent
#' confounder `U`, the natural and policy-shifted outcome probabilities, and the
#' policy-shifted treatment. These extra variables are intended for simulation
#' diagnostics and should not be passed as observed analysis variables.
#'
#' @param n Number of observations.
#' @param spec A specification created by [pmtp_dgp_spec()].
#' @param seed Optional simulation seed. The caller's random-number state is
#'   preserved.
#'
#' @return A data frame containing `Y`, `A`, `L`, `Z`, `W`, `U`, `p_y`,
#'   `p_policy`, `qA`, `target`, and `policy_support`.
#' @export
simulate_pmtp_dgp <- function(n, spec = pmtp_dgp_spec(), seed = NULL) {
  spec <- validate_dgp_spec(spec)
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || n < 1 || n != as.integer(n)) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  generator <- function() {
    beta <- spec$beta
    l <- rtruncnorm_vector(rep(0, n), -3, 3)
    mean_u <- beta["beta1"] * l
    u <- rtruncnorm_vector(mean_u, mean_u - 3, mean_u + 3)
    mean_z <- beta["beta2"] * l + beta["beta3"] * u
    z <- rtruncnorm_vector(mean_z, mean_z - 3, mean_z + 3)
    mean_w <- beta["beta4"] * l + beta["beta5"] * u
    w <- rtruncnorm_vector(mean_w, mean_w - 3, mean_w + 3)
    mean_a <- beta["beta6"] * l + beta["beta7"] * u + beta["beta8"] * z
    a <- rtruncnorm_vector(mean_a, spec$c, spec$d)
    q_a <- pmtp_taper_policy(a, spec)
    target <- pmtp_policy_target(a, spec)

    natural_linear_predictor <- spec$mu + beta["beta9"] * l +
      beta["beta10"] * u + beta["beta11"] * a +
      beta["beta12"] * w + spec$gamma * a^2
    probability <- stats::plogis(natural_linear_predictor)
    policy_probability <- rep(NA_real_, n)
    policy_linear_predictor <- spec$mu + beta["beta9"] * l[target == 1] +
      beta["beta10"] * u[target == 1] + beta["beta11"] * q_a[target == 1] +
      beta["beta12"] * w[target == 1] + spec$gamma * q_a[target == 1]^2
    policy_probability[target == 1] <- stats::plogis(policy_linear_predictor)

    data.frame(
      Y = stats::rbinom(n, 1, probability),
      A = a,
      L = l,
      Z = z,
      W = w,
      U = u,
      p_y = probability,
      p_policy = policy_probability,
      qA = q_a,
      target = target,
      policy_support = pmtp_policy_support(a, spec)
    )
  }
  data <- if (is.null(seed)) generator() else withr::with_seed(seed, generator())
  attr(data, "spec") <- spec
  attr(data, "sample_truth") <-
    mean(data$p_policy[data$target == 1])
  data
}

#' Draw a two-phase sample from a simulated paper-DGP cohort
#'
#' Generates inclusion probabilities depending only on the observed outcome
#' and covariate `L`, as required by the two-phase design considered in
#' Supplement B.2. The defaults deliberately enrich outcome cases so that
#' ignoring the weights is visibly consequential in repeated simulations.
#'
#' @param data A complete phase-one cohort generated by
#'   [simulate_pmtp_dgp()].
#' @param seed Optional sampling seed.
#' @param intercept,y_coefficient,l_coefficient Coefficients of the logistic
#'   inclusion-probability model. When `target_sample_size` is supplied, the
#'   intercept is calibrated and the provided `intercept` is ignored.
#' @param target_sample_size Optional desired expected phase-two sample size.
#'   Bernoulli sampling is still used, so the realized size will fluctuate
#'   around this value.
#' @param probability_bounds Length-two vector used to bound inclusion
#'   probabilities away from zero and one.
#'
#' @return A list containing the complete cohort, observed phase-two sample,
#'   inclusion indicators and probabilities, and inclusion rate.
#' @export
sample_pmtp_two_phase <- function(
    data,
    seed = NULL,
    intercept = -1,
    y_coefficient = 1.4,
    l_coefficient = 0.25,
    target_sample_size = NULL,
    probability_bounds = c(0.005, 0.995)) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  assert_columns(data, c("Y", "L", "A", "Z", "W"), "phase-one data")
  coefficients <- c(intercept, y_coefficient, l_coefficient)
  if (anyNA(coefficients) || any(!is.finite(coefficients))) {
    stop("Inclusion-model coefficients must be finite.", call. = FALSE)
  }
  if (!is.numeric(probability_bounds) || length(probability_bounds) != 2L ||
      anyNA(probability_bounds) || probability_bounds[1] <= 0 ||
      probability_bounds[2] >= 1 || probability_bounds[1] >= probability_bounds[2]) {
    stop("`probability_bounds` must be increasing values strictly inside (0, 1).",
         call. = FALSE)
  }
  bounded_probability <- function(candidate_intercept) {
    probability <- stats::plogis(
      candidate_intercept + y_coefficient * data$Y + l_coefficient * data$L
    )
    pmin(pmax(probability, probability_bounds[1]), probability_bounds[2])
  }
  if (!is.null(target_sample_size)) {
    if (length(target_sample_size) != 1L || !is.numeric(target_sample_size) ||
        is.na(target_sample_size) || target_sample_size <= 0 ||
        target_sample_size >= nrow(data)) {
      stop("`target_sample_size` must be between zero and the phase-one size.",
           call. = FALSE)
    }
    feasible <- nrow(data) * probability_bounds
    if (target_sample_size <= feasible[1] || target_sample_size >= feasible[2]) {
      stop("`target_sample_size` is incompatible with `probability_bounds`.",
           call. = FALSE)
    }
    objective <- function(candidate_intercept) {
      sum(bounded_probability(candidate_intercept)) - target_sample_size
    }
    intercept <- stats::uniroot(objective, c(-30, 30), tol = 1e-12)$root
  }
  inclusion_probability <- bounded_probability(intercept)
  sampler <- function() {
    stats::rbinom(nrow(data), 1, inclusion_probability) == 1
  }
  included <- if (is.null(seed)) sampler() else withr::with_seed(seed, sampler())
  if (!any(included)) {
    stop("The two-phase sample was empty; use a different seed or inclusion model.",
         call. = FALSE)
  }
  phase_two <- data[included, , drop = FALSE]
  phase_two$inclusion_probability <- inclusion_probability[included]
  phase_two$ipw <- 1 / inclusion_probability[included]
  attr(phase_two, "spec") <- attr(data, "spec")
  list(
    phase_one = data,
    phase_two = phase_two,
    included = included,
    inclusion_probability = inclusion_probability,
    inclusion_rate = mean(included),
    intercept = intercept,
    expected_phase_two_n = sum(inclusion_probability)
  )
}

#' Analytic bridge-reference parameters for the paper's DGP
#'
#' The treatment-bridge parameters are exact when `beta8 = 0` and are the
#' Appendix F.2 approximation otherwise. The outcome-bridge parameters are the
#' Appendix F.1 approximation; they are exact in the simple sanity-check case
#' `beta10 = beta12 = 0`.
#'
#' @param spec A specification created by [pmtp_dgp_spec()].
#'
#' @return A list containing `phi`, `eta`, the truncated-normal variance,
#'   indicators describing when the bridges are exact, and explicit reference
#'   labels. In particular, `phi` must not be interpreted as an oracle outcome
#'   bridge unless `outcome_exact` is `TRUE`.
#' @export
pmtp_analytic_bridge_parameters <- function(spec = pmtp_dgp_spec()) {
  spec <- validate_dgp_spec(spec)
  beta <- spec$beta
  truncated_variance <- truncated_standard_variance()
  if (abs(beta["beta5"]) <= sqrt(.Machine$double.eps)) {
    stop("The outcome-bridge approximation requires nonzero beta5.", call. = FALSE)
  }
  denominator <- beta["beta3"] - beta["beta7"] * beta["beta8"]
  if (abs(denominator) <= sqrt(.Machine$double.eps)) {
    stop("The treatment-bridge formula has beta3 - beta7 * beta8 equal to zero.",
         call. = FALSE)
  }
  ratio <- beta["beta12"] + beta["beta10"] / beta["beta5"]
  phi <- c(
    phi0 = spec$mu + log1p(truncated_variance * beta["beta12"]^2 / 2) -
      log1p(truncated_variance * ratio^2),
    phi1 = beta["beta11"],
    phi2 = beta["beta9"] - beta["beta4"] * beta["beta10"] / beta["beta5"],
    phi3 = ratio,
    phi4 = spec$gamma
  )
  eta <- c(
    eta0 = spec$delta * beta["beta3"] / denominator,
    eta1 = -spec$delta *
      (beta["beta3"] * beta["beta6"] - beta["beta2"] * beta["beta7"]) /
      denominator,
    eta2 = -spec$delta *
      (beta["beta7"] + beta["beta3"] * beta["beta8"]) / denominator,
    eta3 = -spec$delta^2 * (beta["beta3"]^2 + beta["beta7"]^2) /
      (2 * denominator^2)
  )
  outcome_exact <- isTRUE(all.equal(unname(beta["beta10"]), 0)) &&
    isTRUE(all.equal(unname(beta["beta12"]), 0))
  treatment_exact <- isTRUE(all.equal(unname(beta["beta8"]), 0))
  list(
    phi = unname(phi),
    eta = unname(eta),
    truncated_variance = truncated_variance,
    outcome_exact = outcome_exact,
    treatment_exact = treatment_exact,
    outcome_reference = if (outcome_exact) {
      "exact_sanity_check"
    } else {
      "appendix_f1_approximation"
    },
    treatment_reference = if (treatment_exact) {
      "exact_beta8_zero"
    } else {
      "appendix_f2_approximation"
    }
  )
}

#' Backward-compatible name for analytic bridge-reference parameters
#'
#' `pmtp_oracle_bridge_parameters()` is retained for existing code. The name
#' predates the distinction between the Appendix F.1 outcome approximation and
#' an exact oracle bridge. New code should use
#' [pmtp_analytic_bridge_parameters()].
#'
#' @inheritParams pmtp_analytic_bridge_parameters
#' @return The value returned by [pmtp_analytic_bridge_parameters()].
#' @export
pmtp_oracle_bridge_parameters <- function(spec = pmtp_dgp_spec()) {
  pmtp_analytic_bridge_parameters(spec)
}

pmtp_parametric_h_value <- function(a, l, w, phi, truncated_variance) {
  if (length(phi) != 5L || anyNA(phi) || any(!is.finite(phi))) {
    stop("`phi` must contain five finite values.", call. = FALSE)
  }
  multiplier <- 1 + truncated_variance * phi[4]^2 / 2
  multiplier * stats::plogis(
    phi[1] + phi[2] * a + phi[3] * l + phi[4] * w + phi[5] * a^2
  )
}

pmtp_parametric_g_value <- function(a, l, z, eta, spec) {
  if (length(eta) != 4L || anyNA(eta) || any(!is.finite(eta))) {
    stop("`eta` must contain four finite values.", call. = FALSE)
  }
  v <- pmtp_policy_v(a, spec)
  base <- pmtp_policy_support(a, spec) +
    spec$delta / spec$epsilon *
      as.numeric(a > spec$d - spec$epsilon & a <= spec$d - spec$r * spec$epsilon)
  normalizer <- (stats::pnorm(3) - stats::pnorm(-3)) /
    (stats::pnorm(3 - eta[3] * v) - stats::pnorm(-3 - eta[3] * v))
  log_bridge <- (eta[1] * a + eta[2] * l + eta[3] * z + eta[4] * v) * v
  value <- numeric(length(a))
  inside <- base > 0
  value[inside] <- base[inside] * normalizer[inside] *
    exp(pmin(pmax(log_bridge[inside], -700), 700))
  value
}

#' Evaluate the paper's analytic bridge-reference functions
#'
#' @param a,l,z,w Numeric vectors of equal length.
#' @param spec A specification created by [pmtp_dgp_spec()].
#' @param parameters Optional parameters returned by
#'   [pmtp_analytic_bridge_parameters()].
#'
#' @return A list containing the Appendix F.1 outcome approximation and the
#'   analytic treatment-bridge values. The latter are exact when `beta8 = 0`.
#' @export
pmtp_analytic_bridges <- function(a, l, z, w, spec = pmtp_dgp_spec(),
                                  parameters = NULL) {
  spec <- validate_dgp_spec(spec)
  lengths <- c(length(a), length(l), length(z), length(w))
  if (length(unique(lengths)) != 1L || any(lengths == 0L)) {
    stop("`a`, `l`, `z`, and `w` must have the same positive length.",
         call. = FALSE)
  }
  if (is.null(parameters)) parameters <- pmtp_analytic_bridge_parameters(spec)
  list(
    h = pmtp_parametric_h_value(
      a, l, w, parameters$phi, parameters$truncated_variance
    ),
    g = pmtp_parametric_g_value(a, l, z, parameters$eta, spec)
  )
}

#' Backward-compatible name for analytic bridge-reference functions
#'
#' @inheritParams pmtp_analytic_bridges
#' @return The value returned by [pmtp_analytic_bridges()].
#' @export
pmtp_oracle_bridges <- function(a, l, z, w, spec = pmtp_dgp_spec(),
                                parameters = NULL) {
  pmtp_analytic_bridges(a, l, z, w, spec, parameters)
}

#' Analytic bridge-reference estimators for a paper-DGP sample
#'
#' Calculates the OR, DQW, and DR estimators with the analytic bridge-reference
#' functions. The outcome function is the Appendix F.1 approximation unless
#' its exactness indicator is true; the treatment function is exact when
#' `beta8 = 0`. Versions using the realized outcome and its latent conditional
#' probability are both returned. The latter removes Bernoulli outcome noise
#' and is useful only for simulation diagnostics.
#'
#' @param data Data generated by [simulate_pmtp_dgp()].
#' @param spec A specification created by [pmtp_dgp_spec()]. By default, uses
#'   the specification attached to `data`.
#' @param weights Optional nonnegative analysis weights.
#'
#' @return A list containing estimator values, bridge values, contributions,
#'   and analytic bridge parameters.
#' @export
pmtp_analytic_estimates <- function(data, spec = attr(data, "spec"),
                                    weights = rep(1, nrow(data))) {
  spec <- validate_dgp_spec(spec)
  assert_columns(
    data, c("Y", "A", "L", "Z", "W", "p_y"),
    "analytic-reference data"
  )
  assert_positive(weights, "weights")
  if (length(weights) != nrow(data)) {
    stop("`weights` must have one value per data row.", call. = FALSE)
  }
  q_a <- pmtp_taper_policy(data$A, spec)
  target <- pmtp_policy_target(data$A, spec)
  support <- pmtp_policy_support(data$A, spec)
  parameters <- pmtp_analytic_bridge_parameters(spec)
  observed <- pmtp_analytic_bridges(
    data$A, data$L, data$Z, data$W, spec, parameters
  )
  shifted_h <- numeric(nrow(data))
  shifted_h[target == 1] <- pmtp_parametric_h_value(
    q_a[target == 1], data$L[target == 1], data$W[target == 1],
    parameters$phi, parameters$truncated_variance
  )
  denominator <- sum(weights * target)
  summarize_for <- function(outcome) {
    contributions <- cbind(
      OR = target * shifted_h,
      DQW = observed$g * outcome,
      DR = target * shifted_h + support * observed$g * (outcome - observed$h)
    )
    list(
      estimates = colSums(weights * contributions) / denominator,
      contributions = contributions
    )
  }
  realized <- summarize_for(data$Y)
  conditional_mean <- summarize_for(data$p_y)
  list(
    realized = realized$estimates,
    conditional_mean = conditional_mean$estimates,
    sample_truth = if ("p_policy" %in% names(data)) {
      sum(weights[target == 1] * data$p_policy[target == 1]) / denominator
    } else {
      NA_real_
    },
    bridges = c(observed, list(hq = shifted_h)),
    contributions = list(
      realized = realized$contributions,
      conditional_mean = conditional_mean$contributions
    ),
    parameters = parameters
  )
}

#' Backward-compatible name for analytic bridge-reference estimators
#'
#' @inheritParams pmtp_analytic_estimates
#' @return The value returned by [pmtp_analytic_estimates()].
#' @export
pmtp_oracle_estimates <- function(data, spec = attr(data, "spec"),
                                  weights = rep(1, nrow(data))) {
  pmtp_analytic_estimates(data, spec, weights)
}
