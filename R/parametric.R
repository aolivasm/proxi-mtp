weighted_column_means <- function(x, weights) {
  colSums(weights * x) / sum(weights)
}

match_parametric_outcome_model <- function(model) {
  match.arg(model, c("correct", "omit_quadratic"))
}

match_parametric_treatment_model <- function(model) {
  match.arg(model, c("correct", "constant_v"))
}

parametric_outcome_dimension <- function(model) {
  if (identical(match_parametric_outcome_model(model), "correct")) 5L else 4L
}

pmtp_parametric_h_model_value <- function(
    a, l, w, phi, truncated_variance, outcome_model) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  if (identical(outcome_model, "correct")) {
    return(pmtp_parametric_h_value(a, l, w, phi, truncated_variance))
  }
  if (length(phi) != 4L || anyNA(phi) || any(!is.finite(phi))) {
    stop(
      "`phi` must contain four finite values for `omit_quadratic`.",
      call. = FALSE
    )
  }
  multiplier <- 1 + truncated_variance * phi[4]^2 / 2
  multiplier * stats::plogis(
    phi[1] + phi[2] * a + phi[3] * l + phi[4] * w
  )
}

pmtp_parametric_g_model_value <- function(
    a, l, z, eta, spec, treatment_model) {
  treatment_model <- match_parametric_treatment_model(treatment_model)
  if (identical(treatment_model, "correct")) {
    return(pmtp_parametric_g_value(a, l, z, eta, spec))
  }
  if (length(eta) != 4L || anyNA(eta) || any(!is.finite(eta))) {
    stop("`eta` must contain four finite values.", call. = FALSE)
  }
  v <- pmtp_policy_v(a, spec)
  upper_component <- if (spec$epsilon == 0) {
    numeric(length(a))
  } else {
    spec$delta / spec$epsilon *
      as.numeric(
        a > spec$d - spec$epsilon &
          a <= spec$d - spec$r * spec$epsilon
      )
  }
  base <- pmtp_policy_support(a, spec) + upper_component
  normalizer <- (stats::pnorm(3) - stats::pnorm(-3)) /
    (stats::pnorm(3 - eta[3] * v) - stats::pnorm(-3 - eta[3] * v))
  log_bridge <- (eta[1] * a + eta[2] * l + eta[3] * z + eta[4]) * v
  value <- numeric(length(a))
  inside <- base > 0
  value[inside] <- base[inside] * normalizer[inside] *
    exp(pmin(pmax(log_bridge[inside], -700), 700))
  value
}

parametric_outcome_instruments <- function(data, outcome_model) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  if (identical(outcome_model, "correct")) {
    cbind(1, data$A, data$L, data$Z, data$A^2)
  } else {
    cbind(1, data$A, data$L, data$Z)
  }
}

outcome_bridge_moments <- function(
    phi, data, weights, truncated_variance,
    outcome_model = c("correct", "omit_quadratic")) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  h_value <- pmtp_parametric_h_model_value(
    data$A, data$L, data$W, phi, truncated_variance, outcome_model
  )
  instruments <- parametric_outcome_instruments(data, outcome_model)
  weighted_column_means(instruments * (data$Y - h_value), weights)
}

outcome_bridge_jacobian <- function(
    phi, data, weights, truncated_variance,
    outcome_model = c("correct", "omit_quadratic")) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  design <- if (identical(outcome_model, "correct")) {
    cbind(1, data$A, data$L, data$W, data$A^2)
  } else {
    cbind(1, data$A, data$L, data$W)
  }
  linear_predictor <- drop(design %*% phi)
  probability <- stats::plogis(linear_predictor)
  multiplier <- 1 + truncated_variance * phi[4]^2 / 2
  derivative <- multiplier * probability * (1 - probability) * design
  derivative[, 4L] <- derivative[, 4L] +
    truncated_variance * phi[4] * probability
  instruments <- parametric_outcome_instruments(data, outcome_model)
  -crossprod(instruments, weights * derivative) / sum(weights)
}

treatment_bridge_moments <- function(
    eta, data, weights, spec,
    treatment_model = c("correct", "constant_v")) {
  treatment_model <- match_parametric_treatment_model(treatment_model)
  target <- pmtp_policy_target(data$A, spec)
  q_a <- pmtp_taper_policy(data$A, spec)
  q_for_moment <- ifelse(target == 1, q_a, 0)
  g_value <- pmtp_parametric_g_model_value(
    data$A, data$L, data$Z, eta, spec, treatment_model
  )
  shifted_instruments <- cbind(target, target * q_for_moment, target * data$L,
                               target * data$W)
  observed_instruments <- cbind(1, data$A, data$L, data$W)
  weighted_column_means(
    shifted_instruments - observed_instruments * g_value,
    weights
  )
}

treatment_bridge_jacobian <- function(
    eta, data, weights, spec,
    treatment_model = c("correct", "constant_v")) {
  treatment_model <- match_parametric_treatment_model(treatment_model)
  policy_v <- pmtp_policy_v(data$A, spec)
  bridge <- pmtp_parametric_g_model_value(
    data$A, data$L, data$Z, eta, spec, treatment_model
  )
  normalizer_denominator <-
    stats::pnorm(3 - eta[3] * policy_v) -
    stats::pnorm(-3 - eta[3] * policy_v)
  normalizer_derivative <- policy_v * (
    stats::dnorm(3 - eta[3] * policy_v) -
      stats::dnorm(-3 - eta[3] * policy_v)
  ) / normalizer_denominator
  log_derivative <- cbind(
    data$A * policy_v,
    data$L * policy_v,
    data$Z * policy_v + normalizer_derivative,
    if (identical(treatment_model, "correct")) policy_v^2 else policy_v
  )
  derivative <- bridge * log_derivative
  observed_instruments <- cbind(1, data$A, data$L, data$W)
  -crossprod(observed_instruments, weights * derivative) / sum(weights)
}

parametric_estimating_system <- function(
    parameters, data, weights, population_size, truncated_variance, spec,
    outcome_model = c("correct", "omit_quadratic"),
    treatment_model = c("correct", "constant_v")) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  treatment_model <- match_parametric_treatment_model(treatment_model)
  outcome_dimension <- parametric_outcome_dimension(outcome_model)
  phi <- parameters[seq_len(outcome_dimension)]
  eta <- parameters[outcome_dimension + seq_len(4L)]
  estimates <- parameters[outcome_dimension + 4L + seq_len(3L)]

  q_a <- pmtp_taper_policy(data$A, spec)
  target <- pmtp_policy_target(data$A, spec)
  support <- pmtp_policy_support(data$A, spec)
  q_for_moment <- ifelse(target == 1, q_a, 0)

  h0 <- pmtp_parametric_h_model_value(
    data$A, data$L, data$W, phi, truncated_variance, outcome_model
  )
  hq <- numeric(nrow(data))
  hq[target == 1] <- pmtp_parametric_h_model_value(
    q_a[target == 1], data$L[target == 1], data$W[target == 1],
    phi, truncated_variance, outcome_model
  )
  g0 <- pmtp_parametric_g_model_value(
    data$A, data$L, data$Z, eta, spec, treatment_model
  )

  outcome_instruments <- parametric_outcome_instruments(data, outcome_model)
  outcome_moments <- outcome_instruments * (data$Y - h0)
  outcome_design <- if (identical(outcome_model, "correct")) {
    cbind(1, data$A, data$L, data$W, data$A^2)
  } else {
    cbind(1, data$A, data$L, data$W)
  }
  outcome_probability <- stats::plogis(drop(outcome_design %*% phi))
  outcome_multiplier <- 1 + truncated_variance * phi[4]^2 / 2
  outcome_derivative <-
    outcome_multiplier * outcome_probability * (1 - outcome_probability) *
    outcome_design
  outcome_derivative[, 4L] <- outcome_derivative[, 4L] +
    truncated_variance * phi[4] * outcome_probability

  shifted_instruments <- cbind(
    target, target * q_for_moment, target * data$L, target * data$W
  )
  observed_instruments <- cbind(1, data$A, data$L, data$W)
  treatment_moments <-
    shifted_instruments - observed_instruments * g0
  policy_v <- pmtp_policy_v(data$A, spec)
  normalizer_denominator <-
    stats::pnorm(3 - eta[3] * policy_v) -
    stats::pnorm(-3 - eta[3] * policy_v)
  normalizer_derivative <- policy_v * (
    stats::dnorm(3 - eta[3] * policy_v) -
      stats::dnorm(-3 - eta[3] * policy_v)
  ) / normalizer_denominator
  treatment_log_derivative <- cbind(
    data$A * policy_v,
    data$L * policy_v,
    data$Z * policy_v + normalizer_derivative,
    if (identical(treatment_model, "correct")) policy_v^2 else policy_v
  )
  treatment_derivative <- g0 * treatment_log_derivative

  weight_total <- sum(weights)
  outcome_mean <- weighted_column_means(outcome_moments, weights)
  treatment_mean <- weighted_column_means(treatment_moments, weights)
  outcome_jacobian <-
    -crossprod(outcome_instruments, weights * outcome_derivative) /
    weight_total
  treatment_jacobian <-
    -crossprod(observed_instruments, weights * treatment_derivative) /
    weight_total
  outcome_stationarity <- drop(crossprod(outcome_jacobian, outcome_mean))
  treatment_stationarity <-
    drop(crossprod(treatment_jacobian, treatment_mean))

  # These are the contamination derivatives of J(P)' m(P), centered at the
  # empirical distribution. They reduce to J' psi when the bridge moments have
  # an exact root, but remain valid for the minimum-distance fallback.
  outcome_scores <-
    -outcome_derivative * drop(outcome_instruments %*% outcome_mean) +
    outcome_moments %*% outcome_jacobian
  outcome_scores <- sweep(
    outcome_scores, 2L, 2 * outcome_stationarity, "-"
  )
  treatment_scores <-
    -treatment_derivative * drop(observed_instruments %*% treatment_mean) +
    treatment_moments %*% treatment_jacobian
  treatment_scores <- sweep(
    treatment_scores, 2L, 2 * treatment_stationarity, "-"
  )

  contributions <- cbind(
    OR = target * hq,
    DQW = g0 * data$Y,
    DR = target * hq + support * g0 * (data$Y - h0)
  )
  estimate_equations <- contributions - outer(target, estimates)

  scores <- cbind(
    outcome_scores, treatment_scores, estimate_equations
  )
  colnames(scores) <- c(
    paste0("h_minimum_distance_", seq_len(outcome_dimension)),
    paste0("g_minimum_distance_", seq_len(4L)),
    paste0("target_", colnames(contributions))
  )
  system <- c(
    weight_total / population_size * outcome_stationarity,
    weight_total / population_size * treatment_stationarity,
    colSums(weights * estimate_equations) / population_size
  )
  names(system) <- colnames(scores)
  list(scores = scores, system = system)
}

finite_difference_jacobian <- function(
    function_value, parameters,
    relative_step = .Machine$double.eps^(1 / 3)) {
  baseline <- function_value(parameters)
  if (!is.numeric(baseline) || any(!is.finite(baseline))) {
    stop("The parametric estimating equations are not finite.", call. = FALSE)
  }
  jacobian <- matrix(
    NA_real_, nrow = length(baseline), ncol = length(parameters)
  )
  for (column in seq_along(parameters)) {
    step <- relative_step * max(1, abs(parameters[column]))
    upper <- lower <- parameters
    upper[column] <- upper[column] + step
    lower[column] <- lower[column] - step
    upper_value <- function_value(upper)
    lower_value <- function_value(lower)
    derivative <- (upper_value - lower_value) / (2 * step)
    if (any(!is.finite(derivative))) {
      stop(
        "A finite-difference derivative for parametric inference is not finite.",
        call. = FALSE
      )
    }
    jacobian[, column] <- derivative
  }
  dimnames(jacobian) <- list(names(baseline), names(parameters))
  jacobian
}

generalized_inverse <- function(matrix) {
  decomposition <- svd(matrix)
  if (!length(decomposition$d) || !all(is.finite(decomposition$d))) {
    stop("The parametric estimating-equation Jacobian is invalid.",
         call. = FALSE)
  }
  largest <- max(decomposition$d)
  cutoff <- max(dim(matrix)) * .Machine$double.eps^0.75 * largest
  keep <- decomposition$d > cutoff
  if (!any(keep)) {
    stop("The parametric estimating-equation Jacobian has rank zero.",
         call. = FALSE)
  }
  inverse <- decomposition$v[, keep, drop = FALSE] %*%
    (t(decomposition$u[, keep, drop = FALSE]) /
       decomposition$d[keep])
  list(
    inverse = inverse,
    rank = sum(keep),
    condition_number = largest / min(decomposition$d[keep]),
    singular_values = decomposition$d,
    cutoff = cutoff
  )
}

parametric_inference <- function(
    phi, eta, estimates, data, weights, population_size,
    truncated_variance, spec,
    outcome_model = c("correct", "omit_quadratic"),
    treatment_model = c("correct", "constant_v")) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  treatment_model <- match_parametric_treatment_model(treatment_model)
  parameters <- c(
    stats::setNames(phi, paste0("phi", seq_along(phi) - 1L)),
    stats::setNames(eta, paste0("eta", 0:3)),
    estimates
  )
  components <- parametric_estimating_system(
    parameters, data, weights, population_size, truncated_variance, spec,
    outcome_model, treatment_model
  )
  system_function <- function(candidate) {
    values <- parametric_estimating_system(
      candidate, data, weights, population_size, truncated_variance, spec,
      outcome_model, treatment_model
    )
    values$system
  }
  jacobian <- finite_difference_jacobian(system_function, parameters)
  inverse <- generalized_inverse(jacobian)
  if (inverse$rank < length(parameters)) {
    warning(
      "The parametric estimating-equation Jacobian is rank deficient; ",
      "inference uses its singular-value generalized inverse.",
      call. = FALSE
    )
  }

  parameter_influence <-
    -components$scores %*% t(inverse$inverse)
  colnames(parameter_influence) <- names(parameters)
  influence <- parameter_influence[, names(estimates), drop = FALSE]
  weighted_influence <- weights * influence
  covariance <- crossprod(weighted_influence) / population_size^2
  dimnames(covariance) <- list(names(estimates), names(estimates))
  standard_error <- sqrt(pmax(diag(covariance), 0))
  names(standard_error) <- names(estimates)

  list(
    standard_error = standard_error,
    covariance = covariance,
    asymptotic_variance = stats::setNames(
      diag(covariance) * population_size, names(estimates)
    ),
    influence_function = influence,
    jacobian = jacobian,
    jacobian_rank = inverse$rank,
    jacobian_condition_number = inverse$condition_number,
    jacobian_singular_values = inverse$singular_values,
    equation_residual = system_function(parameters)
  )
}

unique_starts <- function(starts) {
  keys <- vapply(starts, function(x) paste(signif(x, 12), collapse = "|"), character(1))
  starts[!duplicated(keys)]
}

solve_bridge_moments <- function(
    moment_function, starts, max_iterations, jacobian_function = NULL) {
  starts <- unique_starts(starts)
  candidates <- vector("list", 0L)
  root_tolerance <- 1e-8
  stationarity_tolerance <- 1e-6
  add_candidate <- function(coefficients, method, convergence) {
    moments <- tryCatch(moment_function(coefficients), error = function(e) rep(Inf, length(coefficients)))
    if (length(moments) == length(coefficients) && all(is.finite(moments))) {
      stationarity <- if (is.null(jacobian_function)) {
        rep(NA_real_, length(coefficients))
      } else {
        jacobian <- tryCatch(
          jacobian_function(coefficients),
          error = function(e) matrix(NA_real_, length(coefficients), length(coefficients))
        )
        if (all(is.finite(jacobian))) {
          drop(crossprod(jacobian, moments))
        } else {
          rep(NA_real_, length(coefficients))
        }
      }
      candidates[[length(candidates) + 1L]] <<- list(
        coefficients = coefficients,
        moments = moments,
        residual_norm = sqrt(sum(moments^2)),
        stationarity_norm = sqrt(sum(stationarity^2)),
        method = method,
        convergence = convergence
      )
    }
  }

  try_root <- function(start, method, label) {
    root <- tryCatch(
      nleqslv::nleqslv(
        start,
        moment_function,
        jac = jacobian_function,
        method = method,
        global = "dbldog",
        control = list(
          ftol = 1e-10,
          xtol = 1e-10,
          maxit = as.integer(max_iterations),
          allowSingular = TRUE
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(root)) {
      add_candidate(
        root$x,
        paste0("nleqslv-", tolower(method), "-", label),
        root$termcd
      )
    }
    invisible(root)
  }

  for (start_index in seq_along(starts)) {
    start <- starts[[start_index]]
    try_root(start, "Newton", paste0("start-", start_index))
    if (length(candidates) &&
        candidates[[length(candidates)]]$residual_norm <= root_tolerance) {
      best <- candidates[[length(candidates)]]
      best$converged <- TRUE
      best$solution_type <- "root"
      return(best)
    }
    try_root(start, "Broyden", paste0("start-", start_index))
    if (length(candidates) &&
        candidates[[length(candidates)]]$residual_norm <= root_tolerance) {
      best <- candidates[[length(candidates)]]
      best$converged <- TRUE
      best$solution_type <- "root"
      return(best)
    }
    add_candidate(start, paste0("initial-start-", start_index), NA_integer_)
  }

  if (!length(candidates)) {
    stop("No finite bridge estimating-equation candidate was obtained.", call. = FALSE)
  }
  candidate_order <- order(vapply(
    candidates, `[[`, numeric(1), "residual_norm"
  ))
  optimization_starts <- unique_starts(lapply(
    candidate_order[seq_len(min(4L, length(candidate_order)))],
    function(index) candidates[[index]]$coefficients
  ))
  for (start_index in seq_along(optimization_starts)) {
    optimized <- tryCatch(
      stats::optim(
        optimization_starts[[start_index]],
        function(x) {
          moments <- moment_function(x)
          if (any(!is.finite(moments))) return(.Machine$double.xmax / 100)
          sum(moments^2)
        },
        gr = if (is.null(jacobian_function)) {
          NULL
        } else {
          function(x) {
            moments <- moment_function(x)
            jacobian <- jacobian_function(x)
            2 * drop(crossprod(jacobian, moments))
          }
        },
        method = "BFGS",
        control = list(
          maxit = as.integer(max(2000L, 10L * max_iterations)),
          reltol = 1e-12
        )
      ),
      error = function(e) NULL
    )
    if (is.null(optimized)) next
    add_candidate(
      optimized$par, paste0("optim-BFGS-start-", start_index),
      optimized$convergence
    )
    try_root(
      optimized$par, "Newton",
      paste0("optim-refinement-", start_index)
    )
    if (length(candidates) &&
        candidates[[length(candidates)]]$residual_norm <= root_tolerance) {
      best <- candidates[[length(candidates)]]
      best$converged <- TRUE
      best$solution_type <- "root"
      return(best)
    }
  }
  root_candidates <- which(vapply(
    candidates,
    function(candidate) is.finite(candidate$residual_norm) &&
      candidate$residual_norm <= 1e-6,
    logical(1)
  ))
  stationary_candidates <- which(vapply(
    candidates,
    function(candidate) is.finite(candidate$stationarity_norm) &&
      candidate$stationarity_norm <= stationarity_tolerance,
    logical(1)
  ))
  if (length(root_candidates)) {
    best_index <- root_candidates[which.min(vapply(
      candidates[root_candidates], `[[`, numeric(1), "residual_norm"
    ))]
    solution_type <- "root"
  } else if (length(stationary_candidates)) {
    best_index <- stationary_candidates[which.min(vapply(
      candidates[stationary_candidates], `[[`, numeric(1), "residual_norm"
    ))]
    solution_type <- "minimum_distance"
  } else {
    best_index <- which.min(vapply(
      candidates, `[[`, numeric(1), "stationarity_norm"
    ))
    solution_type <- "failed"
  }
  best <- candidates[[best_index]]
  best$converged <- !identical(solution_type, "failed")
  best$solution_type <- solution_type
  best
}

initial_outcome_coefficients <- function(
    data, outcome_model = c("correct", "omit_quadratic")) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  design <- if (identical(outcome_model, "correct")) {
    cbind(1, data$A, data$L, data$W, data$A^2)
  } else {
    cbind(1, data$A, data$L, data$W)
  }
  fit <- suppressWarnings(tryCatch(
    stats::glm.fit(design, data$Y, family = stats::binomial()),
    error = function(e) NULL
  ))
  if (is.null(fit)) return(rep(0, ncol(design)))
  coefficients <- fit$coefficients
  coefficients[!is.finite(coefficients)] <- 0
  unname(coefficients)
}

prepare_parametric_inputs <- function(
    data, spec, weights, population_size, max_iterations) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  spec <- validate_dgp_spec(spec)
  assert_columns(data, c("Y", "A", "L", "Z", "W"), "parametric data")
  if (!all(vapply(
    data[c("Y", "A", "L", "Z", "W")], is.numeric, logical(1)
  )) || anyNA(data[c("Y", "A", "L", "Z", "W")])) {
    stop(
      "`Y`, `A`, `L`, `Z`, and `W` must be complete numeric columns.",
      call. = FALSE
    )
  }
  analysis_weights <- resolve_vector_argument(
    data, weights, "weights", default = 1
  )
  assert_positive(analysis_weights, "weights")
  if (is.null(population_size)) population_size <- sum(analysis_weights)
  assert_positive(population_size, "population_size")
  if (length(population_size) != 1L) {
    stop("`population_size` must be a single positive value.", call. = FALSE)
  }
  if (length(max_iterations) != 1L || is.na(max_iterations) ||
      max_iterations < 1 || max_iterations != as.integer(max_iterations)) {
    stop("`max_iterations` must be a positive integer.", call. = FALSE)
  }
  list(
    data = data,
    spec = spec,
    weights = analysis_weights,
    population_size = population_size,
    weighted = !is.null(weights),
    max_iterations = as.integer(max_iterations)
  )
}

solve_parametric_outcome_bridge <- function(
    data, weights, truncated_variance, outcome_model, start_h,
    max_iterations) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  dimension <- parametric_outcome_dimension(outcome_model)
  if (!is.null(start_h) &&
      (length(start_h) != dimension || any(!is.finite(start_h)))) {
    stop(
      "`start_h` must contain ", dimension, " finite values for `",
      outcome_model, "`.",
      call. = FALSE
    )
  }
  starts <- list(
    initial_outcome_coefficients(data, outcome_model),
    rep(0, dimension)
  )
  if (!is.null(start_h)) {
    start_h <- as.numeric(start_h)
    starts <- c(list(start_h), starts)
    if (abs(start_h[4L]) >= 1.5) {
      perturbations <- if (dimension == 5L) {
        rbind(
          c(0, 0, -1, -0.5, 0),
          c(-1, 0, -1, -0.5, 0),
          c(0, -0.2, -1, -0.5, 0.2)
        )
      } else {
        rbind(
          c(0, 0, -1, -0.5),
          c(-1, 0, -1, -0.5),
          c(0, -0.2, -1, -0.5)
        )
      }
      perturbed <- lapply(
        seq_len(nrow(perturbations)),
        function(index) start_h + perturbations[index, ]
      )
      starts <- c(list(start_h), perturbed, starts[-1L])
    }
  }
  solve_bridge_moments(
    function(phi) {
      outcome_bridge_moments(
        phi, data, weights, truncated_variance, outcome_model
      )
    },
    starts,
    max_iterations,
    function(phi) {
      outcome_bridge_jacobian(
        phi, data, weights, truncated_variance, outcome_model
      )
    }
  )
}

solve_parametric_treatment_bridge <- function(
    data, weights, spec, treatment_model, start_g, max_iterations) {
  treatment_model <- match_parametric_treatment_model(treatment_model)
  if (!is.null(start_g) &&
      (length(start_g) != 4L || any(!is.finite(start_g)))) {
    stop("`start_g` must contain four finite values.", call. = FALSE)
  }
  starts <- list(rep(0, 4L), c(0.1, 0, -0.1, -0.1))
  if (!is.null(start_g)) starts <- c(list(as.numeric(start_g)), starts)
  solve_bridge_moments(
    function(eta) {
      treatment_bridge_moments(
        eta, data, weights, spec, treatment_model
      )
    },
    starts,
    max_iterations,
    function(eta) {
      treatment_bridge_jacobian(
        eta, data, weights, spec, treatment_model
      )
    }
  )
}

parametric_estimate_table <- function(
    estimates, standard_error, conf_level = 0.95) {
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  data.frame(
    estimator = names(estimates),
    estimate = unname(estimates),
    std_error = unname(standard_error),
    conf_low = unname(estimates - critical * standard_error),
    conf_high = unname(estimates + critical * standard_error),
    row.names = NULL
  )
}

assemble_parametric_fit <- function(
    data, spec, analysis_weights, population_size, weighted,
    bridge_parameters, h_solution, g_solution, outcome_model,
    treatment_model) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  treatment_model <- match_parametric_treatment_model(treatment_model)
  q_a <- pmtp_taper_policy(data$A, spec)
  target <- pmtp_policy_target(data$A, spec)
  support <- pmtp_policy_support(data$A, spec)
  h0 <- pmtp_parametric_h_model_value(
    data$A, data$L, data$W, h_solution$coefficients,
    bridge_parameters$truncated_variance, outcome_model
  )
  hq <- numeric(nrow(data))
  hq[target == 1] <- pmtp_parametric_h_model_value(
    q_a[target == 1], data$L[target == 1], data$W[target == 1],
    h_solution$coefficients, bridge_parameters$truncated_variance,
    outcome_model
  )
  g0 <- pmtp_parametric_g_model_value(
    data$A, data$L, data$Z, g_solution$coefficients, spec,
    treatment_model
  )
  denominator <- sum(analysis_weights * target)
  contributions <- cbind(
    OR = target * hq,
    DQW = g0 * data$Y,
    DR = target * hq + support * g0 * (data$Y - h0)
  )
  estimates <- colSums(analysis_weights * contributions) / denominator
  inference <- parametric_inference(
    h_solution$coefficients, g_solution$coefficients, estimates,
    data, analysis_weights, population_size,
    bridge_parameters$truncated_variance, spec,
    outcome_model, treatment_model
  )
  estimate_table <- parametric_estimate_table(
    estimates, inference$standard_error
  )

  structure(list(
    estimates = estimates,
    standard_error = inference$standard_error,
    estimate_table = estimate_table,
    covariance = inference$covariance,
    asymptotic_variance = inference$asymptotic_variance,
    influence_function = inference$influence_function,
    inference = inference[c(
      "jacobian", "jacobian_rank", "jacobian_condition_number",
      "jacobian_singular_values", "equation_residual"
    )],
    coefficients = list(
      phi = h_solution$coefficients,
      eta = g_solution$coefficients
    ),
    moments = list(h = h_solution$moments, g = g_solution$moments),
    residual_norm = c(
      h = h_solution$residual_norm,
      g = g_solution$residual_norm
    ),
    stationarity_norm = c(
      h = h_solution$stationarity_norm,
      g = g_solution$stationarity_norm
    ),
    solution_type = c(
      h = h_solution$solution_type,
      g = g_solution$solution_type
    ),
    converged = c(h = h_solution$converged, g = g_solution$converged),
    solver = list(
      h = h_solution[c(
        "method", "convergence", "solution_type", "stationarity_norm"
      )],
      g = g_solution[c(
        "method", "convergence", "solution_type", "stationarity_norm"
      )]
    ),
    nuisance = list(h0 = h0, hq = hq, g0 = g0),
    contributions = contributions,
    weights = analysis_weights,
    population_size = population_size,
    target = target,
    policy_support = support,
    target_probability = denominator / population_size,
    weighted = weighted,
    models = list(
      outcome = outcome_model,
      treatment = treatment_model
    ),
    spec = spec
  ), class = "pmtp_parametric_fit")
}

paper_parametric_start_h <- function(spec, fallback) {
  calibrated <- rbind(
    `-2` = c(
      -1.19960935243, -1.48333752394, 0.237030425319,
      0.471832194950, -0.742492117225
    ),
    `-1` = c(
      -1.64333565913, -1.47560415970, 0.00425279455721,
      0.879024663645, -0.736654844877
    ),
    `-0.5` = c(
      -2.99444625370, -1.51880248120, -0.423736151791,
      1.696944153087, -0.755909320377
    ),
    `-0.25` = c(
      -7.4326, -1.8163, -1.3801, 3.6970, -0.86994
    )
  )
  proxy_coefficient <- unname(spec$beta[["beta5"]])
  index <- which(
    abs(as.numeric(rownames(calibrated)) - proxy_coefficient) < 1e-12
  )
  if (length(index) == 1L) unname(calibrated[index, ]) else fallback
}

paper_parametric_start_g_misspecified <- function(spec, fallback) {
  calibrated <- rbind(
    `2` = c(0.41516779, -0.077018686, -0.20322019, -0.099837796),
    `1` = c(0.42167649, -0.036646768, -0.40508250, -0.158808250),
    `0.5` = c(0.44634155, 0.042526931, -0.80144240, -0.386970780),
    `0.25` = c(0.52968742, 0.188863270, -1.53891300, -1.211547000)
  )
  proxy_coefficient <- unname(spec$beta[["beta3"]])
  index <- which(
    abs(as.numeric(rownames(calibrated)) - proxy_coefficient) < 1e-12
  )
  if (length(index) == 1L) unname(calibrated[index, ]) else fallback
}

#' Fit the paper's proximal parametric bridge estimators
#'
#' Fits the estimating equations in Supplement C.3 and returns the proximal
#' outcome-regression, density-quotient-weighted, and doubly robust point
#' estimates. An exact finite-sample root is used when available. Otherwise,
#' the bridge fit is the minimum-distance minimizer of the squared empirical
#' moments. Standard errors use the corresponding stacked influence function,
#' jointly accounting for bridge estimation and the target population
#' denominator. With inverse two-phase sampling weights, the covariance uses
#' squared inverse weights and the supplied phase-one population size.
#'
#' @param data A data frame containing scalar columns `Y`, `A`, `L`, `Z`, and
#'   `W`.
#' @param spec A specification created by [pmtp_dgp_spec()]. Only its known
#'   policy parameters are used by the estimating equations.
#' @param weights Optional positive observation weights or the name of a column
#'   in `data` containing them.
#' @param population_size Full phase-one population size. By default, uses the
#'   sum of the analysis weights. Supplying the known population size is
#'   recommended with inverse two-phase sampling weights.
#' @param outcome_model Outcome-bridge specification. `"correct"` includes the
#'   quadratic treatment term; `"omit_quadratic"` is the misspecified model
#'   used in Supplement C.3.
#' @param treatment_model Treatment-bridge specification. `"correct"` includes
#'   `eta_3 * V(a)` inside the exponential's linear predictor;
#'   `"constant_v"` replaces that term with the constant `eta_3`, as in the
#'   supplement's misspecified model.
#' @param start_h,start_g Optional starting coefficients for the outcome and
#'   treatment bridge estimating equations.
#' @param max_iterations Maximum iterations for each numerical solver.
#'
#' @return An object of class `pmtp_parametric_fit` containing estimates,
#'   influence-function standard errors and covariance, coefficients, moment
#'   residuals, convergence information, and nuisance values.
#' @export
pmtp_parametric <- function(data, spec = pmtp_dgp_spec(), weights = NULL,
                            population_size = NULL,
                            start_h = NULL, start_g = NULL,
                            max_iterations = 500L,
                            outcome_model = c("correct", "omit_quadratic"),
                            treatment_model = c("correct", "constant_v")) {
  outcome_model <- match_parametric_outcome_model(outcome_model)
  treatment_model <- match_parametric_treatment_model(treatment_model)
  inputs <- prepare_parametric_inputs(
    data, spec, weights, population_size, max_iterations
  )
  bridge_parameters <- pmtp_analytic_bridge_parameters(inputs$spec)
  h_solution <- solve_parametric_outcome_bridge(
    inputs$data, inputs$weights, bridge_parameters$truncated_variance,
    outcome_model, start_h, inputs$max_iterations
  )
  g_solution <- solve_parametric_treatment_bridge(
    inputs$data, inputs$weights, inputs$spec, treatment_model, start_g,
    inputs$max_iterations
  )
  assemble_parametric_fit(
    inputs$data, inputs$spec, inputs$weights, inputs$population_size,
    inputs$weighted, bridge_parameters, h_solution, g_solution,
    outcome_model, treatment_model
  )
}

#' Fit the supplement's eight proximal parametric estimators
#'
#' Fits the correctly specified and deliberately misspecified outcome and
#' treatment bridge models from Supplement C.3. Each of the four bridge
#' functions is fitted once. The shared fits are then combined into two OR,
#' two DQW, and four DR estimators with stacked influence-function inference.
#'
#' @inheritParams pmtp_parametric
#' @param start_h_correct,start_h_misspecified Optional starting values for the
#'   correct and misspecified outcome bridge models.
#' @param start_g_correct,start_g_misspecified Optional starting values for the
#'   correct and misspecified treatment bridge models.
#'
#' @return An object of class `pmtp_parametric_suite` containing the eight
#'   estimates, their joint covariance matrix and influence functions, the
#'   four component fits, and the shared bridge solutions. The `solution_type`
#'   component identifies exact roots and minimum-distance solutions.
#' @export
pmtp_parametric_suite <- function(
    data, spec = pmtp_dgp_spec(), weights = NULL, population_size = NULL,
    start_h_correct = NULL, start_h_misspecified = NULL,
    start_g_correct = NULL, start_g_misspecified = NULL,
    max_iterations = 500L) {
  inputs <- prepare_parametric_inputs(
    data, spec, weights, population_size, max_iterations
  )
  bridge_parameters <- pmtp_analytic_bridge_parameters(inputs$spec)
  if (is.null(start_h_correct)) {
    start_h_correct <- paper_parametric_start_h(
      inputs$spec, bridge_parameters$phi
    )
  }
  if (is.null(start_h_misspecified)) {
    start_h_misspecified <- start_h_correct[seq_len(4L)]
  }
  if (is.null(start_g_correct)) {
    start_g_correct <- bridge_parameters$eta
  }
  if (is.null(start_g_misspecified)) {
    start_g_misspecified <- paper_parametric_start_g_misspecified(
      inputs$spec, bridge_parameters$eta
    )
  }

  h_solutions <- list(
    correct = solve_parametric_outcome_bridge(
      inputs$data, inputs$weights, bridge_parameters$truncated_variance,
      "correct", start_h_correct, inputs$max_iterations
    ),
    misspecified = solve_parametric_outcome_bridge(
      inputs$data, inputs$weights, bridge_parameters$truncated_variance,
      "omit_quadratic", start_h_misspecified, inputs$max_iterations
    )
  )
  g_solutions <- list(
    correct = solve_parametric_treatment_bridge(
      inputs$data, inputs$weights, inputs$spec, "correct",
      start_g_correct, inputs$max_iterations
    ),
    misspecified = solve_parametric_treatment_bridge(
      inputs$data, inputs$weights, inputs$spec, "constant_v",
      start_g_misspecified, inputs$max_iterations
    )
  )

  assemble <- function(h_name, g_name) {
    assemble_parametric_fit(
      inputs$data, inputs$spec, inputs$weights, inputs$population_size,
      inputs$weighted, bridge_parameters,
      h_solutions[[h_name]], g_solutions[[g_name]],
      if (identical(h_name, "correct")) "correct" else "omit_quadratic",
      if (identical(g_name, "correct")) "correct" else "constant_v"
    )
  }
  fits <- list(
    h_correct_g_correct = assemble("correct", "correct"),
    h_correct_g_misspecified = assemble("correct", "misspecified"),
    h_misspecified_g_correct = assemble("misspecified", "correct"),
    h_misspecified_g_misspecified =
      assemble("misspecified", "misspecified")
  )

  estimates <- c(
    OR_h_correct =
      fits$h_correct_g_correct$estimates[["OR"]],
    OR_h_misspecified =
      fits$h_misspecified_g_correct$estimates[["OR"]],
    DQW_g_correct =
      fits$h_correct_g_correct$estimates[["DQW"]],
    DQW_g_misspecified =
      fits$h_correct_g_misspecified$estimates[["DQW"]],
    DR_h_correct_g_correct =
      fits$h_correct_g_correct$estimates[["DR"]],
    DR_h_correct_g_misspecified =
      fits$h_correct_g_misspecified$estimates[["DR"]],
    DR_h_misspecified_g_correct =
      fits$h_misspecified_g_correct$estimates[["DR"]],
    DR_h_misspecified_g_misspecified =
      fits$h_misspecified_g_misspecified$estimates[["DR"]]
  )
  influence_function <- cbind(
    OR_h_correct =
      fits$h_correct_g_correct$influence_function[, "OR"],
    OR_h_misspecified =
      fits$h_misspecified_g_correct$influence_function[, "OR"],
    DQW_g_correct =
      fits$h_correct_g_correct$influence_function[, "DQW"],
    DQW_g_misspecified =
      fits$h_correct_g_misspecified$influence_function[, "DQW"],
    DR_h_correct_g_correct =
      fits$h_correct_g_correct$influence_function[, "DR"],
    DR_h_correct_g_misspecified =
      fits$h_correct_g_misspecified$influence_function[, "DR"],
    DR_h_misspecified_g_correct =
      fits$h_misspecified_g_correct$influence_function[, "DR"],
    DR_h_misspecified_g_misspecified =
      fits$h_misspecified_g_misspecified$influence_function[, "DR"]
  )
  weighted_influence <- inputs$weights * influence_function
  covariance <- crossprod(weighted_influence) / inputs$population_size^2
  standard_error <- sqrt(pmax(diag(covariance), 0))
  names(standard_error) <- names(estimates)
  dimnames(covariance) <- list(names(estimates), names(estimates))

  structure(list(
    estimates = estimates,
    standard_error = standard_error,
    estimate_table = parametric_estimate_table(
      estimates, standard_error
    ),
    covariance = covariance,
    asymptotic_variance = stats::setNames(
      diag(covariance) * inputs$population_size, names(estimates)
    ),
    influence_function = influence_function,
    fits = fits,
    bridge_solutions = list(h = h_solutions, g = g_solutions),
    converged = c(
      h_correct = h_solutions$correct$converged,
      h_misspecified = h_solutions$misspecified$converged,
      g_correct = g_solutions$correct$converged,
      g_misspecified = g_solutions$misspecified$converged
    ),
    residual_norm = c(
      h_correct = h_solutions$correct$residual_norm,
      h_misspecified = h_solutions$misspecified$residual_norm,
      g_correct = g_solutions$correct$residual_norm,
      g_misspecified = g_solutions$misspecified$residual_norm
    ),
    stationarity_norm = c(
      h_correct = h_solutions$correct$stationarity_norm,
      h_misspecified = h_solutions$misspecified$stationarity_norm,
      g_correct = g_solutions$correct$stationarity_norm,
      g_misspecified = g_solutions$misspecified$stationarity_norm
    ),
    solution_type = c(
      h_correct = h_solutions$correct$solution_type,
      h_misspecified = h_solutions$misspecified$solution_type,
      g_correct = g_solutions$correct$solution_type,
      g_misspecified = g_solutions$misspecified$solution_type
    ),
    weights = inputs$weights,
    population_size = inputs$population_size,
    weighted = inputs$weighted,
    spec = inputs$spec
  ), class = "pmtp_parametric_suite")
}

#' @export
print.pmtp_parametric_fit <- function(x, ...) {
  cat(
    "Proximal parametric MTP estimates (h: ", x$models$outcome,
    "; g: ", x$models$treatment, ")\n",
    sep = ""
  )
  print(x$estimate_table, row.names = FALSE)
  cat("\nEstimating-equation residual norms\n")
  print(x$residual_norm)
  cat("\nBridge solution types\n")
  print(x$solution_type)
  invisible(x)
}

#' @export
coef.pmtp_parametric_fit <- function(object, ...) {
  object$estimates
}

#' @export
vcov.pmtp_parametric_fit <- function(object, ...) {
  object$covariance
}

#' @export
summary.pmtp_parametric_fit <- function(object, conf_level = 0.95, ...) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must lie strictly between zero and one.", call. = FALSE)
  }
  output <- object$estimate_table
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  output$conf_low <- output$estimate - critical * output$std_error
  output$conf_high <- output$estimate + critical * output$std_error
  attr(output, "conf_level") <- conf_level
  class(output) <- c("summary_pmtp_parametric_fit", class(output))
  output
}

#' @export
print.summary_pmtp_parametric_fit <- function(x, ...) {
  level <- 100 * attr(x, "conf_level")
  cat(
    "Proximal parametric MTP estimates with ", level,
    "% confidence intervals\n\n", sep = ""
  )
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' @export
print.pmtp_parametric_suite <- function(x, ...) {
  cat("Supplement C.3 proximal parametric estimator suite\n")
  print(x$estimate_table, row.names = FALSE)
  cat("\nEstimating-equation residual norms\n")
  print(x$residual_norm)
  cat("\nBridge solution types\n")
  print(x$solution_type)
  invisible(x)
}

#' @export
coef.pmtp_parametric_suite <- function(object, ...) {
  object$estimates
}

#' @export
vcov.pmtp_parametric_suite <- function(object, ...) {
  object$covariance
}

#' @export
summary.pmtp_parametric_suite <- function(
    object, conf_level = 0.95, ...) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must lie strictly between zero and one.", call. = FALSE)
  }
  output <- parametric_estimate_table(
    object$estimates, object$standard_error, conf_level
  )
  attr(output, "conf_level") <- conf_level
  class(output) <- c("summary_pmtp_parametric_suite", class(output))
  output
}

#' @export
print.summary_pmtp_parametric_suite <- function(x, ...) {
  level <- 100 * attr(x, "conf_level")
  cat(
    "Supplement C.3 proximal parametric estimates with ", level,
    "% confidence intervals\n\n", sep = ""
  )
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}
