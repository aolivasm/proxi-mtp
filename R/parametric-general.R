quote_formula_name <- function(name) {
  paste0("`", gsub("`", "\\\\`", name, fixed = TRUE), "`")
}

general_parametric_default_formulas <- function(
    treatment, covariates, negative_control_treatment,
    negative_control_outcome) {
  treatment_term <- quote_formula_name(treatment)
  quadratic_term <- paste0("I(", treatment_term, "^2)")
  covariate_terms <- vapply(
    covariates, quote_formula_name, character(1L)
  )
  z_terms <- vapply(
    negative_control_treatment, quote_formula_name, character(1L)
  )
  w_terms <- vapply(
    negative_control_outcome, quote_formula_name, character(1L)
  )
  list(
    outcome_bridge = stats::reformulate(c(
      treatment_term, quadratic_term, covariate_terms, w_terms
    )),
    outcome_instruments = stats::reformulate(c(
      treatment_term, quadratic_term, covariate_terms, z_terms
    )),
    treatment_bridge = stats::reformulate(c(
      treatment_term, covariate_terms, z_terms
    )),
    treatment_instruments = stats::reformulate(c(
      treatment_term, covariate_terms, w_terms
    ))
  )
}

normalize_general_parametric_formula <- function(formula, default, name) {
  if (is.null(formula)) formula <- default
  if (is.character(formula) && length(formula) == 1L) {
    formula <- stats::as.formula(formula, env = parent.frame())
  }
  if (!inherits(formula, "formula")) {
    stop("`", name, "` must be a formula or one formula string.",
         call. = FALSE)
  }
  terms <- stats::terms(formula)
  if (attr(terms, "response") > 0L) {
    terms <- stats::delete.response(terms)
  }
  terms
}

weighted_matrix_center_scale <- function(matrix, weights, standardize) {
  number_columns <- ncol(matrix)
  center <- rep(0, number_columns)
  scale <- rep(1, number_columns)
  if (!standardize || !number_columns) {
    return(list(center = center, scale = scale))
  }
  weight_total <- sum(weights)
  means <- colSums(weights * matrix) / weight_total
  variances <- colSums(weights * (matrix - rep(
    means, each = nrow(matrix)
  ))^2) / weight_total
  standard_deviations <- sqrt(pmax(variances, 0))
  tolerance <- sqrt(.Machine$double.eps) *
    pmax(1, apply(abs(matrix), 2L, max))
  nonconstant <- is.finite(standard_deviations) &
    standard_deviations > tolerance
  center[nonconstant] <- means[nonconstant]
  scale[nonconstant] <- standard_deviations[nonconstant]
  list(center = center, scale = scale)
}

build_general_formula_blueprint <- function(
    formula, data, weights, standardize, name) {
  frame <- tryCatch(
    stats::model.frame(
      formula, data = data, na.action = stats::na.pass,
      drop.unused.levels = FALSE
    ),
    error = function(error) {
      stop(
        "Could not evaluate `", name, "`: ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (nrow(frame) != nrow(data) || anyNA(frame)) {
    stop(
      "`", name, "` must be evaluable without missing values on every ",
      "complete analysis row.",
      call. = FALSE
    )
  }
  matrix <- stats::model.matrix(formula, data = frame)
  storage.mode(matrix) <- "double"
  if (!ncol(matrix) || any(!is.finite(matrix))) {
    stop("`", name, "` produced an empty or non-finite model matrix.",
         call. = FALSE)
  }
  scaling <- weighted_matrix_center_scale(matrix, weights, standardize)
  transformed <- sweep(matrix, 2L, scaling$center, "-")
  transformed <- sweep(transformed, 2L, scaling$scale, "/")
  colnames(transformed) <- colnames(matrix)
  list(
    terms = formula,
    contrasts = attr(matrix, "contrasts"),
    columns = colnames(matrix),
    center = scaling$center,
    scale = scaling$scale,
    matrix = transformed,
    rank = qr(transformed)$rank,
    name = name
  )
}

predict_general_formula_blueprint <- function(blueprint, data) {
  frame <- stats::model.frame(
    blueprint$terms, data = data, na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  if (nrow(frame) != nrow(data) || anyNA(frame)) {
    stop(
      "`", blueprint$name, "` is not finite after evaluating the policy.",
      call. = FALSE
    )
  }
  matrix <- stats::model.matrix(
    blueprint$terms, data = frame,
    contrasts.arg = blueprint$contrasts
  )
  if (!identical(colnames(matrix), blueprint$columns)) {
    stop(
      "`", blueprint$name, "` changed columns after evaluating the policy.",
      call. = FALSE
    )
  }
  storage.mode(matrix) <- "double"
  matrix <- sweep(matrix, 2L, blueprint$center, "-")
  sweep(matrix, 2L, blueprint$scale, "/")
}

validate_general_log_g_bounds <- function(bounds) {
  if (!is.numeric(bounds) || length(bounds) != 2L ||
      anyNA(bounds) || any(!is.finite(bounds)) ||
      bounds[1L] >= bounds[2L]) {
    stop("`log_g_bounds` must contain two increasing finite numbers.",
         call. = FALSE)
  }
  as.numeric(bounds)
}

general_outcome_bridge_evaluation <- function(
    coefficients, matrix, link = c("logit", "identity")) {
  link <- match.arg(link)
  linear_predictor <- drop(matrix %*% coefficients)
  if (identical(link, "identity")) {
    value <- linear_predictor
    derivative <- matrix
  } else {
    value <- stats::plogis(linear_predictor)
    derivative <- value * (1 - value) * matrix
  }
  list(
    value = value,
    derivative = derivative,
    linear_predictor = linear_predictor
  )
}

general_treatment_bridge_evaluation <- function(
    coefficients, matrix, log_g_bounds) {
  linear_predictor <- drop(matrix %*% coefficients)
  bounded_predictor <- pmin(
    pmax(linear_predictor, log_g_bounds[1L]),
    log_g_bounds[2L]
  )
  value <- exp(bounded_predictor)
  active <- linear_predictor > log_g_bounds[1L] &
    linear_predictor < log_g_bounds[2L]
  derivative <- value * active * matrix
  list(
    value = value,
    derivative = derivative,
    linear_predictor = linear_predictor,
    bounded = !active
  )
}

initial_general_outcome_coefficients <- function(
    design, outcome, weights, link) {
  fit <- suppressWarnings(tryCatch(
    if (identical(link, "identity")) {
      stats::lm.wfit(design, outcome, weights)
    } else {
      stats::glm.fit(
        design, outcome, weights = weights,
        family = stats::binomial()
      )
    },
    error = function(error) NULL
  ))
  if (is.null(fit)) return(rep(0, ncol(design)))
  coefficients <- unname(fit$coefficients)
  coefficients[!is.finite(coefficients)] <- 0
  coefficients
}

normalize_general_treatment_starts <- function(
    start_g, policies, dimension) {
  number_policies <- length(policies)
  if (is.null(start_g)) {
    output <- rep(list(NULL), number_policies)
  } else if (is.numeric(start_g)) {
    output <- rep(list(start_g), number_policies)
  } else if (is.list(start_g)) {
    if (!is.null(names(start_g)) &&
        all(policies %in% names(start_g))) {
      output <- start_g[policies]
    } else {
      output <- start_g
    }
  } else {
    stop("`start_g` must be NULL, a numeric vector, or a list.",
         call. = FALSE)
  }
  if (length(output) != number_policies) {
    stop("`start_g` must provide one vector per policy.", call. = FALSE)
  }
  invalid <- vapply(output, function(value) {
    !is.null(value) &&
      (!is.numeric(value) || length(value) != dimension ||
       anyNA(value) || any(!is.finite(value)))
  }, logical(1))
  if (any(invalid)) {
    stop(
      "Every non-NULL `start_g` vector must contain ", dimension,
      " finite values.",
      call. = FALSE
    )
  }
  lapply(output, function(value) {
    if (is.null(value)) NULL else as.numeric(value)
  })
}

prepare_general_parametric_inputs <- function(
    data, treatment, outcome, covariates,
    negative_control_treatment, negative_control_outcome,
    policy, weights, target, policy_support, population_size,
    outcome_bridge, outcome_instruments,
    treatment_bridge, treatment_instruments,
    standardize, log_g_bounds) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  for (name in c("treatment", "outcome")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L) {
      stop("`", name, "` must be one column name.", call. = FALSE)
    }
  }
  for (name in c(
    "covariates", "negative_control_treatment",
    "negative_control_outcome"
  )) {
    value <- get(name)
    if (!is.character(value) || !length(value) || anyNA(value)) {
      stop("`", name, "` must contain one or more column names.",
           call. = FALSE)
    }
  }
  assert_columns(
    data,
    unique(c(
      treatment, outcome, covariates,
      negative_control_treatment, negative_control_outcome
    )),
    "general parametric analysis"
  )
  assert_flag(standardize, "standardize")
  log_g_bounds <- validate_general_log_g_bounds(log_g_bounds)
  policies <- normalize_policies(policy)
  built <- build_analysis_data(
    data, treatment, outcome, covariates,
    negative_control_treatment, negative_control_outcome,
    weights, target, policies, policy_support
  )
  analysis_data <- data[built$complete_rows, , drop = FALSE]
  if (is.null(population_size)) population_size <- sum(built$weight)
  assert_positive(population_size, "population_size")
  if (length(population_size) != 1L) {
    stop("`population_size` must be one positive value.", call. = FALSE)
  }
  defaults <- general_parametric_default_formulas(
    treatment, covariates,
    negative_control_treatment, negative_control_outcome
  )
  formulas <- list(
    outcome_bridge = normalize_general_parametric_formula(
      outcome_bridge, defaults$outcome_bridge, "outcome_bridge"
    ),
    outcome_instruments = normalize_general_parametric_formula(
      outcome_instruments, defaults$outcome_instruments,
      "outcome_instruments"
    ),
    treatment_bridge = normalize_general_parametric_formula(
      treatment_bridge, defaults$treatment_bridge, "treatment_bridge"
    ),
    treatment_instruments = normalize_general_parametric_formula(
      treatment_instruments, defaults$treatment_instruments,
      "treatment_instruments"
    )
  )
  blueprints <- list(
    outcome_bridge = build_general_formula_blueprint(
      formulas$outcome_bridge, analysis_data, built$weight,
      standardize, "outcome_bridge"
    ),
    outcome_instruments = build_general_formula_blueprint(
      formulas$outcome_instruments, analysis_data, built$weight,
      standardize, "outcome_instruments"
    ),
    treatment_bridge = build_general_formula_blueprint(
      formulas$treatment_bridge, analysis_data, built$weight,
      standardize, "treatment_bridge"
    ),
    treatment_instruments = build_general_formula_blueprint(
      formulas$treatment_instruments, analysis_data, built$weight,
      standardize, "treatment_instruments"
    )
  )
  if (ncol(blueprints$outcome_instruments$matrix) <
      ncol(blueprints$outcome_bridge$matrix)) {
    stop(
      "The outcome bridge has more coefficients than outcome instruments.",
      call. = FALSE
    )
  }
  if (ncol(blueprints$treatment_instruments$matrix) <
      ncol(blueprints$treatment_bridge$matrix)) {
    stop(
      "The treatment bridge has more coefficients than treatment instruments.",
      call. = FALSE
    )
  }
  overidentified <-
    ncol(blueprints$outcome_instruments$matrix) >
      ncol(blueprints$outcome_bridge$matrix) ||
    ncol(blueprints$treatment_instruments$matrix) >
      ncol(blueprints$treatment_bridge$matrix)
  if (standardize && overidentified) {
    stop(
      "Overidentified formula systems require `standardize = FALSE` so ",
      "the minimum-distance instrument metric is fixed rather than ",
      "estimated from the same sample.",
      call. = FALSE
    )
  }
  shifted_data <- lapply(built$q, function(shifted_treatment) {
    value <- analysis_data
    value[[treatment]] <- shifted_treatment
    value
  })
  matrices <- list(
    h = blueprints$outcome_bridge$matrix,
    h_instrument = blueprints$outcome_instruments$matrix,
    g = blueprints$treatment_bridge$matrix,
    g_instrument = blueprints$treatment_instruments$matrix,
    hq = lapply(shifted_data, function(value) {
      predict_general_formula_blueprint(
        blueprints$outcome_bridge, value
      )
    }),
    g_instrument_q = lapply(shifted_data, function(value) {
      predict_general_formula_blueprint(
        blueprints$treatment_instruments, value
      )
    })
  )
  names(matrices$hq) <- names(policies)
  names(matrices$g_instrument_q) <- names(policies)
  denominator <- sum(built$weight * built$target)
  if (!is.finite(denominator) || denominator <= 0) {
    stop("The weighted target population is empty.", call. = FALSE)
  }
  list(
    data = analysis_data,
    treatment = treatment,
    outcome = built$y,
    weights = built$weight,
    target = built$target,
    q = built$q,
    support = built$policy_support,
    policies = names(policies),
    complete_rows = built$complete_rows,
    original_n = built$original_n,
    population_size = as.numeric(population_size),
    denominator = denominator,
    target_probability = denominator / population_size,
    formulas = formulas,
    blueprints = blueprints,
    matrices = matrices,
    standardize = standardize,
    log_g_bounds = log_g_bounds
  )
}

general_outcome_moment_components <- function(
    coefficients, prepared, link) {
  evaluation <- general_outcome_bridge_evaluation(
    coefficients, prepared$matrices$h, link
  )
  instruments <- prepared$matrices$h_instrument
  moments <- instruments * (prepared$outcome - evaluation$value)
  mean_moments <- weighted_column_means(moments, prepared$weights)
  jacobian <- -crossprod(
    instruments, prepared$weights * evaluation$derivative
  ) / sum(prepared$weights)
  stationarity <- drop(crossprod(jacobian, mean_moments))
  scores <- -evaluation$derivative *
    drop(instruments %*% mean_moments) +
    moments %*% jacobian
  scores <- sweep(scores, 2L, 2 * stationarity, "-")
  list(
    moments = mean_moments,
    jacobian = jacobian,
    stationarity = stationarity,
    scores = scores,
    value = evaluation$value,
    derivative = evaluation$derivative
  )
}

general_treatment_moment_components <- function(
    coefficients, prepared, policy_index) {
  evaluation <- general_treatment_bridge_evaluation(
    coefficients, prepared$matrices$g, prepared$log_g_bounds
  )
  observed_instruments <- prepared$matrices$g_instrument
  shifted_instruments <-
    prepared$matrices$g_instrument_q[[policy_index]]
  support <- prepared$support[[policy_index]]
  moments <- prepared$target * shifted_instruments -
    support * evaluation$value * observed_instruments
  mean_moments <- weighted_column_means(moments, prepared$weights)
  jacobian <- -crossprod(
    observed_instruments,
    prepared$weights * support * evaluation$derivative
  ) / sum(prepared$weights)
  stationarity <- drop(crossprod(jacobian, mean_moments))
  scores <- -support * evaluation$derivative *
    drop(observed_instruments %*% mean_moments) +
    moments %*% jacobian
  scores <- sweep(scores, 2L, 2 * stationarity, "-")
  list(
    moments = mean_moments,
    jacobian = jacobian,
    stationarity = stationarity,
    scores = scores,
    value = evaluation$value,
    derivative = evaluation$derivative,
    bounded = evaluation$bounded,
    linear_predictor = evaluation$linear_predictor
  )
}

solve_general_outcome_bridge <- function(
    prepared, link, start_h, max_iterations) {
  dimension <- ncol(prepared$matrices$h)
  if (!is.null(start_h) &&
      (!is.numeric(start_h) || length(start_h) != dimension ||
       anyNA(start_h) || any(!is.finite(start_h)))) {
    stop(
      "`start_h` must contain ", dimension, " finite values.",
      call. = FALSE
    )
  }
  initial <- initial_general_outcome_coefficients(
    prepared$matrices$h, prepared$outcome, prepared$weights, link
  )
  starts <- list(initial, rep(0, dimension))
  if (!is.null(start_h)) starts <- c(list(as.numeric(start_h)), starts)
  solve_bridge_moments(
    function(coefficients) {
      general_outcome_moment_components(
        coefficients, prepared, link
      )$moments
    },
    starts = starts,
    max_iterations = max_iterations,
    jacobian_function = function(coefficients) {
      general_outcome_moment_components(
        coefficients, prepared, link
      )$jacobian
    }
  )
}

solve_general_treatment_bridge <- function(
    prepared, policy_index, start_g, max_iterations) {
  dimension <- ncol(prepared$matrices$g)
  zero <- rep(0, dimension)
  intercept <- match(
    "(Intercept)", prepared$blueprints$treatment_bridge$columns
  )
  if (!is.na(intercept)) {
    numerator <- sum(prepared$weights * prepared$target)
    denominator <- sum(
      prepared$weights * prepared$support[[policy_index]]
    )
    if (denominator > 0) zero[intercept] <- log(numerator / denominator)
  }
  starts <- list(zero, rep(0, dimension))
  if (!is.null(start_g)) starts <- c(list(start_g), starts)
  solve_bridge_moments(
    function(coefficients) {
      general_treatment_moment_components(
        coefficients, prepared, policy_index
      )$moments
    },
    starts = starts,
    max_iterations = max_iterations,
    jacobian_function = function(coefficients) {
      general_treatment_moment_components(
        coefficients, prepared, policy_index
      )$jacobian
    }
  )
}

general_parametric_estimate_names <- function(policies) {
  unlist(lapply(policies, function(policy) {
    paste(policy, c("OR", "DQW", "DR"), sep = "::")
  }), use.names = FALSE)
}

general_parametric_contributions <- function(
    outcome_coefficients, treatment_coefficients,
    prepared, outcome_link) {
  h0 <- general_outcome_bridge_evaluation(
    outcome_coefficients, prepared$matrices$h, outcome_link
  )$value
  hq <- lapply(prepared$matrices$hq, function(matrix) {
    general_outcome_bridge_evaluation(
      outcome_coefficients, matrix, outcome_link
    )$value
  })
  g0 <- lapply(seq_along(prepared$policies), function(index) {
    general_treatment_bridge_evaluation(
      treatment_coefficients[[index]], prepared$matrices$g,
      prepared$log_g_bounds
    )$value
  })
  names(g0) <- prepared$policies
  contributions <- do.call(cbind, lapply(
    seq_along(prepared$policies),
    function(index) {
      target_h <- prepared$target * hq[[index]]
      weighted_y <- prepared$support[[index]] *
        g0[[index]] * prepared$outcome
      correction <- prepared$support[[index]] *
        g0[[index]] * (prepared$outcome - h0)
      cbind(
        OR = target_h,
        DQW = weighted_y,
        DR = target_h + correction
      )
    }
  ))
  colnames(contributions) <- general_parametric_estimate_names(
    prepared$policies
  )
  list(
    contributions = contributions,
    h0 = h0,
    hq = do.call(cbind, hq),
    g0 = do.call(cbind, g0)
  )
}

general_parametric_parameter_map <- function(prepared) {
  outcome_dimension <- ncol(prepared$matrices$h)
  treatment_dimension <- ncol(prepared$matrices$g)
  number_policies <- length(prepared$policies)
  outcome_index <- seq_len(outcome_dimension)
  treatment_index <- lapply(seq_len(number_policies), function(index) {
    outcome_dimension +
      (index - 1L) * treatment_dimension +
      seq_len(treatment_dimension)
  })
  estimate_index <- outcome_dimension +
    number_policies * treatment_dimension +
    seq_len(3L * number_policies)
  list(
    outcome = outcome_index,
    treatment = treatment_index,
    estimates = estimate_index,
    number_parameters = max(estimate_index)
  )
}

general_parametric_estimating_system <- function(
    parameters, prepared, outcome_link) {
  map <- general_parametric_parameter_map(prepared)
  outcome_coefficients <- parameters[map$outcome]
  treatment_coefficients <- lapply(
    map$treatment, function(index) parameters[index]
  )
  estimates <- parameters[map$estimates]
  outcome <- general_outcome_moment_components(
    outcome_coefficients, prepared, outcome_link
  )
  treatment <- lapply(seq_along(prepared$policies), function(index) {
    general_treatment_moment_components(
      treatment_coefficients[[index]], prepared, index
    )
  })
  nuisance <- general_parametric_contributions(
    outcome_coefficients, treatment_coefficients,
    prepared, outcome_link
  )
  estimate_equations <- nuisance$contributions -
    outer(prepared$target, estimates)
  scores <- do.call(cbind, c(
    list(outcome$scores),
    lapply(treatment, `[[`, "scores"),
    list(estimate_equations)
  ))
  score_names <- c(
    paste0("h_minimum_distance_", seq_along(map$outcome)),
    unlist(lapply(seq_along(prepared$policies), function(index) {
      paste0(
        prepared$policies[index], "_g_minimum_distance_",
        seq_along(map$treatment[[index]])
      )
    })),
    paste0("target_", colnames(nuisance$contributions))
  )
  colnames(scores) <- score_names
  weight_scale <- sum(prepared$weights) / prepared$population_size
  system <- c(
    weight_scale * outcome$stationarity,
    unlist(lapply(treatment, function(value) {
      weight_scale * value$stationarity
    }), use.names = FALSE),
    colSums(prepared$weights * estimate_equations) /
      prepared$population_size
  )
  names(system) <- score_names
  list(
    scores = scores,
    system = system,
    nuisance = nuisance,
    outcome = outcome,
    treatment = treatment
  )
}

general_parametric_inference <- function(
    outcome_solution, treatment_solutions, estimates,
    prepared, outcome_link) {
  outcome_names <- paste0(
    "h:", prepared$blueprints$outcome_bridge$columns
  )
  treatment_names <- unlist(lapply(
    seq_along(prepared$policies), function(index) {
      paste0(
        prepared$policies[index], ":g:",
        prepared$blueprints$treatment_bridge$columns
      )
    }
  ))
  parameters <- c(
    stats::setNames(outcome_solution$coefficients, outcome_names),
    stats::setNames(
      unlist(lapply(treatment_solutions, `[[`, "coefficients")),
      treatment_names
    ),
    estimates
  )
  components <- general_parametric_estimating_system(
    parameters, prepared, outcome_link
  )
  system_function <- function(candidate) {
    general_parametric_estimating_system(
      candidate, prepared, outcome_link
    )$system
  }
  jacobian <- finite_difference_jacobian(
    system_function, parameters
  )
  inverse <- generalized_inverse(jacobian)
  if (inverse$rank < length(parameters)) {
    warning(
      "The general parametric estimating-equation Jacobian is rank ",
      "deficient; inference uses its singular-value generalized inverse.",
      call. = FALSE
    )
  }
  parameter_influence <-
    -components$scores %*% t(inverse$inverse)
  colnames(parameter_influence) <- names(parameters)
  influence <- parameter_influence[, names(estimates), drop = FALSE]
  weighted_influence <- prepared$weights * influence
  covariance <- crossprod(weighted_influence) /
    prepared$population_size^2
  dimnames(covariance) <- list(names(estimates), names(estimates))
  standard_error <- sqrt(pmax(diag(covariance), 0))
  names(standard_error) <- names(estimates)
  list(
    standard_error = standard_error,
    covariance = covariance,
    asymptotic_variance = stats::setNames(
      diag(covariance) * prepared$population_size,
      names(estimates)
    ),
    influence_function = influence,
    jacobian = jacobian,
    jacobian_rank = inverse$rank,
    jacobian_condition_number = inverse$condition_number,
    jacobian_singular_values = inverse$singular_values,
    equation_residual = system_function(parameters),
    nuisance = components$nuisance,
    treatment = components$treatment
  )
}

general_parametric_diagnostics <- function(
    outcome_solution, treatment_solutions, inference, prepared) {
  diagnostics <- rbind(
    data.frame(
      component = "outcome",
      policy = NA_character_,
      converged = outcome_solution$converged,
      solution_type = outcome_solution$solution_type,
      residual_norm = outcome_solution$residual_norm,
      stationarity_norm = outcome_solution$stationarity_norm,
      method = outcome_solution$method,
      bound_fraction = NA_real_,
      stringsAsFactors = FALSE
    ),
    do.call(rbind, lapply(seq_along(treatment_solutions), function(index) {
      solution <- treatment_solutions[[index]]
      bound_fraction <- if (is.null(inference)) {
        evaluation <- general_treatment_bridge_evaluation(
          solution$coefficients, prepared$matrices$g,
          prepared$log_g_bounds
        )
        mean(evaluation$bounded)
      } else {
        mean(inference$treatment[[index]]$bounded)
      }
      data.frame(
        component = "treatment",
        policy = prepared$policies[index],
        converged = solution$converged,
        solution_type = solution$solution_type,
        residual_norm = solution$residual_norm,
        stationarity_norm = solution$stationarity_norm,
        method = solution$method,
        bound_fraction = bound_fraction,
        stringsAsFactors = FALSE
      )
    }))
  )
  diagnostics
}

#' Fit general parametric proximal MTP estimators
#'
#' Fits fully parametric proximal outcome and treatment bridge models for one
#' or more modified treatment policies. Unlike [pmtp_parametric()], this
#' function is not tied to the paper's simulation DGP: it accepts multiple
#' measured covariates, arbitrary monotone policies, formula-defined bridge
#' models and instruments, and inverse two-phase sampling weights.
#'
#' The default outcome bridge is logistic and contains main effects of
#' treatment, measured covariates and negative-control outcomes, together with
#' a quadratic treatment term. Its default instruments replace the
#' negative-control outcomes by the negative-control treatments. The default
#' treatment bridge is positive and log-linear in treatment, measured
#' covariates and negative-control treatments. Its moment instruments replace
#' the negative-control treatments by the negative-control outcomes.
#'
#' For every policy, the treatment bridge solves moments of the form
#' `P[target * v(q(A), L, W) - support * g(A,L,Z) * v(A,L,W)] = 0`.
#' The OR, DQW and DR estimators are divided by the weighted target-population
#' mass. Standard errors use the stacked influence function and squared inverse
#' sampling weights. Exact empirical roots are used when they exist; otherwise
#' a stationary minimum-distance solution is used and labeled explicitly.
#'
#' @param data A data frame.
#' @param treatment,outcome Single treatment and outcome column names.
#' @param covariates Measured-covariate column names.
#' @param negative_control_treatment Negative-control treatment column names.
#' @param negative_control_outcome Negative-control outcome column names.
#' @param policy A policy function or named list of policy functions, with the
#'   same interface as [pmtp()].
#' @param weights Optional positive inverse-inclusion-weight column name or
#'   numeric vector. Do not normalize sampling weights to mean one.
#' @param target Optional target-population column name or 0/1 vector.
#' @param policy_support Optional support indicator or list of indicators. See
#'   [pmtp()].
#' @param population_size Known phase-one population size. When omitted, uses
#'   the sum of the supplied weights.
#' @param outcome_bridge Formula for the outcome bridge design
#'   `h(A,L,W)`. The default contains treatment, squared treatment, measured
#'   covariates, and negative-control outcomes.
#' @param outcome_instruments Formula for the outcome bridge moments. The
#'   number of resulting columns must be at least the number of outcome-bridge
#'   coefficients.
#' @param treatment_bridge Formula for the positive log-linear treatment
#'   bridge `g(A,L,Z)`.
#' @param treatment_instruments Formula defining `v(A,L,W)` in the treatment
#'   bridge moments. The number of resulting columns must be at least the
#'   number of treatment-bridge coefficients.
#' @param outcome_link Outcome-bridge inverse link, either `"logit"` or
#'   `"identity"`. A logistic link is convenient for binary outcomes but is a
#'   working bridge restriction, not an ordinary outcome regression.
#' @param standardize Logical; standardize nonconstant model-matrix columns
#'   using weighted means and standard deviations. The same transformation is
#'   reused at policy-shifted treatment values. Overidentified systems require
#'   `standardize = FALSE` so their minimum-distance metric is fixed.
#' @param log_g_bounds Two finite bounds used to prevent numerical overflow in
#'   the log-linear treatment bridge. A positive `bound_fraction` diagnostic
#'   indicates that the fitted model reached a bound.
#' @param start_h Optional outcome-bridge starting coefficients.
#' @param start_g Optional treatment-bridge starting coefficients, or a list
#'   with one vector per policy.
#' @param max_iterations Maximum numerical-solver iterations.
#' @param conf_level Confidence level for reported intervals.
#' @param strict Logical; if `TRUE`, stop when any bridge has neither an exact
#'   root nor a stationary minimum-distance solution. If `FALSE`, return point
#'   estimates and convergence diagnostics with unavailable standard errors.
#'
#' @return An object of class `pmtp_parametric_general_fit` containing the
#'   three estimators per policy, bridge coefficients and solution diagnostics,
#'   model-matrix specifications, nuisance predictions, influence functions,
#'   and the joint covariance matrix across policies and estimators.
#' @export
pmtp_parametric_general <- function(
    data,
    treatment,
    outcome,
    covariates,
    negative_control_treatment,
    negative_control_outcome,
    policy,
    weights = NULL,
    target = NULL,
    policy_support = NULL,
    population_size = NULL,
    outcome_bridge = NULL,
    outcome_instruments = NULL,
    treatment_bridge = NULL,
    treatment_instruments = NULL,
    outcome_link = c("logit", "identity"),
    standardize = TRUE,
    log_g_bounds = c(-30, 30),
    start_h = NULL,
    start_g = NULL,
    max_iterations = 500L,
    conf_level = 0.95,
    strict = TRUE) {
  outcome_link <- match.arg(outcome_link)
  assert_flag(strict, "strict")
  if (length(max_iterations) != 1L || is.na(max_iterations) ||
      max_iterations < 1L || max_iterations != as.integer(max_iterations)) {
    stop("`max_iterations` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must lie strictly between zero and one.",
         call. = FALSE)
  }
  prepared <- prepare_general_parametric_inputs(
    data, treatment, outcome, covariates,
    negative_control_treatment, negative_control_outcome,
    policy, weights, target, policy_support, population_size,
    outcome_bridge, outcome_instruments,
    treatment_bridge, treatment_instruments,
    standardize, log_g_bounds
  )
  treatment_starts <- normalize_general_treatment_starts(
    start_g, prepared$policies, ncol(prepared$matrices$g)
  )
  outcome_solution <- solve_general_outcome_bridge(
    prepared, outcome_link, start_h, as.integer(max_iterations)
  )
  treatment_solutions <- lapply(
    seq_along(prepared$policies), function(index) {
      solve_general_treatment_bridge(
        prepared, index, treatment_starts[[index]],
        as.integer(max_iterations)
      )
    }
  )
  names(treatment_solutions) <- prepared$policies
  all_converged <- outcome_solution$converged &&
    all(vapply(
      treatment_solutions, `[[`, logical(1L), "converged"
    ))
  if (!all_converged && strict) {
    failed <- c(
      if (!outcome_solution$converged) "outcome" else character(),
      prepared$policies[!vapply(
        treatment_solutions, `[[`, logical(1L), "converged"
      )]
    )
    stop(
      "No valid root or stationary minimum-distance solution for: ",
      paste(failed, collapse = ", "),
      ". Refit with `strict = FALSE` to inspect the failed solutions.",
      call. = FALSE
    )
  }
  treatment_coefficients <- lapply(
    treatment_solutions, `[[`, "coefficients"
  )
  nuisance <- general_parametric_contributions(
    outcome_solution$coefficients, treatment_coefficients,
    prepared, outcome_link
  )
  estimates <- colSums(
    prepared$weights * nuisance$contributions
  ) / prepared$denominator
  inference <- if (all_converged) {
    general_parametric_inference(
      outcome_solution, treatment_solutions, estimates,
      prepared, outcome_link
    )
  } else {
    NULL
  }
  standard_error <- if (is.null(inference)) {
    stats::setNames(rep(NA_real_, length(estimates)), names(estimates))
  } else {
    inference$standard_error
  }
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  split_names <- strsplit(names(estimates), "::", fixed = TRUE)
  estimate_table <- data.frame(
    policy = vapply(split_names, `[[`, character(1L), 1L),
    estimator = vapply(split_names, `[[`, character(1L), 2L),
    estimate = unname(estimates),
    std_error = unname(standard_error),
    conf_low = unname(estimates - critical * standard_error),
    conf_high = unname(estimates + critical * standard_error),
    valid = all_converged,
    stringsAsFactors = FALSE
  )
  diagnostics <- general_parametric_diagnostics(
    outcome_solution, treatment_solutions, inference, prepared
  )
  if (any(diagnostics$bound_fraction > 0, na.rm = TRUE)) {
    warning(
      "At least one treatment bridge reached `log_g_bounds`; inspect ",
      "`diagnostics$bound_fraction` and consider revising the model.",
      call. = FALSE
    )
  }
  coefficients <- list(
    outcome = stats::setNames(
      outcome_solution$coefficients,
      prepared$blueprints$outcome_bridge$columns
    ),
    treatment = lapply(treatment_solutions, function(solution) {
      stats::setNames(
        solution$coefficients,
        prepared$blueprints$treatment_bridge$columns
      )
    })
  )
  covariance <- if (is.null(inference)) {
    matrix(
      NA_real_, nrow = length(estimates), ncol = length(estimates),
      dimnames = list(names(estimates), names(estimates))
    )
  } else {
    inference$covariance
  }
  structure(list(
    estimates = estimate_table,
    estimate_vector = estimates,
    standard_error = standard_error,
    covariance = covariance,
    asymptotic_variance = if (is.null(inference)) {
      stats::setNames(rep(NA_real_, length(estimates)), names(estimates))
    } else {
      inference$asymptotic_variance
    },
    influence_function = if (is.null(inference)) NULL else {
      inference$influence_function
    },
    inference = if (is.null(inference)) NULL else inference[c(
      "jacobian", "jacobian_rank", "jacobian_condition_number",
      "jacobian_singular_values", "equation_residual"
    )],
    coefficients = coefficients,
    bridge_solutions = list(
      outcome = outcome_solution,
      treatment = treatment_solutions
    ),
    diagnostics = diagnostics,
    nuisance = list(
      h0 = nuisance$h0,
      hq = nuisance$hq,
      g0 = nuisance$g0
    ),
    contributions = nuisance$contributions,
    formulas = prepared$formulas,
    model_matrices = lapply(prepared$blueprints, function(value) {
      value[c(
        "columns", "center", "scale", "rank", "name"
      )]
    }),
    outcome_link = outcome_link,
    weights = prepared$weights,
    target = prepared$target,
    policy_support = prepared$support,
    complete_rows = prepared$complete_rows,
    population_size = prepared$population_size,
    target_probability = prepared$target_probability,
    weighted = !is.null(weights),
    standardize = standardize,
    log_g_bounds = prepared$log_g_bounds,
    valid = all_converged,
    call = match.call()
  ), class = "pmtp_parametric_general_fit")
}

#' @export
print.pmtp_parametric_general_fit <- function(x, ...) {
  cat("General parametric proximal MTP estimates\n")
  print(x$estimates, row.names = FALSE)
  cat("\nBridge-solution diagnostics\n")
  print(x$diagnostics, row.names = FALSE)
  invisible(x)
}

#' @export
coef.pmtp_parametric_general_fit <- function(object, ...) {
  object$estimate_vector
}

#' @export
vcov.pmtp_parametric_general_fit <- function(object, ...) {
  object$covariance
}

#' @export
summary.pmtp_parametric_general_fit <- function(
    object, conf_level = 0.95, ...) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must lie strictly between zero and one.",
         call. = FALSE)
  }
  output <- object$estimates
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  output$conf_low <- output$estimate - critical * output$std_error
  output$conf_high <- output$estimate + critical * output$std_error
  attr(output, "conf_level") <- conf_level
  attr(output, "diagnostics") <- object$diagnostics
  class(output) <- c(
    "summary_pmtp_parametric_general_fit", class(output)
  )
  output
}

#' @export
print.summary_pmtp_parametric_general_fit <- function(x, ...) {
  level <- format(100 * attr(x, "conf_level"), trim = TRUE)
  cat("General parametric proximal MTP estimates with ", level, "% CIs\n",
      sep = "")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}
