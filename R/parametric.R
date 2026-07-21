weighted_column_means <- function(x, weights) {
  colSums(weights * x) / sum(weights)
}

outcome_bridge_moments <- function(phi, data, weights, truncated_variance) {
  h_value <- pmtp_parametric_h_value(
    data$A, data$L, data$W, phi, truncated_variance
  )
  instruments <- cbind(1, data$A, data$L, data$Z, data$A^2)
  weighted_column_means(instruments * (data$Y - h_value), weights)
}

treatment_bridge_moments <- function(eta, data, weights, spec) {
  target <- pmtp_policy_target(data$A, spec)
  q_a <- pmtp_taper_policy(data$A, spec)
  q_for_moment <- ifelse(target == 1, q_a, 0)
  g_value <- pmtp_parametric_g_value(data$A, data$L, data$Z, eta, spec)
  shifted_instruments <- cbind(target, target * q_for_moment, target * data$L,
                               target * data$W)
  observed_instruments <- cbind(1, data$A, data$L, data$W)
  weighted_column_means(
    shifted_instruments - observed_instruments * g_value,
    weights
  )
}

unique_starts <- function(starts) {
  keys <- vapply(starts, function(x) paste(signif(x, 12), collapse = "|"), character(1))
  starts[!duplicated(keys)]
}

solve_bridge_moments <- function(moment_function, starts, max_iterations) {
  starts <- unique_starts(starts)
  candidates <- vector("list", 0L)
  add_candidate <- function(coefficients, method, convergence) {
    moments <- tryCatch(moment_function(coefficients), error = function(e) rep(Inf, length(coefficients)))
    if (length(moments) == length(coefficients) && all(is.finite(moments))) {
      candidates[[length(candidates) + 1L]] <<- list(
        coefficients = coefficients,
        moments = moments,
        residual_norm = sqrt(sum(moments^2)),
        method = method,
        convergence = convergence
      )
    }
  }

  for (start_index in seq_along(starts)) {
    start <- starts[[start_index]]
    root <- tryCatch(
      nleqslv::nleqslv(
        start,
        moment_function,
        method = "Broyden",
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
      add_candidate(root$x, paste0("nleqslv-start-", start_index), root$termcd)
    }
    add_candidate(start, paste0("initial-start-", start_index), NA_integer_)
  }

  if (!length(candidates)) {
    stop("No finite bridge estimating-equation candidate was obtained.", call. = FALSE)
  }
  best_index <- which.min(vapply(candidates, `[[`, numeric(1), "residual_norm"))
  best <- candidates[[best_index]]
  optimized <- tryCatch(
    stats::optim(
      best$coefficients,
      function(x) {
        moments <- moment_function(x)
        if (any(!is.finite(moments))) return(.Machine$double.xmax / 100)
        sum(moments^2)
      },
      method = "BFGS",
      control = list(maxit = as.integer(max_iterations), reltol = 1e-12)
    ),
    error = function(e) NULL
  )
  if (!is.null(optimized)) {
    add_candidate(optimized$par, "optim-BFGS", optimized$convergence)
    best_index <- which.min(vapply(candidates, `[[`, numeric(1), "residual_norm"))
    best <- candidates[[best_index]]
  }
  best$converged <- is.finite(best$residual_norm) && best$residual_norm <= 1e-6
  best
}

initial_outcome_coefficients <- function(data) {
  design <- cbind(1, data$A, data$L, data$W, data$A^2)
  fit <- suppressWarnings(tryCatch(
    stats::glm.fit(design, data$Y, family = stats::binomial()),
    error = function(e) NULL
  ))
  if (is.null(fit)) return(rep(0, 5L))
  coefficients <- fit$coefficients
  coefficients[!is.finite(coefficients)] <- 0
  unname(coefficients)
}

#' Fit the paper's proximal parametric bridge estimators
#'
#' Solves the estimating equations in Supplement C.3 and returns the proximal
#' outcome-regression, density-quotient-weighted, and doubly robust point
#' estimates. This function currently provides point-estimation diagnostics;
#' influence-function standard errors for estimated parametric nuisances will
#' be added after the estimating equations have been validated against the
#' paper DGP.
#'
#' @param data A data frame containing scalar columns `Y`, `A`, `L`, `Z`, and
#'   `W`.
#' @param spec A specification created by [pmtp_dgp_spec()]. Only its known
#'   policy parameters are used by the estimating equations.
#' @param weights Optional positive observation weights or the name of a column
#'   in `data` containing them.
#' @param start_h,start_g Optional starting coefficients for the outcome and
#'   treatment bridge estimating equations.
#' @param max_iterations Maximum iterations for each numerical solver.
#'
#' @return An object of class `pmtp_parametric_fit` containing estimates,
#'   coefficients, moment residuals, convergence information, and nuisance
#'   values.
#' @export
pmtp_parametric <- function(data, spec = pmtp_dgp_spec(), weights = NULL,
                            start_h = NULL, start_g = NULL,
                            max_iterations = 500L) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  spec <- validate_dgp_spec(spec)
  assert_columns(data, c("Y", "A", "L", "Z", "W"), "parametric data")
  if (!all(vapply(data[c("Y", "A", "L", "Z", "W")], is.numeric, logical(1))) ||
      anyNA(data[c("Y", "A", "L", "Z", "W")])) {
    stop("`Y`, `A`, `L`, `Z`, and `W` must be complete numeric columns.",
         call. = FALSE)
  }
  analysis_weights <- resolve_vector_argument(data, weights, "weights", default = 1)
  assert_positive(analysis_weights, "weights")
  if (length(max_iterations) != 1L || is.na(max_iterations) ||
      max_iterations < 1 || max_iterations != as.integer(max_iterations)) {
    stop("`max_iterations` must be a positive integer.", call. = FALSE)
  }
  bridge_parameters <- pmtp_oracle_bridge_parameters(spec)

  if (!is.null(start_h) && (length(start_h) != 5L || any(!is.finite(start_h)))) {
    stop("`start_h` must contain five finite values.", call. = FALSE)
  }
  if (!is.null(start_g) && (length(start_g) != 4L || any(!is.finite(start_g)))) {
    stop("`start_g` must contain four finite values.", call. = FALSE)
  }
  h_starts <- list(initial_outcome_coefficients(data), rep(0, 5L))
  g_starts <- list(rep(0, 4L), c(0.1, 0, -0.1, -0.1))
  if (!is.null(start_h)) h_starts <- c(list(as.numeric(start_h)), h_starts)
  if (!is.null(start_g)) g_starts <- c(list(as.numeric(start_g)), g_starts)

  h_solution <- solve_bridge_moments(
    function(phi) outcome_bridge_moments(
      phi, data, analysis_weights, bridge_parameters$truncated_variance
    ),
    h_starts,
    max_iterations
  )
  g_solution <- solve_bridge_moments(
    function(eta) treatment_bridge_moments(eta, data, analysis_weights, spec),
    g_starts,
    max_iterations
  )

  q_a <- pmtp_taper_policy(data$A, spec)
  target <- pmtp_policy_target(data$A, spec)
  support <- pmtp_policy_support(data$A, spec)
  h0 <- pmtp_parametric_h_value(
    data$A, data$L, data$W, h_solution$coefficients,
    bridge_parameters$truncated_variance
  )
  hq <- numeric(nrow(data))
  hq[target == 1] <- pmtp_parametric_h_value(
    q_a[target == 1], data$L[target == 1], data$W[target == 1],
    h_solution$coefficients, bridge_parameters$truncated_variance
  )
  g0 <- pmtp_parametric_g_value(
    data$A, data$L, data$Z, g_solution$coefficients, spec
  )
  denominator <- sum(analysis_weights * target)
  contributions <- cbind(
    OR = target * hq,
    DQW = g0 * data$Y,
    DR = target * hq + support * g0 * (data$Y - h0)
  )

  structure(list(
    estimates = colSums(analysis_weights * contributions) / denominator,
    coefficients = list(phi = h_solution$coefficients, eta = g_solution$coefficients),
    moments = list(h = h_solution$moments, g = g_solution$moments),
    residual_norm = c(h = h_solution$residual_norm, g = g_solution$residual_norm),
    converged = c(h = h_solution$converged, g = g_solution$converged),
    solver = list(h = h_solution[c("method", "convergence")],
                  g = g_solution[c("method", "convergence")]),
    nuisance = list(h0 = h0, hq = hq, g0 = g0),
    contributions = contributions,
    weights = analysis_weights,
    spec = spec
  ), class = "pmtp_parametric_fit")
}

#' @export
print.pmtp_parametric_fit <- function(x, ...) {
  cat("Proximal parametric MTP estimates\n")
  print(x$estimates)
  cat("\nEstimating-equation residual norms\n")
  print(x$residual_norm)
  invisible(x)
}

#' @export
coef.pmtp_parametric_fit <- function(object, ...) {
  object$estimates
}
