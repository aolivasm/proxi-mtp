#' Control tuning and numerical behavior for proximal MTP estimation
#'
#' The default grids reproduce the scaling-factor grids described in the
#' preprint. They can be expensive because every configuration is evaluated in
#' nested cross-validation. Use [pmtp_control_fast()] for examples and smoke
#' tests, but not as the final tuning strategy for a scientific analysis.
#'
#' @param outer_folds Number of outer cross-fitting folds. Must be at least 2.
#' @param inner_folds Number of inner tuning folds within each outer-training
#'   sample. Must be at least 2 when `tune = TRUE`; ignored otherwise.
#' @param tune Logical; perform inner cross-validation. When `FALSE`, every
#'   tuning field must be scalar and those values are used directly.
#' @param lambda_h,lambda_g Scaling factors for the outer RKHS penalties for the
#'   outcome and treatment bridges. The base rate is computed from the number
#'   of observed rows in the relevant fold, including under two-phase sampling;
#'   it is not computed from the sum of inverse-probability weights.
#' @param lambda_gp,lambda_hp Scaling factors for the inner adversarial RKHS
#'   penalties for the outcome and treatment bridge objectives, using the same
#'   observed-fold sample-size convention.
#' @param bandwidth_h,bandwidth_g Multipliers for the median-heuristic external
#'   kernel variances.
#' @param bandwidth_gp,bandwidth_hp Multipliers for the median-heuristic
#'   adversarial kernel variances.
#' @param risk_bandwidth Multiplier for validation-risk kernel variances. The
#'   base bandwidth is recomputed on each validation fold.
#' @param risk_penalty Scaling factor for `log(n) / n` in validation risks.
#' @param max_norm_h,max_norm_g Optional RKHS norm bounds. Use `Inf` to disable.
#' @param seed Random seed used for fold creation without permanently changing
#'   the caller's random-number state.
#' @param jitter Initial scale-aware diagonal jitter for numerical solves.
#' @param max_solve_tries Maximum number of increasing-jitter Cholesky attempts.
#' @param keep_cv Logical; retain candidate-level cross-validation results.
#' @param progress Logical; print fold progress.
#'
#' @return A list of class `pmtp_control`.
#' @export
pmtp_control <- function(
    outer_folds = 3L,
    inner_folds = 3L,
    tune = TRUE,
    lambda_h = 10^seq(-5, -1, by = 1),
    lambda_gp = 10^seq(-1, 2, by = 1),
    lambda_g = 10^seq(-5, -1, by = 1),
    lambda_hp = 10^seq(-1, 2, by = 1),
    bandwidth_h = c(1 / 4, 1 / 2, 1, 2, 4),
    bandwidth_gp = 1 / 4,
    bandwidth_g = c(1 / 4, 1 / 2, 1, 2, 4),
    bandwidth_hp = 1 / 4,
    risk_bandwidth = 1 / 4,
    risk_penalty = 1,
    max_norm_h = 50,
    max_norm_g = 50,
    seed = 1234L,
    jitter = 1e-8,
    max_solve_tries = 8L,
    keep_cv = TRUE,
    progress = FALSE) {
  assert_flag(tune, "tune")
  if (outer_folds < 2L || (tune && inner_folds < 2L)) {
    stop("Outer folds must be at least 2, and inner folds must be at least 2 when tuning.",
         call. = FALSE)
  }
  if (!tune) {
    tuning_values <- list(
      lambda_h, lambda_gp, lambda_g, lambda_hp,
      bandwidth_h, bandwidth_gp, bandwidth_g, bandwidth_hp
    )
    if (any(lengths(tuning_values) != 1L)) {
      stop("All penalty and bandwidth fields must be scalar when `tune = FALSE`.",
           call. = FALSE)
    }
  }
  for (name in c(
    "lambda_h", "lambda_gp", "lambda_g", "lambda_hp",
    "bandwidth_h", "bandwidth_gp", "bandwidth_g", "bandwidth_hp",
    "risk_bandwidth", "risk_penalty", "jitter"
  )) {
    assert_positive(get(name), name)
  }
  assert_positive(max_norm_h, "max_norm_h", allow_inf = TRUE)
  assert_positive(max_norm_g, "max_norm_g", allow_inf = TRUE)
  assert_flag(keep_cv, "keep_cv")
  assert_flag(progress, "progress")

  structure(list(
    outer_folds = as.integer(outer_folds),
    inner_folds = as.integer(inner_folds),
    tune = tune,
    lambda_h = lambda_h,
    lambda_gp = lambda_gp,
    lambda_g = lambda_g,
    lambda_hp = lambda_hp,
    bandwidth_h = bandwidth_h,
    bandwidth_gp = bandwidth_gp,
    bandwidth_g = bandwidth_g,
    bandwidth_hp = bandwidth_hp,
    risk_bandwidth = risk_bandwidth,
    risk_penalty = risk_penalty,
    max_norm_h = max_norm_h,
    max_norm_g = max_norm_g,
    seed = as.integer(seed),
    jitter = jitter,
    max_solve_tries = as.integer(max_solve_tries),
    keep_cv = keep_cv,
    progress = progress
  ), class = "pmtp_control")
}

#' Fixed-hyperparameter control for diagnostic runs
#'
#' Skips inner cross-validation and uses one specified set of penalty and
#' bandwidth multipliers. This is intended for algebraic debugging and timing
#' comparisons, not final scientific analyses.
#'
#' @inheritParams pmtp_control
#' @return A list of class `pmtp_control`.
#' @export
pmtp_control_fixed <- function(
    outer_folds = 2L,
    lambda_h = 1e-3,
    lambda_gp = 1,
    lambda_g = 1e-3,
    lambda_hp = 1,
    bandwidth_h = 1,
    bandwidth_gp = 1 / 4,
    bandwidth_g = 1,
    bandwidth_hp = 1 / 4,
    risk_bandwidth = 1 / 4,
    risk_penalty = 1,
    max_norm_h = 50,
    max_norm_g = 50,
    seed = 1234L,
    jitter = 1e-8,
    max_solve_tries = 8L,
    keep_cv = FALSE,
    progress = FALSE) {
  pmtp_control(
    outer_folds = outer_folds,
    inner_folds = 1L,
    tune = FALSE,
    lambda_h = lambda_h,
    lambda_gp = lambda_gp,
    lambda_g = lambda_g,
    lambda_hp = lambda_hp,
    bandwidth_h = bandwidth_h,
    bandwidth_gp = bandwidth_gp,
    bandwidth_g = bandwidth_g,
    bandwidth_hp = bandwidth_hp,
    risk_bandwidth = risk_bandwidth,
    risk_penalty = risk_penalty,
    max_norm_h = max_norm_h,
    max_norm_g = max_norm_g,
    seed = seed,
    jitter = jitter,
    max_solve_tries = max_solve_tries,
    keep_cv = keep_cv,
    progress = progress
  )
}

#' A small control grid for examples and smoke tests
#'
#' This grid is intentionally too small for a final scientific analysis.
#'
#' @inheritParams pmtp_control
#' @return A list of class `pmtp_control`.
#' @export
pmtp_control_fast <- function(outer_folds = 3L, inner_folds = 2L,
                              seed = 1234L, progress = FALSE) {
  pmtp_control(
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    lambda_h = c(1e-3, 1e-2),
    lambda_gp = c(1, 10),
    lambda_g = c(1e-3, 1e-2),
    lambda_hp = c(1, 10),
    bandwidth_h = c(1 / 2, 1),
    bandwidth_gp = 1 / 4,
    bandwidth_g = c(1 / 2, 1),
    bandwidth_hp = 1 / 4,
    seed = seed,
    progress = progress
  )
}
