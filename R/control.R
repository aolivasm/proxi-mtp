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
#' @param inner_repeats Number of independently randomized inner-CV partitions.
#'   The default of one reproduces ordinary inner cross-validation. Values
#'   greater than one average validation risk over all repeat-by-fold
#'   evaluations while preserving the training fraction of `inner_folds`.
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
#' @param risk_penalty Scaling factor for the squared critical-radius rate in
#'   validation risks.
#' @param critical_radius_rule Penalty-rate rule. `"legacy_d1"` preserves the
#'   preprint implementation, using `log(n)` as though every Gaussian RKHS had
#'   dimension one. The `"gaussian_dimension"` rule used in the revised paper uses
#'   `log(n)^d`, where `d` is the number of columns in the relevant bridge or
#'   adversarial-kernel input. The number of observed rows in each fold remains
#'   the sample-size argument. The `"matern_sobolev"` rule uses the polynomial
#'   critical-radius rate determined by the input dimension and `sobolev_l`.
#' @param selection_rule Cross-validation selection rule. `"minimum"` selects
#'   the smallest mean validation risk. The experimental
#'   `"one_se_regularized"` rule forms a one-standard-error set using the
#'   foldwise risk standard error at the minimum, then lexicographically
#'   prefers larger outer and inner penalties followed by larger external and
#'   adversarial kernel bandwidths; remaining ties use the smallest mean risk.
#'   The experimental
#'   `"one_se_interior"`
#'   rule instead prefers the fewest boundary coordinates within the same
#'   one-standard-error set, then uses the smallest mean risk.
#' @param max_norm_h,max_norm_g Optional RKHS norm bounds. Use `Inf` to disable.
#' @param seed Random seed used for fold creation without permanently changing
#'   the caller's random-number state.
#' @param jitter Initial scale-aware diagonal jitter for numerical solves.
#' @param max_solve_tries Maximum number of increasing-jitter Cholesky attempts.
#' @param keep_cv Logical; retain candidate-level cross-validation results.
#' @param progress Logical; print fold progress.
#' @param kernel_approximation Kernel backend. `"exact"` uses the complete
#'   fold-specific Gram matrices. `"nystrom"` uses fold-specific Nystrom
#'   features for bridge fitting and validation risks.
#' @param nystrom_rank Positive fixed rank or a function of the number of rows
#'   in the current fold. The default [pmtp_nystrom_rank()] grows as
#'   `2 * n^(2/3)`, subject to a minimum of 30 and the fold size.
#' @param nystrom_landmarks Landmark sampling scheme. `"uniform"` samples rows
#'   uniformly; `"weighted"` samples in proportion to the analysis weights.
#' @param cache_kernel_features Logical; reuse fold-specific Nystrom feature
#'   maps and shared minimax systems across tuning candidates with the same
#'   bandwidths and inner penalty. This changes only computation, not
#'   landmarks, risks, or fitted estimators. It is ignored by the exact
#'   backend.
#' @param kernel_family RKHS kernel family. `"gaussian"` uses the Gaussian RBF
#'   kernel. `"matern_sobolev"` uses a Matérn kernel whose native space is
#'   equivalent to a Sobolev space.
#' @param matern_smoothness Positive Matérn smoothness parameter.
#' @param sobolev_l Positive spectral-order parameter used by the
#'   `"matern_sobolev"` critical-radius rule. It equals twice
#'   `matern_smoothness`.
#' @param weighted_loss_normalization Normalization for weighted empirical
#'   bridge losses. `"hajek"` divides each fold-specific weighted sum by the
#'   sum of its inverse-probability weights and preserves earlier analyses.
#'   `"horvitz_thompson"`, used in the revised paper, divides by the
#'   corresponding phase-one fold size, matching the estimator displayed in
#'   the paper. The latter requires `population_size` in [pmtp()] for an exact
#'   full-sample normalization; fold sizes are allocated in proportion to the
#'   number of observed phase-two rows in each split.
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
    critical_radius_rule = c(
      "legacy_d1", "gaussian_dimension", "matern_sobolev"
    ),
    max_norm_h = 50,
    max_norm_g = 50,
    seed = 1234L,
    jitter = 1e-8,
    max_solve_tries = 8L,
    keep_cv = TRUE,
    progress = FALSE,
    selection_rule = c(
      "minimum", "one_se_regularized", "one_se_interior"
    ),
    inner_repeats = 1L,
    kernel_approximation = c("exact", "nystrom"),
    nystrom_rank = pmtp_nystrom_rank(),
    nystrom_landmarks = c("uniform", "weighted"),
    cache_kernel_features = TRUE,
    kernel_family = c("gaussian", "matern_sobolev"),
    matern_smoothness = 2,
    sobolev_l = 4,
    weighted_loss_normalization = c("hajek", "horvitz_thompson")) {
  assert_flag(tune, "tune")
  selection_rule <- match.arg(selection_rule)
  critical_radius_rule <- match.arg(critical_radius_rule)
  kernel_approximation <- match.arg(kernel_approximation)
  nystrom_landmarks <- match.arg(nystrom_landmarks)
  kernel_family <- match.arg(kernel_family)
  weighted_loss_normalization <- match.arg(weighted_loss_normalization)
  if (length(inner_repeats) != 1L || is.na(inner_repeats) ||
      inner_repeats < 1L || inner_repeats != as.integer(inner_repeats)) {
    stop("`inner_repeats` must be a positive integer.", call. = FALSE)
  }
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
  assert_positive(matern_smoothness, "matern_smoothness")
  assert_positive(sobolev_l, "sobolev_l")
  if (length(matern_smoothness) != 1L || length(sobolev_l) != 1L) {
    stop("Kernel smoothness parameters must be scalar.", call. = FALSE)
  }
  assert_flag(keep_cv, "keep_cv")
  assert_flag(progress, "progress")
  assert_flag(cache_kernel_features, "cache_kernel_features")
  resolve_nystrom_rank(nystrom_rank, 10L)

  structure(list(
    outer_folds = as.integer(outer_folds),
    inner_folds = as.integer(inner_folds),
    inner_repeats = as.integer(inner_repeats),
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
    critical_radius_rule = critical_radius_rule,
    selection_rule = selection_rule,
    max_norm_h = max_norm_h,
    max_norm_g = max_norm_g,
    seed = as.integer(seed),
    jitter = jitter,
    max_solve_tries = as.integer(max_solve_tries),
    keep_cv = keep_cv,
    progress = progress,
    kernel_approximation = kernel_approximation,
    nystrom_rank = nystrom_rank,
    nystrom_landmarks = nystrom_landmarks,
    cache_kernel_features = cache_kernel_features,
    kernel_family = kernel_family,
    matern_smoothness = matern_smoothness,
    sobolev_l = sobolev_l,
    weighted_loss_normalization = weighted_loss_normalization
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
    critical_radius_rule = c(
      "legacy_d1", "gaussian_dimension", "matern_sobolev"
    ),
    max_norm_h = 50,
    max_norm_g = 50,
    seed = 1234L,
    jitter = 1e-8,
    max_solve_tries = 8L,
    keep_cv = FALSE,
    progress = FALSE,
    selection_rule = "minimum",
    kernel_approximation = c("exact", "nystrom"),
    nystrom_rank = pmtp_nystrom_rank(),
    nystrom_landmarks = c("uniform", "weighted"),
    cache_kernel_features = TRUE,
    kernel_family = c("gaussian", "matern_sobolev"),
    matern_smoothness = 2,
    sobolev_l = 4,
    weighted_loss_normalization = c("hajek", "horvitz_thompson")) {
  kernel_approximation <- match.arg(kernel_approximation)
  nystrom_landmarks <- match.arg(nystrom_landmarks)
  critical_radius_rule <- match.arg(critical_radius_rule)
  kernel_family <- match.arg(kernel_family)
  weighted_loss_normalization <- match.arg(weighted_loss_normalization)
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
    critical_radius_rule = critical_radius_rule,
    selection_rule = selection_rule,
    max_norm_h = max_norm_h,
    max_norm_g = max_norm_g,
    seed = seed,
    jitter = jitter,
    max_solve_tries = max_solve_tries,
    keep_cv = keep_cv,
    progress = progress,
    kernel_approximation = kernel_approximation,
    nystrom_rank = nystrom_rank,
    nystrom_landmarks = nystrom_landmarks,
    cache_kernel_features = cache_kernel_features,
    kernel_family = kernel_family,
    matern_smoothness = matern_smoothness,
    sobolev_l = sobolev_l,
    weighted_loss_normalization = weighted_loss_normalization
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
                              seed = 1234L, progress = FALSE,
                              selection_rule = "minimum",
                              inner_repeats = 1L,
                              critical_radius_rule = c(
                                "legacy_d1", "gaussian_dimension",
                                "matern_sobolev"
                              ),
                              kernel_approximation = c("exact", "nystrom"),
                              nystrom_rank = pmtp_nystrom_rank(),
                              nystrom_landmarks = c("uniform", "weighted"),
                              cache_kernel_features = TRUE,
                              kernel_family = c(
                                "gaussian", "matern_sobolev"
                              ),
                              matern_smoothness = 2,
                              sobolev_l = 4,
                              weighted_loss_normalization = c(
                                "hajek", "horvitz_thompson"
                              )) {
  kernel_approximation <- match.arg(kernel_approximation)
  nystrom_landmarks <- match.arg(nystrom_landmarks)
  critical_radius_rule <- match.arg(critical_radius_rule)
  kernel_family <- match.arg(kernel_family)
  weighted_loss_normalization <- match.arg(weighted_loss_normalization)
  pmtp_control(
    outer_folds = outer_folds,
    inner_folds = inner_folds,
    inner_repeats = inner_repeats,
    lambda_h = c(1e-3, 1e-2),
    lambda_gp = c(1, 10),
    lambda_g = c(1e-3, 1e-2),
    lambda_hp = c(1, 10),
    bandwidth_h = c(1 / 2, 1),
    bandwidth_gp = 1 / 4,
    bandwidth_g = c(1 / 2, 1),
    bandwidth_hp = 1 / 4,
    selection_rule = selection_rule,
    critical_radius_rule = critical_radius_rule,
    seed = seed,
    progress = progress,
    kernel_approximation = kernel_approximation,
    nystrom_rank = nystrom_rank,
    nystrom_landmarks = nystrom_landmarks,
    cache_kernel_features = cache_kernel_features,
    kernel_family = kernel_family,
    matern_smoothness = matern_smoothness,
    sobolev_l = sobolev_l,
    weighted_loss_normalization = weighted_loss_normalization
  )
}
