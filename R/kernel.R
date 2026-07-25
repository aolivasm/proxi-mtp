#' Gaussian radial basis function kernel
#'
#' @param x Numeric matrix with observations in rows.
#' @param y Optional numeric matrix. If omitted, computes the Gram matrix of
#'   `x`.
#' @param sigma2 Positive kernel variance in
#'   `exp(-squared_distance / (2 * sigma2))`.
#'
#' @return A numeric kernel matrix.
#' @export
gaussian_kernel <- function(x, y = NULL, sigma2 = 1) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!nrow(x) || !ncol(x) || anyNA(x) || any(!is.finite(x))) {
    stop("`x` must be a finite, nonempty numeric matrix.", call. = FALSE)
  }
  assert_positive(sigma2, "sigma2")

  if (is.null(y)) {
    y <- x
    same_input <- TRUE
  } else {
    y <- as.matrix(y)
    storage.mode(y) <- "double"
    same_input <- FALSE
    if (!nrow(y) || ncol(y) != ncol(x) || anyNA(y) || any(!is.finite(y))) {
      stop("`y` must be finite and have the same columns as `x`.", call. = FALSE)
    }
  }

  x_norm <- rowSums(x^2)
  y_norm <- rowSums(y^2)
  distance2 <- outer(x_norm, y_norm, "+") - 2 * tcrossprod(x, y)
  distance2[distance2 < 0 & distance2 > -1e-10] <- 0
  if (any(distance2 < 0)) {
    distance2 <- pmax(distance2, 0)
  }
  kernel <- exp(-distance2 / (2 * sigma2))
  if (same_input) {
    kernel <- (kernel + t(kernel)) / 2
    diag(kernel) <- 1
  }
  kernel
}

#' Matérn radial basis function kernel
#'
#' @param x Numeric matrix with observations in rows.
#' @param y Optional numeric matrix. If omitted, computes the Gram matrix of
#'   `x`.
#' @param sigma2 Positive squared length scale.
#' @param smoothness Positive Matérn smoothness parameter. On
#'   `d`-dimensional Euclidean space, the associated native space is
#'   equivalent to a Sobolev space of order `smoothness + d / 2`.
#'
#' @return A numeric kernel matrix.
#' @export
matern_kernel <- function(x, y = NULL, sigma2 = 1, smoothness = 2) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!nrow(x) || !ncol(x) || anyNA(x) || any(!is.finite(x))) {
    stop("`x` must be a finite, nonempty numeric matrix.", call. = FALSE)
  }
  assert_positive(sigma2, "sigma2")
  assert_positive(smoothness, "smoothness")
  if (length(smoothness) != 1L || !is.finite(smoothness)) {
    stop("`smoothness` must be a positive finite scalar.", call. = FALSE)
  }

  if (is.null(y)) {
    y <- x
    same_input <- TRUE
  } else {
    y <- as.matrix(y)
    storage.mode(y) <- "double"
    same_input <- FALSE
    if (!nrow(y) || ncol(y) != ncol(x) || anyNA(y) || any(!is.finite(y))) {
      stop("`y` must be finite and have the same columns as `x`.", call. = FALSE)
    }
  }

  x_norm <- rowSums(x^2)
  y_norm <- rowSums(y^2)
  distance2 <- outer(x_norm, y_norm, "+") - 2 * tcrossprod(x, y)
  distance2 <- pmax(distance2, 0)
  scaled_distance <- sqrt(
    2 * smoothness * distance2 / sigma2
  )
  kernel <- matrix(1, nrow(x), nrow(y))
  positive <- scaled_distance > sqrt(.Machine$double.eps)
  if (any(positive)) {
    distance <- scaled_distance[positive]
    log_kernel <- (1 - smoothness) * log(2) -
      lgamma(smoothness) +
      smoothness * log(distance) +
      log(besselK(
        distance, nu = smoothness, expon.scaled = TRUE
      )) -
      distance
    kernel[positive] <- exp(log_kernel)
  }
  kernel[!is.finite(kernel)] <- 0
  kernel <- pmin(pmax(kernel, 0), 1)
  if (same_input) {
    kernel <- (kernel + t(kernel)) / 2
    diag(kernel) <- 1
  }
  kernel
}

kernel_matrix <- function(x, y = NULL, sigma2 = 1,
                          kernel_family = c(
                            "gaussian", "matern_sobolev"
                          ),
                          matern_smoothness = 2) {
  kernel_family <- match.arg(kernel_family)
  if (identical(kernel_family, "gaussian")) {
    return(gaussian_kernel(x, y, sigma2 = sigma2))
  }
  matern_kernel(
    x, y, sigma2 = sigma2, smoothness = matern_smoothness
  )
}

#' Weighted median-heuristic kernel variance
#'
#' Defines the base variance as one half of the weighted median of pairwise
#' squared Euclidean distances. Pair weights are products `w[i] * w[j]`. This
#' convention is explicit so that kernel code and manuscript notation cannot
#' silently disagree about distances versus squared distances.
#'
#' @param x Numeric matrix with observations in rows.
#' @param w Positive observation weights.
#'
#' @return A positive scalar kernel variance.
#' @export
median_bandwidth <- function(x, w = rep(1, nrow(x))) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) < 2L || !ncol(x) || anyNA(x) || any(!is.finite(x))) {
    stop("At least two finite observations are required for a bandwidth.", call. = FALSE)
  }
  assert_positive(w, "w")
  if (length(w) != nrow(x)) {
    stop("`w` must have one value per row of `x`.", call. = FALSE)
  }

  distance2 <- as.numeric(stats::dist(x))^2
  pair_weights <- outer(w, w, "*")[lower.tri(matrix(FALSE, nrow(x), nrow(x)))]
  sigma2 <- weighted_median(distance2, pair_weights) / 2

  if (!is.finite(sigma2) || sigma2 <= .Machine$double.eps) {
    positive <- distance2[distance2 > .Machine$double.eps]
    if (!length(positive)) {
      stop(
        "All observations are identical; a positive kernel bandwidth cannot be computed.",
        call. = FALSE
      )
    }
    sigma2 <- stats::median(positive) / 2
  }
  sigma2
}

#' A sample-size-dependent Nystrom rank rule
#'
#' Creates a rank rule of the form
#' `ceiling(multiplier * n^exponent)`, truncated to the requested lower and
#' upper bounds and to the number of rows available in the current training
#' fold. Using an exponent strictly between zero and one makes the rank grow
#' with the fold sample size while remaining sublinear.
#'
#' @param exponent Power applied to the number of rows in the current fold.
#' @param multiplier Positive multiplier for the power rule.
#' @param min_rank Smallest requested rank before truncation to the fold size.
#' @param max_rank Largest requested rank. Use `Inf` for no fixed upper bound.
#'
#' @return A function that maps a positive fold sample size to a Nystrom rank.
#' @export
pmtp_nystrom_rank <- function(exponent = 2 / 3, multiplier = 2,
                              min_rank = 30L, max_rank = Inf) {
  if (!is.numeric(exponent) || length(exponent) != 1L || is.na(exponent) ||
      exponent <= 0 || exponent >= 1) {
    stop("`exponent` must lie strictly between zero and one.", call. = FALSE)
  }
  assert_positive(multiplier, "multiplier")
  assert_positive(min_rank, "min_rank")
  assert_positive(max_rank, "max_rank", allow_inf = TRUE)
  if (length(multiplier) != 1L || length(min_rank) != 1L ||
      length(max_rank) != 1L) {
    stop("Rank-rule arguments must be scalar.", call. = FALSE)
  }
  if (min_rank > max_rank) {
    stop("`min_rank` cannot exceed `max_rank`.", call. = FALSE)
  }

  rule <- function(n) {
    assert_positive(n, "n")
    if (length(n) != 1L || n != as.integer(n)) {
      stop("`n` must be a positive integer.", call. = FALSE)
    }
    requested <- ceiling(multiplier * n^exponent)
    as.integer(min(n, max(min_rank, min(max_rank, requested))))
  }
  structure(
    rule,
    class = c("pmtp_nystrom_rank_rule", "function"),
    exponent = exponent,
    multiplier = multiplier,
    min_rank = min_rank,
    max_rank = max_rank
  )
}

resolve_nystrom_rank <- function(rank, n) {
  requested <- if (is.function(rank)) rank(n) else rank
  if (!is.numeric(requested) || length(requested) != 1L ||
      is.na(requested) || !is.finite(requested) || requested <= 0) {
    stop(
      "`nystrom_rank` must be a positive number or a function returning one.",
      call. = FALSE
    )
  }
  as.integer(min(n, ceiling(requested)))
}

select_nystrom_landmarks <- function(n, rank, weights, method, seed) {
  if (rank >= n) return(seq_len(n))
  probability <- NULL
  if (identical(method, "weighted")) {
    probability <- weights / sum(weights)
  }
  withr::with_seed(seed, sample.int(
    n, size = rank, replace = FALSE, prob = probability
  ))
}

fit_nystrom_map <- function(x, sigma2, weights, control, seed_offset = 0L,
                            kernel_family = "gaussian",
                            matern_smoothness = 2) {
  kernel_family <- match.arg(
    kernel_family, c("gaussian", "matern_sobolev")
  )
  x <- as.matrix(x)
  n <- nrow(x)
  rank <- resolve_nystrom_rank(control$nystrom_rank, n)
  landmark_indices <- select_nystrom_landmarks(
    n = n,
    rank = rank,
    weights = weights,
    method = control$nystrom_landmarks,
    seed = control$seed + seed_offset
  )
  landmarks <- x[landmark_indices, , drop = FALSE]
  landmark_kernel <- kernel_matrix(
    landmarks,
    sigma2 = sigma2,
    kernel_family = kernel_family,
    matern_smoothness = matern_smoothness
  )
  decomposition <- eigen(
    (landmark_kernel + t(landmark_kernel)) / 2,
    symmetric = TRUE
  )
  cutoff <- max(decomposition$values) *
    .Machine$double.eps * max(1, length(landmark_indices))
  keep <- which(decomposition$values > cutoff)
  if (!length(keep)) {
    stop("The Nystrom landmark kernel has zero numerical rank.", call. = FALSE)
  }
  transform <- sweep(
    decomposition$vectors[, keep, drop = FALSE],
    2L,
    sqrt(decomposition$values[keep]),
    "/"
  )
  training_features <- kernel_matrix(
    x,
    landmarks,
    sigma2 = sigma2,
    kernel_family = kernel_family,
    matern_smoothness = matern_smoothness
  ) %*% transform

  structure(list(
    landmarks = landmarks,
    transform = transform,
    training_features = training_features,
    sigma2 = sigma2,
    kernel_family = kernel_family,
    matern_smoothness = matern_smoothness,
    n_observations = n,
    requested_rank = rank,
    effective_rank = length(keep),
    landmark_indices = landmark_indices
  ), class = "pmtp_nystrom_map")
}

predict_nystrom_features <- function(object, new_arguments) {
  kernel_family <- if (is.null(object$kernel_family)) {
    "gaussian"
  } else {
    object$kernel_family
  }
  matern_smoothness <- if (is.null(object$matern_smoothness)) {
    2
  } else {
    object$matern_smoothness
  }
  kernel_matrix(
    new_arguments,
    object$landmarks,
    sigma2 = object$sigma2,
    kernel_family = kernel_family,
    matern_smoothness = matern_smoothness
  ) %*% object$transform
}
