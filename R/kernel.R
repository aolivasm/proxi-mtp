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
