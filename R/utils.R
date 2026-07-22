`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

assert_positive <- function(x, name, allow_inf = FALSE) {
  ok <- is.numeric(x) && length(x) >= 1L && !anyNA(x) &&
    all(x > 0) && (allow_inf || all(is.finite(x)))
  if (!ok) {
    stop("`", name, "` must contain positive numeric values.", call. = FALSE)
  }
  invisible(x)
}

assert_columns <- function(data, columns, argument) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      "Columns supplied through `", argument, "` are missing: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(columns)
}

as_numeric_matrix <- function(data, columns, argument) {
  assert_columns(data, columns, argument)
  bad <- columns[!vapply(data[columns], is.numeric, logical(1))]
  if (length(bad)) {
    stop(
      "All variables supplied through `", argument,
      "` must be numeric. Non-numeric variables: ",
      paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }
  out <- as.matrix(data[columns])
  storage.mode(out) <- "double"
  out
}

safe_solve <- function(a, b, jitter = 1e-8, max_tries = 8L,
                       symmetric = TRUE) {
  a <- as.matrix(a)
  if (nrow(a) != ncol(a)) {
    stop("The matrix supplied to `safe_solve()` must be square.", call. = FALSE)
  }
  if (symmetric) {
    a <- (a + t(a)) / 2
  }

  b_was_vector <- is.null(dim(b))
  b <- as.matrix(b)
  n <- nrow(a)
  if (nrow(b) != n) {
    stop("The right-hand side has incompatible dimensions.", call. = FALSE)
  }

  matrix_scale <- mean(abs(diag(a)))
  if (!is.finite(matrix_scale) || matrix_scale <= 0) {
    matrix_scale <- max(abs(a))
  }
  if (!is.finite(matrix_scale) || matrix_scale <= 0) {
    matrix_scale <- 1
  }

  for (attempt in 0:max_tries) {
    ridge <- if (attempt == 0L) 0 else jitter * 10^(attempt - 1L)
    a_try <- a
    diag(a_try) <- diag(a_try) + ridge * matrix_scale
    factor <- tryCatch(chol(a_try), error = function(e) NULL)
    if (!is.null(factor)) {
      solution <- backsolve(factor, forwardsolve(t(factor), b))
      if (all(is.finite(solution))) {
        if (b_was_vector) return(drop(solution))
        return(solution)
      }
    }
  }

  eig <- eigen(a, symmetric = TRUE)
  cutoff <- max(abs(eig$values)) * max(jitter, .Machine$double.eps^0.75)
  inverse_values <- ifelse(eig$values > cutoff, 1 / eig$values, 0)
  solution <- eig$vectors %*%
    (inverse_values * crossprod(eig$vectors, b))

  if (!all(is.finite(solution))) {
    stop("A regularized linear system could not be solved.", call. = FALSE)
  }
  if (b_was_vector) drop(solution) else solution
}

weighted_mean <- function(x, w) {
  sum(w * x) / sum(w)
}

weighted_variance <- function(x, w) {
  center <- weighted_mean(x, w)
  sum(w * (x - center)^2) / sum(w)
}

weighted_median <- function(x, w) {
  if (!length(x) || length(x) != length(w)) {
    stop("`x` and `w` must have the same positive length.", call. = FALSE)
  }
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x)) {
    stop("No finite values with positive weight were available.", call. = FALSE)
  }
  order_x <- order(x)
  x <- x[order_x]
  w <- w[order_x]
  cumulative <- cumsum(w)
  half <- sum(w) / 2
  index <- which(cumulative >= half)[1L]
  if (index < length(x) && isTRUE(all.equal(cumulative[index], half))) {
    return(mean(x[c(index, index + 1L)]))
  }
  x[index]
}

make_folds <- function(a, y, k, seed) {
  n <- length(a)
  if (k < 2L || k > n) {
    stop("The number of folds must be between 2 and the sample size.", call. = FALSE)
  }

  breaks <- unique(stats::quantile(
    a,
    probs = seq(0, 1, length.out = min(6L, n + 1L)),
    na.rm = TRUE,
    type = 8
  ))
  a_group <- if (length(breaks) >= 3L) {
    cut(a, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  } else {
    rep.int(1L, n)
  }
  y_group <- if (length(unique(y)) <= 5L) match(y, unique(y)) else 1L
  strata <- interaction(a_group, y_group, drop = TRUE)

  withr::with_seed(seed, {
    fold_id <- integer(n)
    for (indices in split(seq_len(n), strata)) {
      indices <- sample(indices)
      fold_id[indices] <- rep(seq_len(k), length.out = length(indices))
    }
    if (any(tabulate(fold_id, nbins = k) == 0L)) {
      fold_id <- sample(rep(seq_len(k), length.out = n))
    }
    fold_id
  })
}

compact_grid <- function(outer_lambda, inner_lambda,
                         outer_bandwidth, inner_bandwidth) {
  grid <- expand.grid(
    outer_lambda_scale = outer_lambda,
    inner_lambda_scale = inner_lambda,
    outer_bandwidth_scale = outer_bandwidth,
    inner_bandwidth_scale = inner_bandwidth,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(grid) <- NULL
  grid
}

critical_radius_squared <- function(population_size, dimension = 1L,
                                    rule = c(
                                      "legacy_d1", "gaussian_dimension"
                                    )) {
  rule <- match.arg(rule)
  if (length(population_size) != 1L || !is.numeric(population_size) ||
      is.na(population_size) || !is.finite(population_size) ||
      population_size <= 1) {
    stop("`population_size` must be a finite number greater than one.",
         call. = FALSE)
  }
  if (length(dimension) != 1L || !is.numeric(dimension) ||
      is.na(dimension) || !is.finite(dimension) || dimension < 1L ||
      dimension != as.integer(dimension)) {
    stop("`dimension` must be a positive integer.", call. = FALSE)
  }
  exponent <- if (identical(rule, "legacy_d1")) 1L else as.integer(dimension)
  log(population_size)^exponent / population_size
}

actual_outer_lambda <- function(scale, population_size, dimension = 1L,
                                rule = "legacy_d1") {
  scale * sqrt(critical_radius_squared(
    population_size, dimension = dimension, rule = rule
  ))
}

actual_inner_lambda <- function(scale, population_size, dimension = 1L,
                                rule = "legacy_d1") {
  scale * critical_radius_squared(
    population_size, dimension = dimension, rule = rule
  )
}

actual_risk_lambda <- function(scale, population_size, dimension = 1L,
                               rule = "legacy_d1") {
  actual_inner_lambda(
    scale, population_size, dimension = dimension, rule = rule
  )
}

is_grid_boundary <- function(value, candidates) {
  candidates <- sort(unique(candidates))
  length(candidates) > 1L && value %in% range(candidates)
}
