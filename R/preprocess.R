fit_weighted_scaler <- function(x, w) {
  x <- as.matrix(x)
  centers <- vapply(seq_len(ncol(x)), function(j) weighted_mean(x[, j], w), numeric(1))
  scales <- vapply(seq_len(ncol(x)), function(j) {
    sqrt(weighted_variance(x[, j], w))
  }, numeric(1))
  bad <- !is.finite(scales) | scales <= sqrt(.Machine$double.eps)
  if (any(bad)) {
    stop(
      "Variables with zero or near-zero weighted variance: ",
      paste(colnames(x)[bad], collapse = ", "), ".",
      call. = FALSE
    )
  }
  list(center = centers, scale = scales, names = colnames(x))
}

apply_weighted_scaler <- function(x, scaler) {
  x <- as.matrix(x)
  if (!identical(colnames(x), scaler$names)) {
    stop("Preprocessing columns do not match their training columns.", call. = FALSE)
  }
  sweep(sweep(x, 2L, scaler$center, "-"), 2L, scaler$scale, "/")
}

make_core_matrix <- function(dat, indices) {
  out <- cbind(
    A = dat$a[indices],
    dat$l[indices, , drop = FALSE],
    dat$z[indices, , drop = FALSE],
    dat$w_proxy[indices, , drop = FALSE]
  )
  colnames(out) <- dat$core_names
  out
}

make_components <- function(dat, indices, scaled_core, scaler, policy_index = NULL) {
  h <- scaled_core[, c(dat$index$a, dat$index$l, dat$index$w), drop = FALSE]
  g <- scaled_core[, c(dat$index$a, dat$index$l, dat$index$z), drop = FALSE]
  out <- list(h = h, g = g)

  if (!is.null(policy_index)) {
    q_scaled <- (dat$q[[policy_index]][indices] - scaler$center[dat$index$a]) /
      scaler$scale[dat$index$a]
    out$hq <- cbind(
      q_scaled,
      scaled_core[, dat$index$l, drop = FALSE],
      scaled_core[, dat$index$w, drop = FALSE]
    )
    colnames(out$hq) <- colnames(h)
  }
  out
}

prepare_fold_data <- function(dat, train, validation, policy_index = NULL) {
  train_core <- make_core_matrix(dat, train)
  validation_core <- make_core_matrix(dat, validation)
  scaler <- fit_weighted_scaler(train_core, dat$weight[train])
  train_scaled <- apply_weighted_scaler(train_core, scaler)
  validation_scaled <- apply_weighted_scaler(validation_core, scaler)
  list(
    train = make_components(dat, train, train_scaled, scaler, policy_index),
    validation = make_components(
      dat, validation, validation_scaled, scaler, policy_index
    ),
    scaler = scaler
  )
}
