tuning_grid_h <- function(control) {
  compact_grid(
    control$lambda_h,
    control$lambda_gp,
    control$bandwidth_h,
    control$bandwidth_gp
  )
}

tuning_grid_g <- function(control) {
  compact_grid(
    control$lambda_g,
    control$lambda_hp,
    control$bandwidth_g,
    control$bandwidth_hp
  )
}

fixed_outcome_tuning <- function(control) {
  selected <- tuning_grid_h(control)
  selected$mean_risk <- NA_real_
  selected$grid_index <- NA_integer_
  list(selected = selected, results = NULL, fold_id = NULL, fixed = TRUE)
}

fixed_treatment_tuning <- function(control) {
  selected <- tuning_grid_g(control)
  selected$mean_risk <- NA_real_
  selected$grid_index <- NA_integer_
  list(selected = selected, results = NULL, fold_id = NULL, fixed = TRUE)
}

summarize_cv_grid <- function(grid, risks) {
  finite_count <- rowSums(is.finite(risks))
  mean_risk <- vapply(seq_len(nrow(risks)), function(i) {
    if (finite_count[i] == 0L) Inf else mean(risks[i, is.finite(risks[i, ])])
  }, numeric(1))
  out <- cbind(
    grid,
    mean_risk = mean_risk,
    finite_folds = finite_count
  )
  for (j in seq_len(ncol(risks))) {
    out[[paste0("fold", j, "_risk")]] <- risks[, j]
  }
  out
}

select_cv_row <- function(results, grid, bridge_name) {
  if (!any(is.finite(results$mean_risk))) {
    stop(
      "Every ", bridge_name, " bridge tuning configuration failed. ",
      "Inspect variable scaling, policy support, and tuning ranges.",
      call. = FALSE
    )
  }
  selected_index <- which.min(results$mean_risk)
  selected <- grid[selected_index, , drop = FALSE]
  selected$mean_risk <- results$mean_risk[selected_index]
  selected$grid_index <- selected_index
  selected
}

warn_grid_boundary <- function(selected, grid, bridge_name) {
  fields <- c(
    "outer_lambda_scale", "inner_lambda_scale",
    "outer_bandwidth_scale", "inner_bandwidth_scale"
  )
  boundary <- vapply(fields, function(field) {
    is_grid_boundary(selected[[field]], grid[[field]])
  }, logical(1))
  if (any(boundary)) {
    warning(
      bridge_name, " tuning selected a grid boundary for: ",
      paste(fields[boundary], collapse = ", "),
      ". Consider expanding the grid.",
      call. = FALSE
    )
  }
}

tune_outcome_bridge <- function(dat, indices, control, seed_offset = 0L) {
  grid <- tuning_grid_h(control)
  fold_id <- make_folds(
    dat$a[indices], dat$y[indices], control$inner_folds,
    control$seed + seed_offset
  )
  risks <- matrix(Inf, nrow(grid), control$inner_folds)

  for (fold in seq_len(control$inner_folds)) {
    training <- indices[fold_id != fold]
    validation <- indices[fold_id == fold]
    prepared <- prepare_fold_data(dat, training, validation)
    # Penalty rates are governed by the number of observed phase-two rows,
    # not by the Horvitz-Thompson estimate of the phase-one population size.
    # The weighted empirical objectives themselves remain normalized by the
    # sum of weights inside the bridge-fitting and risk functions.
    n_training <- length(training)
    n_validation <- length(validation)
    base_h <- median_bandwidth(prepared$train$h, dat$weight[training])
    base_gp <- median_bandwidth(prepared$train$g, dat$weight[training])
    risk_base <- median_bandwidth(
      prepared$validation$g,
      dat$weight[validation]
    )
    risk_lambda <- control$risk_penalty * log(n_validation) / n_validation

    for (candidate in seq_len(nrow(grid))) {
      risks[candidate, fold] <- tryCatch({
        fit <- fit_outcome_bridge(
          h = prepared$train$h,
          gp = prepared$train$g,
          y = dat$y[training],
          weights = dat$weight[training],
          lambda_h = actual_outer_lambda(
            grid$outer_lambda_scale[candidate], n_training
          ),
          lambda_gp = actual_inner_lambda(
            grid$inner_lambda_scale[candidate], n_training
          ),
          sigma2_h = base_h * grid$outer_bandwidth_scale[candidate],
          sigma2_gp = base_gp * grid$inner_bandwidth_scale[candidate],
          max_norm = control$max_norm_h,
          control = control
        )
        prediction <- predict_outcome_bridge(fit, prepared$validation$h)
        outcome_validation_risk(
          residual = dat$y[validation] - prediction,
          adversary_arguments = prepared$validation$g,
          weights = dat$weight[validation],
          sigma2 = risk_base * control$risk_bandwidth,
          lambda = risk_lambda,
          control = control
        )
      }, error = function(e) Inf)
    }
  }

  results <- summarize_cv_grid(grid, risks)
  selected <- select_cv_row(results, grid, "Outcome")
  warn_grid_boundary(selected, grid, "Outcome bridge")
  list(selected = selected, results = results, fold_id = fold_id)
}

tune_treatment_bridge <- function(dat, indices, policy_index, control,
                                  seed_offset = 0L) {
  grid <- tuning_grid_g(control)
  fold_id <- make_folds(
    dat$a[indices], dat$y[indices], control$inner_folds,
    control$seed + seed_offset
  )
  risks <- matrix(Inf, nrow(grid), control$inner_folds)

  for (fold in seq_len(control$inner_folds)) {
    training <- indices[fold_id != fold]
    validation <- indices[fold_id == fold]
    prepared <- prepare_fold_data(dat, training, validation, policy_index)
    n_training <- length(training)
    n_validation <- length(validation)
    base_g <- median_bandwidth(prepared$train$g, dat$weight[training])
    base_hp <- median_bandwidth(prepared$train$h, dat$weight[training])
    risk_base <- median_bandwidth(
      prepared$validation$h,
      dat$weight[validation]
    )
    risk_lambda <- control$risk_penalty * log(n_validation) / n_validation

    for (candidate in seq_len(nrow(grid))) {
      risks[candidate, fold] <- tryCatch({
        fit <- fit_treatment_bridge(
          g = prepared$train$g,
          hp = prepared$train$h,
          hp_q = prepared$train$hq,
          weights = dat$weight[training],
          target = dat$target[training],
          policy_support = dat$policy_support[[policy_index]][training],
          lambda_g = actual_outer_lambda(
            grid$outer_lambda_scale[candidate], n_training
          ),
          lambda_hp = actual_inner_lambda(
            grid$inner_lambda_scale[candidate], n_training
          ),
          sigma2_g = base_g * grid$outer_bandwidth_scale[candidate],
          sigma2_hp = base_hp * grid$inner_bandwidth_scale[candidate],
          max_norm = control$max_norm_g,
          control = control
        )
        prediction <- predict_treatment_bridge(fit, prepared$validation$g)
        treatment_validation_risk(
          g_value = prediction,
          adversary_arguments = prepared$validation$h,
          adversary_policy_arguments = prepared$validation$hq,
          weights = dat$weight[validation],
          target = dat$target[validation],
          policy_support = dat$policy_support[[policy_index]][validation],
          sigma2 = risk_base * control$risk_bandwidth,
          lambda = risk_lambda,
          control = control
        )
      }, error = function(e) Inf)
    }
  }

  results <- summarize_cv_grid(grid, risks)
  selected <- select_cv_row(results, grid, "Treatment")
  warn_grid_boundary(selected, grid, "Treatment bridge")
  list(selected = selected, results = results, fold_id = fold_id)
}

fit_selected_outcome <- function(dat, training, validation, tuning, control) {
  prepared <- prepare_fold_data(dat, training, validation)
  n_training <- length(training)
  base_h <- median_bandwidth(prepared$train$h, dat$weight[training])
  base_gp <- median_bandwidth(prepared$train$g, dat$weight[training])
  fit <- fit_outcome_bridge(
    h = prepared$train$h,
    gp = prepared$train$g,
    y = dat$y[training],
    weights = dat$weight[training],
    lambda_h = actual_outer_lambda(tuning$outer_lambda_scale, n_training),
    lambda_gp = actual_inner_lambda(tuning$inner_lambda_scale, n_training),
    sigma2_h = base_h * tuning$outer_bandwidth_scale,
    sigma2_gp = base_gp * tuning$inner_bandwidth_scale,
    max_norm = control$max_norm_h,
    control = control
  )
  list(
    fit = fit,
    h0 = predict_outcome_bridge(fit, prepared$validation$h),
    prepared = prepared,
    base_bandwidth = c(h = base_h, gp = base_gp)
  )
}

fit_selected_treatment <- function(dat, training, validation, policy_index,
                                   tuning, prepared, control) {
  n_training <- length(training)
  base_g <- median_bandwidth(prepared$train$g, dat$weight[training])
  base_hp <- median_bandwidth(prepared$train$h, dat$weight[training])
  fit <- fit_treatment_bridge(
    g = prepared$train$g,
    hp = prepared$train$h,
    hp_q = prepared$train$hq,
    weights = dat$weight[training],
    target = dat$target[training],
    policy_support = dat$policy_support[[policy_index]][training],
    lambda_g = actual_outer_lambda(tuning$outer_lambda_scale, n_training),
    lambda_hp = actual_inner_lambda(tuning$inner_lambda_scale, n_training),
    sigma2_g = base_g * tuning$outer_bandwidth_scale,
    sigma2_hp = base_hp * tuning$inner_bandwidth_scale,
    max_norm = control$max_norm_g,
    control = control
  )
  list(
    fit = fit,
    g0 = predict_treatment_bridge(fit, prepared$validation$g),
    base_bandwidth = c(g = base_g, hp = base_hp)
  )
}
