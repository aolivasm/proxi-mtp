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

use_nystrom_feature_cache <- function(control) {
  identical(control$kernel_approximation, "nystrom") &&
    isTRUE(control$cache_kernel_features)
}

cache_nystrom_maps <- function(arguments, sigma2, weights, control,
                               seed_offset) {
  lapply(sigma2, function(value) {
    fit_nystrom_map(
      arguments, value, weights, control, seed_offset = seed_offset
    )
  })
}

nystrom_penalty_path_groups <- function(grid) {
  split(
    seq_len(nrow(grid)),
    interaction(
      grid$inner_lambda_scale,
      grid$outer_bandwidth_scale,
      grid$inner_bandwidth_scale,
      drop = TRUE,
      lex.order = TRUE
    )
  )
}

evaluate_outcome_nystrom_grid <- function(
    prepared, training_y, validation_y, training_weights, validation_weights, grid,
    n_training, base_h, base_gp, risk_base, risk_lambda,
    feature_cache, control) {
  risks <- rep(Inf, nrow(grid))
  groups <- nystrom_penalty_path_groups(grid)
  for (indices in groups) {
    first <- indices[1L]
    outer_index <- match(
      grid$outer_bandwidth_scale[first], feature_cache$outer_scales
    )
    inner_index <- match(
      grid$inner_bandwidth_scale[first], feature_cache$inner_scales
    )
    system <- tryCatch(
      prepare_outcome_bridge_nystrom_system(
        feature_cache$outer[[outer_index]],
        feature_cache$inner[[inner_index]],
        training_y,
        training_weights,
        actual_inner_lambda(
          grid$inner_lambda_scale[first], n_training
        ),
        control
      ),
      error = function(e) NULL
    )
    if (is.null(system)) next
    for (candidate in indices) {
      risks[candidate] <- tryCatch({
        fit <- finish_outcome_bridge_nystrom_system(
          system = system,
          lambda_h = actual_outer_lambda(
            grid$outer_lambda_scale[candidate], n_training
          ),
          sigma2_h = base_h * grid$outer_bandwidth_scale[candidate],
          sigma2_gp = base_gp * grid$inner_bandwidth_scale[candidate],
          max_norm = control$max_norm_h,
          control = control
        )
        prediction <- predict_outcome_bridge(fit, prepared$validation$h)
        outcome_validation_risk(
          residual = validation_y - prediction,
          adversary_arguments = prepared$validation$g,
          weights = validation_weights,
          sigma2 = risk_base * control$risk_bandwidth,
          lambda = risk_lambda,
          control = control,
          feature_map = feature_cache$risk
        )
      }, error = function(e) Inf)
    }
  }
  risks
}

evaluate_treatment_nystrom_grid <- function(
    prepared, training_weights, validation_weights, training_target,
    validation_target, training_support, validation_support, grid,
    n_training, base_g, base_hp, risk_base, risk_lambda,
    feature_cache, control) {
  risks <- rep(Inf, nrow(grid))
  groups <- nystrom_penalty_path_groups(grid)
  for (indices in groups) {
    first <- indices[1L]
    outer_index <- match(
      grid$outer_bandwidth_scale[first], feature_cache$outer_scales
    )
    inner_index <- match(
      grid$inner_bandwidth_scale[first], feature_cache$inner_scales
    )
    system <- tryCatch(
      prepare_treatment_bridge_nystrom_system(
        feature_cache$outer[[outer_index]],
        feature_cache$inner[[inner_index]],
        feature_cache$policy_inner_features[[inner_index]],
        training_weights,
        training_target,
        training_support,
        actual_inner_lambda(
          grid$inner_lambda_scale[first], n_training
        ),
        control
      ),
      error = function(e) NULL
    )
    if (is.null(system)) next
    for (candidate in indices) {
      risks[candidate] <- tryCatch({
        fit <- finish_treatment_bridge_nystrom_system(
          system = system,
          lambda_g = actual_outer_lambda(
            grid$outer_lambda_scale[candidate], n_training
          ),
          sigma2_g = base_g * grid$outer_bandwidth_scale[candidate],
          sigma2_hp = base_hp * grid$inner_bandwidth_scale[candidate],
          max_norm = control$max_norm_g,
          control = control
        )
        prediction <- predict_treatment_bridge(
          fit, prepared$validation$g
        )
        treatment_validation_risk(
          g_value = prediction,
          adversary_arguments = prepared$validation$h,
          adversary_policy_arguments = prepared$validation$hq,
          weights = validation_weights,
          target = validation_target,
          policy_support = validation_support,
          sigma2 = risk_base * control$risk_bandwidth,
          lambda = risk_lambda,
          control = control,
          feature_map = feature_cache$risk,
          policy_features = feature_cache$risk_policy_features
        )
      }, error = function(e) Inf)
    }
  }
  risks
}

fixed_outcome_tuning <- function(control) {
  selected <- tuning_grid_h(control)
  selected$mean_risk <- NA_real_
  selected$grid_index <- NA_integer_
  list(
    selected = selected, results = NULL, fold_id = NULL,
    risk_folds = NULL, fixed = TRUE
  )
}

fixed_treatment_tuning <- function(control) {
  selected <- tuning_grid_g(control)
  selected$mean_risk <- NA_real_
  selected$grid_index <- NA_integer_
  list(
    selected = selected, results = NULL, fold_id = NULL,
    risk_folds = NULL, fixed = TRUE
  )
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

make_repeated_inner_folds <- function(dat, indices, control, seed_offset) {
  fold_id <- vapply(seq_len(control$inner_repeats), function(inner_repeat) {
    make_folds(
      dat$a[indices], dat$y[indices], control$inner_folds,
      control$seed + seed_offset + 100000L * (inner_repeat - 1L)
    )
  }, integer(length(indices)))
  if (is.null(dim(fold_id))) {
    fold_id <- matrix(fold_id, ncol = 1L)
  }
  colnames(fold_id) <- paste0("repeat", seq_len(control$inner_repeats))
  risk_folds <- expand.grid(
    inner_fold = seq_len(control$inner_folds),
    inner_repeat = seq_len(control$inner_repeats),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  risk_folds$risk_index <- seq_len(nrow(risk_folds))
  list(fold_id = fold_id, risk_folds = risk_folds)
}

select_cv_row <- function(results, grid, bridge_name,
                          selection_rule = "minimum") {
  selection_rule <- match.arg(
    selection_rule,
    c("minimum", "one_se_regularized", "one_se_interior")
  )
  if (!any(is.finite(results$mean_risk))) {
    stop(
      "Every ", bridge_name, " bridge tuning configuration failed. ",
      "Inspect variable scaling, policy support, and tuning ranges.",
      call. = FALSE
    )
  }
  minimum_index <- which.min(results$mean_risk)
  fold_columns <- grep("^fold[0-9]+_risk$", names(results), value = TRUE)
  minimum_fold_risks <- unlist(
    results[minimum_index, fold_columns, drop = FALSE],
    use.names = FALSE
  )
  minimum_fold_risks <- minimum_fold_risks[is.finite(minimum_fold_risks)]
  minimum_risk_se <- if (length(minimum_fold_risks) > 1L) {
    stats::sd(minimum_fold_risks) / sqrt(length(minimum_fold_risks))
  } else {
    0
  }
  one_se_threshold <- results$mean_risk[minimum_index] + minimum_risk_se
  one_se_indices <- which(
    is.finite(results$mean_risk) & results$mean_risk <= one_se_threshold
  )
  selected_index <- minimum_index
  if (identical(selection_rule, "one_se_regularized")) {
    largest_penalty <- max(grid$outer_lambda_scale[one_se_indices])
    eligible <- one_se_indices[
      grid$outer_lambda_scale[one_se_indices] == largest_penalty
    ]
    largest_bandwidth <- max(grid$outer_bandwidth_scale[eligible])
    eligible <- eligible[
      grid$outer_bandwidth_scale[eligible] == largest_bandwidth
    ]
    selected_index <- eligible[which.min(results$mean_risk[eligible])]
  } else if (identical(selection_rule, "one_se_interior")) {
    boundary_fields <- c(
      "outer_lambda_scale", "inner_lambda_scale",
      "outer_bandwidth_scale", "inner_bandwidth_scale"
    )
    boundary_count <- rowSums(vapply(boundary_fields, function(field) {
      candidates <- sort(unique(grid[[field]]))
      length(candidates) > 1L & grid[[field]] %in% range(candidates)
    }, logical(nrow(grid))))
    fewest_boundaries <- min(boundary_count[one_se_indices])
    eligible <- one_se_indices[
      boundary_count[one_se_indices] == fewest_boundaries
    ]
    selected_index <- eligible[which.min(results$mean_risk[eligible])]
  }
  selected <- grid[selected_index, , drop = FALSE]
  selected$mean_risk <- results$mean_risk[selected_index]
  selected$grid_index <- selected_index
  selected$selection_rule <- selection_rule
  selected$minimum_mean_risk <- results$mean_risk[minimum_index]
  selected$minimum_grid_index <- minimum_index
  selected$minimum_risk_se <- minimum_risk_se
  selected$one_se_threshold <- one_se_threshold
  selected$one_se_candidate_count <- length(one_se_indices)
  selected$selected_risk_gap <-
    results$mean_risk[selected_index] - results$mean_risk[minimum_index]
  selected
}

warn_grid_boundary <- function(selected, results, grid, bridge_name) {
  fields <- c(
    "outer_lambda_scale", "inner_lambda_scale",
    "outer_bandwidth_scale", "inner_bandwidth_scale"
  )
  minimum <- grid[selected$minimum_grid_index, , drop = FALSE]
  boundary <- vapply(fields, function(field) {
    is_grid_boundary(minimum[[field]], grid[[field]])
  }, logical(1))
  within_one_se <- is.finite(results$mean_risk) &
    results$mean_risk <= selected$one_se_threshold
  has_interior_alternative <- vapply(fields, function(field) {
    candidates <- sort(unique(grid[[field]]))
    if (length(candidates) <= 1L) return(TRUE)
    any(within_one_se & !grid[[field]] %in% range(candidates))
  }, logical(1))
  unsupported_boundary <- boundary & !has_interior_alternative
  if (any(unsupported_boundary)) {
    warning(
      bridge_name,
      " selected a grid boundary with no interior candidate within one ",
      "cross-validation standard error for: ",
      paste(fields[unsupported_boundary], collapse = ", "),
      ". Consider expanding that grid dimension.",
      call. = FALSE
    )
  }
  invisible(list(
    boundary = stats::setNames(boundary, fields),
    has_interior_alternative = stats::setNames(
      has_interior_alternative, fields
    ),
    unsupported_boundary = stats::setNames(unsupported_boundary, fields)
  ))
}

tune_outcome_bridge <- function(dat, indices, control, seed_offset = 0L) {
  grid <- tuning_grid_h(control)
  repeated_folds <- make_repeated_inner_folds(
    dat, indices, control, seed_offset
  )
  risks <- matrix(Inf, nrow(grid), nrow(repeated_folds$risk_folds))

  for (risk_index in seq_len(nrow(repeated_folds$risk_folds))) {
    fold <- repeated_folds$risk_folds$inner_fold[risk_index]
    inner_repeat <- repeated_folds$risk_folds$inner_repeat[risk_index]
    fold_id <- repeated_folds$fold_id[, inner_repeat]
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
    feature_cache <- NULL
    if (use_nystrom_feature_cache(control)) {
      outer_scales <- unique(grid$outer_bandwidth_scale)
      inner_scales <- unique(grid$inner_bandwidth_scale)
      feature_cache <- list(
        outer_scales = outer_scales,
        inner_scales = inner_scales,
        outer = cache_nystrom_maps(
          prepared$train$h,
          base_h * outer_scales,
          dat$weight[training],
          control,
          seed_offset = 101L
        ),
        inner = cache_nystrom_maps(
          prepared$train$g,
          base_gp * inner_scales,
          dat$weight[training],
          control,
          seed_offset = 102L
        ),
        risk = fit_nystrom_map(
          prepared$validation$g,
          risk_base * control$risk_bandwidth,
          dat$weight[validation],
          control,
          seed_offset = 301L
        )
      )
    }
    if (!is.null(feature_cache)) {
      risks[, risk_index] <- evaluate_outcome_nystrom_grid(
        prepared = prepared,
        training_y = dat$y[training],
        validation_y = dat$y[validation],
        training_weights = dat$weight[training],
        validation_weights = dat$weight[validation],
        grid = grid,
        n_training = n_training,
        base_h = base_h,
        base_gp = base_gp,
        risk_base = risk_base,
        risk_lambda = risk_lambda,
        feature_cache = feature_cache,
        control = control
      )
      next
    }

    for (candidate in seq_len(nrow(grid))) {
      risks[candidate, risk_index] <- tryCatch({
        candidate_maps <- if (is.null(feature_cache)) NULL else list(
          outer = feature_cache$outer[[match(
            grid$outer_bandwidth_scale[candidate],
            feature_cache$outer_scales
          )]],
          inner = feature_cache$inner[[match(
            grid$inner_bandwidth_scale[candidate],
            feature_cache$inner_scales
          )]]
        )
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
          control = control,
          feature_maps = candidate_maps
        )
        prediction <- predict_outcome_bridge(fit, prepared$validation$h)
        outcome_validation_risk(
          residual = dat$y[validation] - prediction,
          adversary_arguments = prepared$validation$g,
          weights = dat$weight[validation],
          sigma2 = risk_base * control$risk_bandwidth,
          lambda = risk_lambda,
          control = control,
          feature_map = feature_cache$risk %||% NULL
        )
      }, error = function(e) Inf)
    }
  }

  results <- summarize_cv_grid(grid, risks)
  selected <- select_cv_row(
    results, grid, "Outcome", control$selection_rule
  )
  warn_grid_boundary(
    selected,
    results,
    grid,
    "Outcome bridge risk minimum"
  )
  stored_fold_id <- if (control$inner_repeats == 1L) {
    repeated_folds$fold_id[, 1L]
  } else {
    repeated_folds$fold_id
  }
  list(
    selected = selected,
    results = results,
    fold_id = stored_fold_id,
    risk_folds = repeated_folds$risk_folds
  )
}

tune_treatment_bridge <- function(dat, indices, policy_index, control,
                                  seed_offset = 0L) {
  grid <- tuning_grid_g(control)
  repeated_folds <- make_repeated_inner_folds(
    dat, indices, control, seed_offset
  )
  risks <- matrix(Inf, nrow(grid), nrow(repeated_folds$risk_folds))

  for (risk_index in seq_len(nrow(repeated_folds$risk_folds))) {
    fold <- repeated_folds$risk_folds$inner_fold[risk_index]
    inner_repeat <- repeated_folds$risk_folds$inner_repeat[risk_index]
    fold_id <- repeated_folds$fold_id[, inner_repeat]
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
    feature_cache <- NULL
    if (use_nystrom_feature_cache(control)) {
      outer_scales <- unique(grid$outer_bandwidth_scale)
      inner_scales <- unique(grid$inner_bandwidth_scale)
      inner_maps <- cache_nystrom_maps(
        prepared$train$h,
        base_hp * inner_scales,
        dat$weight[training],
        control,
        seed_offset = 202L
      )
      risk_map <- fit_nystrom_map(
        prepared$validation$h,
        risk_base * control$risk_bandwidth,
        dat$weight[validation],
        control,
        seed_offset = 302L
      )
      feature_cache <- list(
        outer_scales = outer_scales,
        inner_scales = inner_scales,
        outer = cache_nystrom_maps(
          prepared$train$g,
          base_g * outer_scales,
          dat$weight[training],
          control,
          seed_offset = 201L
        ),
        inner = inner_maps,
        policy_inner_features = lapply(inner_maps, function(map) {
          predict_nystrom_features(map, prepared$train$hq)
        }),
        risk = risk_map,
        risk_policy_features = predict_nystrom_features(
          risk_map, prepared$validation$hq
        )
      )
    }
    if (!is.null(feature_cache)) {
      risks[, risk_index] <- evaluate_treatment_nystrom_grid(
        prepared = prepared,
        training_weights = dat$weight[training],
        validation_weights = dat$weight[validation],
        training_target = dat$target[training],
        validation_target = dat$target[validation],
        training_support =
          dat$policy_support[[policy_index]][training],
        validation_support =
          dat$policy_support[[policy_index]][validation],
        grid = grid,
        n_training = n_training,
        base_g = base_g,
        base_hp = base_hp,
        risk_base = risk_base,
        risk_lambda = risk_lambda,
        feature_cache = feature_cache,
        control = control
      )
      next
    }

    for (candidate in seq_len(nrow(grid))) {
      risks[candidate, risk_index] <- tryCatch({
        candidate_maps <- if (is.null(feature_cache)) NULL else {
          inner_index <- match(
            grid$inner_bandwidth_scale[candidate],
            feature_cache$inner_scales
          )
          list(
            outer = feature_cache$outer[[match(
              grid$outer_bandwidth_scale[candidate],
              feature_cache$outer_scales
            )]],
            inner = feature_cache$inner[[inner_index]],
            policy_inner_features =
              feature_cache$policy_inner_features[[inner_index]]
          )
        }
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
          control = control,
          feature_maps = candidate_maps
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
          control = control,
          feature_map = feature_cache$risk %||% NULL,
          policy_features = feature_cache$risk_policy_features %||% NULL
        )
      }, error = function(e) Inf)
    }
  }

  results <- summarize_cv_grid(grid, risks)
  selected <- select_cv_row(
    results, grid, "Treatment", control$selection_rule
  )
  warn_grid_boundary(
    selected,
    results,
    grid,
    "Treatment bridge risk minimum"
  )
  stored_fold_id <- if (control$inner_repeats == 1L) {
    repeated_folds$fold_id[, 1L]
  } else {
    repeated_folds$fold_id
  }
  list(
    selected = selected,
    results = results,
    fold_id = stored_fold_id,
    risk_folds = repeated_folds$risk_folds
  )
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
