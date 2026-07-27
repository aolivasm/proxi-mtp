rkhs_norm <- function(coef, kernel) {
  value <- drop(crossprod(coef, kernel %*% coef))
  sqrt(max(value, 0))
}

control_kernel_family <- function(control) {
  control$kernel_family %||% "gaussian"
}

control_matern_smoothness <- function(control) {
  control$matern_smoothness %||% 2
}

controlled_kernel_matrix <- function(x, y = NULL, sigma2, control) {
  kernel_matrix(
    x, y,
    sigma2 = sigma2,
    kernel_family = control_kernel_family(control),
    matern_smoothness = control_matern_smoothness(control)
  )
}

solve_with_norm_bound <- function(a, b, kernel, bound, control) {
  solve_at <- function(multiplier) {
    safe_solve(
      a + 4 * multiplier * kernel,
      b,
      jitter = control$jitter,
      max_tries = control$max_solve_tries
    )
  }
  coef <- solve_at(0)
  norm <- rkhs_norm(coef, kernel)
  if (!is.finite(bound) || norm <= bound) {
    return(list(coef = coef, norm = norm, constraint = 0))
  }

  objective <- function(multiplier) {
    rkhs_norm(solve_at(multiplier), kernel) - bound
  }
  upper <- 1
  while (objective(upper) > 0 && upper < 1e12) {
    upper <- upper * 10
  }
  if (upper >= 1e12 && objective(upper) > 0) {
    stop("The RKHS norm constraint could not be bracketed.", call. = FALSE)
  }
  multiplier <- stats::uniroot(objective, c(0, upper), tol = 1e-8)$root
  coef <- solve_at(multiplier)
  list(
    coef = coef,
    norm = rkhs_norm(coef, kernel),
    constraint = multiplier
  )
}

feature_norm <- function(coef) {
  sqrt(sum(coef^2))
}

solve_feature_with_norm_bound <- function(a, b, bound, control) {
  identity_matrix <- diag(nrow(a))
  solve_at <- function(multiplier) {
    safe_solve(
      a + 4 * multiplier * identity_matrix,
      b,
      jitter = control$jitter,
      max_tries = control$max_solve_tries
    )
  }
  coef <- solve_at(0)
  norm <- feature_norm(coef)
  if (!is.finite(bound) || norm <= bound) {
    return(list(coef = coef, norm = norm, constraint = 0))
  }

  objective <- function(multiplier) {
    feature_norm(solve_at(multiplier)) - bound
  }
  upper <- 1
  while (objective(upper) > 0 && upper < 1e12) {
    upper <- upper * 10
  }
  if (upper >= 1e12 && objective(upper) > 0) {
    stop("The RKHS norm constraint could not be bracketed.", call. = FALSE)
  }
  multiplier <- stats::uniroot(objective, c(0, upper), tol = 1e-8)$root
  coef <- solve_at(multiplier)
  list(
    coef = coef,
    norm = feature_norm(coef),
    constraint = multiplier
  )
}

nystrom_map_summary <- function(object) {
  list(
    observations = object$n_observations,
    requested_rank = object$requested_rank,
    effective_rank = object$effective_rank
  )
}

prepare_outcome_bridge_nystrom_system <- function(
    h_map, gp_map, y, weights, lambda_gp, control) {
  n_population <- sum(weights)
  phi_h <- h_map$training_features
  phi_gp <- gp_map$training_features
  weighted_phi_gp <- sweep(phi_gp, 1L, weights, "*")
  q_inner <- crossprod(phi_gp, weighted_phi_gp) / n_population +
    lambda_gp * diag(ncol(phi_gp))
  c_inner <- crossprod(phi_gp, sweep(phi_h, 1L, weights, "*")) /
    n_population
  d_inner <- drop(crossprod(phi_gp, weights * y) / n_population)
  q_solution <- safe_solve(
    q_inner,
    cbind(c_inner, d_inner),
    jitter = control$jitter,
    max_tries = control$max_solve_tries
  )
  q_solve_c <- q_solution[, seq_len(ncol(c_inner)), drop = FALSE]
  q_solve_d <- q_solution[, ncol(q_solution), drop = TRUE]
  list(
    base_matrix = crossprod(c_inner, q_solve_c),
    rhs = drop(crossprod(c_inner, q_solve_d)),
    outer_map = h_map,
    inner_map = gp_map,
    lambda_inner = lambda_gp
  )
}

finish_outcome_bridge_nystrom_system <- function(
    system, lambda_h, sigma2_h, sigma2_gp, max_norm, control) {
  a <- system$base_matrix +
    4 * lambda_h * diag(ncol(system$outer_map$training_features))
  solution <- solve_feature_with_norm_bound(
    a, system$rhs, max_norm, control
  )
  approximation <- list(
    method = "nystrom",
    outer = nystrom_map_summary(system$outer_map),
    inner = nystrom_map_summary(system$inner_map),
    landmark_sampling = control$nystrom_landmarks
  )
  prediction_map <- system$outer_map
  prediction_map$training_features <- NULL
  structure(list(
    coefficients = solution$coef,
    feature_map = prediction_map,
    kernel_approximation = "nystrom",
    approximation = approximation,
    sigma2 = sigma2_h,
    norm = solution$norm,
    constraint = solution$constraint,
    lambda_outer = lambda_h,
    lambda_inner = system$lambda_inner,
    sigma2_inner = sigma2_gp
  ), class = "pmtp_outcome_bridge")
}

fit_outcome_bridge_nystrom <- function(
    h, gp, y, weights, lambda_h, lambda_gp,
    sigma2_h, sigma2_gp, max_norm, control, feature_maps = NULL) {
  h_map <- feature_maps$outer %||% fit_nystrom_map(
    h, sigma2_h, weights, control, seed_offset = 101L
  )
  gp_map <- feature_maps$inner %||% fit_nystrom_map(
    gp, sigma2_gp, weights, control, seed_offset = 102L
  )
  system <- prepare_outcome_bridge_nystrom_system(
    h_map, gp_map, y, weights, lambda_gp, control
  )
  finish_outcome_bridge_nystrom_system(
    system, lambda_h, sigma2_h, sigma2_gp, max_norm, control
  )
}

prepare_treatment_bridge_nystrom_system <- function(
    g_map, hp_map, policy_inner_features, weights, target, policy_support,
    lambda_hp, control) {
  n_population <- sum(weights)
  phi_g <- g_map$training_features
  phi_hp <- hp_map$training_features
  weighted_support <- weights * policy_support
  q_inner <- crossprod(
    phi_hp, sweep(phi_hp, 1L, weighted_support, "*")
  ) / n_population + lambda_hp * diag(ncol(phi_hp))
  a_inner <- drop(crossprod(policy_inner_features, weights * target) /
    n_population)
  b_inner <- crossprod(
    phi_hp, sweep(phi_g, 1L, weighted_support, "*")
  ) / n_population
  q_solution <- safe_solve(
    q_inner,
    cbind(b_inner, a_inner),
    jitter = control$jitter,
    max_tries = control$max_solve_tries
  )
  q_solve_b <- q_solution[, seq_len(ncol(b_inner)), drop = FALSE]
  q_solve_a <- q_solution[, ncol(q_solution), drop = TRUE]
  list(
    base_matrix = crossprod(b_inner, q_solve_b),
    rhs = drop(crossprod(b_inner, q_solve_a)),
    outer_map = g_map,
    inner_map = hp_map,
    lambda_inner = lambda_hp
  )
}

finish_treatment_bridge_nystrom_system <- function(
    system, lambda_g, sigma2_g, sigma2_hp, max_norm, control) {
  a <- system$base_matrix +
    4 * lambda_g * diag(ncol(system$outer_map$training_features))
  solution <- solve_feature_with_norm_bound(
    a, system$rhs, max_norm, control
  )
  approximation <- list(
    method = "nystrom",
    outer = nystrom_map_summary(system$outer_map),
    inner = nystrom_map_summary(system$inner_map),
    landmark_sampling = control$nystrom_landmarks
  )
  prediction_map <- system$outer_map
  prediction_map$training_features <- NULL
  structure(list(
    coefficients = solution$coef,
    feature_map = prediction_map,
    kernel_approximation = "nystrom",
    approximation = approximation,
    sigma2 = sigma2_g,
    norm = solution$norm,
    constraint = solution$constraint,
    lambda_outer = lambda_g,
    lambda_inner = system$lambda_inner,
    sigma2_inner = sigma2_hp
  ), class = "pmtp_treatment_bridge")
}

fit_treatment_bridge_nystrom <- function(
    g, hp, hp_q, weights, target, policy_support,
    lambda_g, lambda_hp, sigma2_g, sigma2_hp, max_norm, control,
    feature_maps = NULL) {
  g_map <- feature_maps$outer %||% fit_nystrom_map(
    g, sigma2_g, weights, control, seed_offset = 201L
  )
  hp_map <- feature_maps$inner %||% fit_nystrom_map(
    hp, sigma2_hp, weights, control, seed_offset = 202L
  )
  policy_inner_features <- feature_maps$policy_inner_features %||%
    predict_nystrom_features(hp_map, hp_q)
  system <- prepare_treatment_bridge_nystrom_system(
    g_map, hp_map, policy_inner_features, weights, target, policy_support,
    lambda_hp, control
  )
  finish_treatment_bridge_nystrom_system(
    system, lambda_g, sigma2_g, sigma2_hp, max_norm, control
  )
}

fit_outcome_bridge <- function(h, gp, y, weights, lambda_h, lambda_gp,
                               sigma2_h, sigma2_gp, max_norm, control,
                               feature_maps = NULL) {
  if (identical(control$kernel_approximation, "nystrom")) {
    return(fit_outcome_bridge_nystrom(
      h, gp, y, weights, lambda_h, lambda_gp,
      sigma2_h, sigma2_gp, max_norm, control, feature_maps
    ))
  }
  n_population <- sum(weights)
  k_h <- controlled_kernel_matrix(h, sigma2 = sigma2_h, control = control)
  k_gp <- controlled_kernel_matrix(
    gp, sigma2 = sigma2_gp, control = control
  )

  k_gp_d <- sweep(k_gp, 2L, weights, "*")
  q_inner <- k_gp_d %*% k_gp / n_population + lambda_gp * k_gp
  c_inner <- k_gp_d / n_population
  inner_solution <- safe_solve(
    q_inner,
    c_inner,
    jitter = control$jitter,
    max_tries = control$max_solve_tries
  )
  minimax_metric <- crossprod(c_inner, inner_solution)
  minimax_metric <- (minimax_metric + t(minimax_metric)) / 2

  a <- k_h %*% minimax_metric %*% k_h + 4 * lambda_h * k_h
  b <- k_h %*% minimax_metric %*% y
  solution <- solve_with_norm_bound(a, b, k_h, max_norm, control)

  structure(list(
    coefficients = solution$coef,
    training_arguments = h,
    kernel_approximation = "exact",
    approximation = list(
      method = "exact",
      outer = list(observations = nrow(h), requested_rank = nrow(h),
                   effective_rank = nrow(h)),
      inner = list(observations = nrow(gp), requested_rank = nrow(gp),
                   effective_rank = nrow(gp))
    ),
    sigma2 = sigma2_h,
    kernel_family = control_kernel_family(control),
    matern_smoothness = control_matern_smoothness(control),
    norm = solution$norm,
    constraint = solution$constraint,
    lambda_outer = lambda_h,
    lambda_inner = lambda_gp,
    sigma2_inner = sigma2_gp
  ), class = "pmtp_outcome_bridge")
}

fit_treatment_bridge <- function(g, hp, hp_q, weights, target, policy_support,
                                 lambda_g, lambda_hp, sigma2_g, sigma2_hp,
                                 max_norm, control, feature_maps = NULL) {
  if (identical(control$kernel_approximation, "nystrom")) {
    return(fit_treatment_bridge_nystrom(
      g, hp, hp_q, weights, target, policy_support,
      lambda_g, lambda_hp, sigma2_g, sigma2_hp, max_norm, control,
      feature_maps
    ))
  }
  n_population <- sum(weights)
  k_g <- controlled_kernel_matrix(g, sigma2 = sigma2_g, control = control)
  k_hp <- controlled_kernel_matrix(hp, sigma2 = sigma2_hp, control = control)
  k_hp_q <- controlled_kernel_matrix(
    hp_q, hp, sigma2 = sigma2_hp, control = control
  )
  weighted_support <- weights * policy_support

  q_inner <- sweep(k_hp, 2L, weighted_support, "*") %*% k_hp /
    n_population + lambda_hp * k_hp
  a_inner <- crossprod(k_hp_q, weights * target) / n_population
  b_inner <- sweep(k_hp, 2L, weighted_support, "*") %*% k_g /
    n_population
  q_solve_a <- safe_solve(
    q_inner,
    a_inner,
    jitter = control$jitter,
    max_tries = control$max_solve_tries
  )
  q_solve_b <- safe_solve(
    q_inner,
    b_inner,
    jitter = control$jitter,
    max_tries = control$max_solve_tries
  )

  a <- crossprod(b_inner, q_solve_b) + 4 * lambda_g * k_g
  b <- crossprod(b_inner, q_solve_a)
  solution <- solve_with_norm_bound(a, b, k_g, max_norm, control)

  structure(list(
    coefficients = solution$coef,
    training_arguments = g,
    kernel_approximation = "exact",
    approximation = list(
      method = "exact",
      outer = list(observations = nrow(g), requested_rank = nrow(g),
                   effective_rank = nrow(g)),
      inner = list(observations = nrow(hp), requested_rank = nrow(hp),
                   effective_rank = nrow(hp))
    ),
    sigma2 = sigma2_g,
    kernel_family = control_kernel_family(control),
    matern_smoothness = control_matern_smoothness(control),
    norm = solution$norm,
    constraint = solution$constraint,
    lambda_outer = lambda_g,
    lambda_inner = lambda_hp,
    sigma2_inner = sigma2_hp
  ), class = "pmtp_treatment_bridge")
}

predict_outcome_bridge <- function(object, new_arguments) {
  if (identical(object$kernel_approximation, "nystrom")) {
    return(drop(
      predict_nystrom_features(object$feature_map, new_arguments) %*%
        object$coefficients
    ))
  }
  drop(kernel_matrix(
    new_arguments,
    object$training_arguments,
    sigma2 = object$sigma2,
    kernel_family = object$kernel_family %||% "gaussian",
    matern_smoothness = object$matern_smoothness %||% 2
  ) %*% object$coefficients)
}

predict_treatment_bridge <- function(object, new_arguments) {
  if (identical(object$kernel_approximation, "nystrom")) {
    return(drop(
      predict_nystrom_features(object$feature_map, new_arguments) %*%
        object$coefficients
    ))
  }
  drop(kernel_matrix(
    new_arguments,
    object$training_arguments,
    sigma2 = object$sigma2,
    kernel_family = object$kernel_family %||% "gaussian",
    matern_smoothness = object$matern_smoothness %||% 2
  ) %*% object$coefficients)
}

outcome_validation_risk <- function(residual, adversary_arguments, weights,
                                    sigma2, lambda, control,
                                    feature_map = NULL) {
  n_population <- sum(weights)
  if (identical(control$kernel_approximation, "nystrom")) {
    feature_map <- feature_map %||% fit_nystrom_map(
      adversary_arguments, sigma2, weights, control, seed_offset = 301L
    )
    features <- feature_map$training_features
    q <- crossprod(features, sweep(features, 1L, weights, "*")) /
      n_population + lambda * diag(ncol(features))
    b <- drop(crossprod(features, weights * residual) / n_population)
    solution <- safe_solve(
      q, b,
      jitter = control$jitter,
      max_tries = control$max_solve_tries
    )
    return(max(drop(crossprod(b, solution)) / 4, 0))
  }
  kernel <- controlled_kernel_matrix(
    adversary_arguments, sigma2 = sigma2, control = control
  )
  q <- sweep(kernel, 2L, weights, "*") %*% kernel /
    n_population + lambda * kernel
  b <- kernel %*% (weights * residual) / n_population
  solution <- safe_solve(
    q, b,
    jitter = control$jitter,
    max_tries = control$max_solve_tries
  )
  max(drop(crossprod(b, solution)) / 4, 0)
}

treatment_validation_risk <- function(g_value, adversary_arguments,
                                      adversary_policy_arguments, weights,
                                      target, policy_support, sigma2, lambda,
                                      control, feature_map = NULL,
                                      policy_features = NULL) {
  n_population <- sum(weights)
  if (identical(control$kernel_approximation, "nystrom")) {
    feature_map <- feature_map %||% fit_nystrom_map(
      adversary_arguments, sigma2, weights, control, seed_offset = 302L
    )
    features <- feature_map$training_features
    policy_features <- policy_features %||% predict_nystrom_features(
      feature_map, adversary_policy_arguments
    )
    weighted_support <- weights * policy_support
    q <- crossprod(
      features, sweep(features, 1L, weighted_support, "*")
    ) / n_population + lambda * diag(ncol(features))
    b <- (
      drop(crossprod(policy_features, weights * target)) -
        drop(crossprod(features, weighted_support * g_value))
    ) / n_population
    solution <- safe_solve(
      q, b,
      jitter = control$jitter,
      max_tries = control$max_solve_tries
    )
    return(max(drop(crossprod(b, solution)) / 4, 0))
  }
  kernel <- controlled_kernel_matrix(
    adversary_arguments, sigma2 = sigma2, control = control
  )
  kernel_q <- controlled_kernel_matrix(
    adversary_policy_arguments,
    adversary_arguments,
    sigma2 = sigma2,
    control = control
  )
  weighted_support <- weights * policy_support
  q <- sweep(kernel, 2L, weighted_support, "*") %*% kernel /
    n_population + lambda * kernel
  b <- (
    crossprod(kernel_q, weights * target) -
      kernel %*% (weighted_support * g_value)
  ) / n_population
  solution <- safe_solve(
    q, b,
    jitter = control$jitter,
    max_tries = control$max_solve_tries
  )
  max(drop(crossprod(b, solution)) / 4, 0)
}
