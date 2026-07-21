rkhs_norm <- function(coef, kernel) {
  value <- drop(crossprod(coef, kernel %*% coef))
  sqrt(max(value, 0))
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

fit_outcome_bridge <- function(h, gp, y, weights, lambda_h, lambda_gp,
                               sigma2_h, sigma2_gp, max_norm, control) {
  n_population <- sum(weights)
  k_h <- gaussian_kernel(h, sigma2 = sigma2_h)
  k_gp <- gaussian_kernel(gp, sigma2 = sigma2_gp)

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
    sigma2 = sigma2_h,
    norm = solution$norm,
    constraint = solution$constraint,
    lambda_outer = lambda_h,
    lambda_inner = lambda_gp,
    sigma2_inner = sigma2_gp
  ), class = "pmtp_outcome_bridge")
}

fit_treatment_bridge <- function(g, hp, hp_q, weights, target, policy_support,
                                 lambda_g, lambda_hp, sigma2_g, sigma2_hp,
                                 max_norm, control) {
  n_population <- sum(weights)
  k_g <- gaussian_kernel(g, sigma2 = sigma2_g)
  k_hp <- gaussian_kernel(hp, sigma2 = sigma2_hp)
  k_hp_q <- gaussian_kernel(hp_q, hp, sigma2 = sigma2_hp)
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
    sigma2 = sigma2_g,
    norm = solution$norm,
    constraint = solution$constraint,
    lambda_outer = lambda_g,
    lambda_inner = lambda_hp,
    sigma2_inner = sigma2_hp
  ), class = "pmtp_treatment_bridge")
}

predict_outcome_bridge <- function(object, new_arguments) {
  drop(gaussian_kernel(
    new_arguments,
    object$training_arguments,
    sigma2 = object$sigma2
  ) %*% object$coefficients)
}

predict_treatment_bridge <- function(object, new_arguments) {
  drop(gaussian_kernel(
    new_arguments,
    object$training_arguments,
    sigma2 = object$sigma2
  ) %*% object$coefficients)
}

outcome_validation_risk <- function(residual, adversary_arguments, weights,
                                    sigma2, lambda, control) {
  n_population <- sum(weights)
  kernel <- gaussian_kernel(adversary_arguments, sigma2 = sigma2)
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
                                      control) {
  n_population <- sum(weights)
  kernel <- gaussian_kernel(adversary_arguments, sigma2 = sigma2)
  kernel_q <- gaussian_kernel(
    adversary_policy_arguments,
    adversary_arguments,
    sigma2 = sigma2
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
