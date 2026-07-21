minimal_control <- function(seed = 19L) {
  pmtp_control(
    outer_folds = 2L,
    inner_folds = 2L,
    lambda_h = 1e-2,
    lambda_gp = 1,
    lambda_g = 1e-2,
    lambda_hp = 1,
    bandwidth_h = 1,
    bandwidth_gp = 1 / 2,
    bandwidth_g = 1,
    bandwidth_hp = 1 / 2,
    max_norm_h = Inf,
    max_norm_g = Inf,
    seed = seed,
    keep_cv = TRUE
  )
}

make_test_data <- function(n = 36L, seed = 7L) {
  withr::with_seed(seed, {
    l <- stats::rnorm(n)
    u <- stats::rnorm(n)
    z <- 0.5 * l + u + stats::rnorm(n)
    w <- -0.2 * l + u + stats::rnorm(n)
    a <- 0.3 * l + u + stats::rnorm(n)
    probability <- stats::plogis(-0.5 + 0.2 * l + 0.3 * a + 0.2 * u)
    data.frame(
      Y = stats::rbinom(n, 1, probability),
      A = a,
      L = l,
      Z = z,
      W = w,
      sampling_weight = seq(1.1, 3, length.out = n),
      target = as.numeric(a <= stats::quantile(a, 0.75))
    )
  })
}
