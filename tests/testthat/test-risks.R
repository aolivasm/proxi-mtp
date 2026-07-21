test_that("outcome risk equals direct weighted adversarial maximum", {
  x <- cbind(c(-1, 0, 1), c(0, 1, 0))
  residual <- c(0.2, -0.4, 0.1)
  weights <- c(1, 2, 3)
  sigma2 <- 0.7
  lambda <- 0.2
  control <- minimal_control()

  kernel <- gaussian_kernel(x, sigma2 = sigma2)
  n_population <- sum(weights)
  q <- sweep(kernel, 2, weights, "*") %*% kernel / n_population +
    lambda * kernel
  b <- kernel %*% (weights * residual) / n_population
  expected <- drop(crossprod(b, solve(q + diag(1e-8, 3), b))) / 4

  observed <- proximtp:::outcome_validation_risk(
    residual, x, weights, sigma2, lambda, control
  )
  expect_equal(observed, expected, tolerance = 1e-6)
})

test_that("treatment risk equals direct weighted adversarial maximum", {
  x <- cbind(c(-1, 0, 1), c(0, 1, 0))
  xq <- x
  xq[, 1] <- xq[, 1] + 0.2
  g <- c(0.8, 1.1, 0.7)
  weights <- c(1, 2, 3)
  target <- c(1, 1, 0)
  support <- c(1, 1, 1)
  sigma2 <- 0.7
  lambda <- 0.2
  control <- minimal_control()

  kernel <- gaussian_kernel(x, sigma2 = sigma2)
  kernel_q <- gaussian_kernel(xq, x, sigma2 = sigma2)
  n_population <- sum(weights)
  q <- sweep(kernel, 2, weights * support, "*") %*% kernel /
    n_population + lambda * kernel
  b <- (
    crossprod(kernel_q, weights * target) -
      kernel %*% (weights * support * g)
  ) / n_population
  expected <- drop(crossprod(b, solve(q + diag(1e-8, 3), b))) / 4

  observed <- proximtp:::treatment_validation_risk(
    g, x, xq, weights, target, support, sigma2, lambda, control
  )
  expect_equal(observed, expected, tolerance = 1e-6)
})
