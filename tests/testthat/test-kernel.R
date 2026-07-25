test_that("Gaussian kernel uses sigma squared explicitly", {
  x <- matrix(c(0, 1), ncol = 1)
  kernel <- gaussian_kernel(x, sigma2 = 2)
  expect_equal(kernel[1, 2], exp(-1 / 4))
  expect_equal(diag(kernel), c(1, 1))
})

test_that("Matérn kernel is normalized and decreases with distance", {
  x <- matrix(c(0, 1, 2), ncol = 1)
  kernel <- matern_kernel(x, sigma2 = 1, smoothness = 2)

  expect_equal(diag(kernel), rep(1, 3))
  expect_equal(kernel, t(kernel), tolerance = 1e-14)
  expect_gt(kernel[1, 2], kernel[1, 3])
  expect_true(all(eigen(kernel, symmetric = TRUE)$values > 0))
})

test_that("weighted median bandwidth agrees with unweighted all-one weights", {
  x <- cbind(c(-1, 0, 2, 4), c(0, 1, 1, 3))
  expect_equal(
    median_bandwidth(x),
    median_bandwidth(x, rep(1, nrow(x)))
  )
  expected <- stats::median(as.numeric(stats::dist(x))^2) / 2
  expect_equal(median_bandwidth(x), expected)
})
