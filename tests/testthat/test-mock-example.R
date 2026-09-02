test_that("the original simulated teaching dataset is installed", {
  path <- system.file("extdata", "sim_trial_data.csv", package = "proximtp")
  expect_true(nzchar(path))
  data <- utils::read.csv(path)
  expect_identical(names(data), c("Y", "L1", "L2", "L3", "A", "Z", "W", "wt"))
  expect_equal(nrow(data), 1000L)
  expect_true(all(vapply(data, is.numeric, logical(1))))
  expect_false(anyNA(data))
  expect_true(all(data$Y %in% c(0, 1)))
  expect_true(all(data$wt > 0))
  expect_true(all(data$A >= 0.39 & data$A <= 3.50))
})

test_that("the worked example is installed and parses", {
  path <- system.file("examples", "example_pmtp.R", package = "proximtp")
  expect_true(nzchar(path))
  expect_no_error(parse(path))
})
