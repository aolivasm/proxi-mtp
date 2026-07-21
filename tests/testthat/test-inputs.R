test_that("nonmonotone policies require an explicit support indicator", {
  data <- make_test_data()
  nonmonotone <- function(a) sin(a)
  expect_error(
    pmtp(data, policy = nonmonotone, control = minimal_control()),
    "not monotone"
  )

  support <- function(a) rep(1, length(a))
  expect_no_error(
    pmtp(
      data,
      policy = nonmonotone,
      policy_support = support,
      control = minimal_control()
    )
  )
})

test_that("fold-specific weighted scaling centers training variables", {
  data <- make_test_data()
  policies <- proximtp:::normalize_policies(list(identity = function(a) a))
  dat <- proximtp:::build_analysis_data(
    data, "A", "Y", "L", "Z", "W", NULL, NULL,
    policies, NULL
  )
  prepared <- proximtp:::prepare_fold_data(dat, 1:24, 25:36)
  core <- proximtp:::make_core_matrix(dat, 1:24)
  scaled <- proximtp:::apply_weighted_scaler(core, prepared$scaler)
  expect_equal(unname(colMeans(scaled)), rep(0, ncol(scaled)), tolerance = 1e-12)
})
