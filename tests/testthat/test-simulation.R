test_that("simulation grids retain the paper settings and seeds", {
  counts <- c(main = 36, secondary = 51, weighted = 12, sobolev = 12,
              nonproximal = 36, nonproximal_weighted = 12,
              parametric = 64, parametric_misspecified = 36)
  for (suite in names(counts)) {
    grid <- pmtp_simulation_grid(suite)
    expect_equal(nrow(grid), unname(counts[suite]))
    expect_false(anyDuplicated(grid) > 0)
    expect_true(all(names(grid) %in% names(formals(pmtp_simulate))))
  }
  expect_equal(unique(pmtp_simulation_grid("main")$seed), 20280722)
  expect_equal(unique(pmtp_simulation_grid("secondary")$seed), 20260722)
  expect_equal(unique(pmtp_simulation_grid("sobolev")$seed), 20290722)
  d <- pmtp_simulation_grid("weighted")
  expect_true(all(d$weighted & d$beta_z == -d$beta_w))
})

test_that("simulation controls preserve the final tuning grids", {
  gaussian <- proximtp:::simulation_control("gaussian", 123L, FALSE)
  sobolev <- proximtp:::simulation_control("sobolev", 123L, FALSE)
  expect_equal(gaussian$outer_folds, 3)
  expect_equal(gaussian$inner_folds, 3)
  expect_equal(gaussian$inner_repeats, 2)
  expect_equal(gaussian$lambda_gp, 10^(-1:2))
  expect_equal(gaussian$bandwidth_h, 2^(-2:2))
  expect_equal(sobolev$outer_folds, 2)
  expect_equal(sobolev$inner_folds, 2)
  expect_equal(sobolev$inner_repeats, 2)
  expect_equal(sobolev$lambda_gp, 10^(-3:4))
  expect_equal(sobolev$bandwidth_h, 2^(-2:8))
  expect_equal(sobolev$kernel_family, "matern_sobolev")
  expect_equal(sobolev$nystrom_rank(1000), 300)
  expect_equal(proximtp:::simulation_control("gaussian", 123L, TRUE)$
                 weighted_loss_normalization, "horvitz_thompson")
})

test_that("simulation data retain the generation and fitting seed offsets", {
  main <- pmtp_paper_scenario("main")
  for (estimator in c("gaussian", "sobolev", "parametric", "nonproximal")) {
    dat <- proximtp:::simulation_data(80, main, estimator, FALSE, 101L)
    expect_equal(dat$data, simulate_pmtp_dgp(80, main$spec, seed = 101L))
    expect_equal(dat$fit_seed, if (estimator == "nonproximal") 200101L else 101L)
  }
  dat <- proximtp:::simulation_data(80, main, "gaussian", TRUE, 101L)
  population <- simulate_pmtp_dgp(1200, main$spec, seed = 20101L)
  expected <- sample_pmtp_two_phase(population, target_sample_size = 80,
                                    seed = 102L)$phase_two
  expect_equal(dat$data, expected)
  expect_equal(dat$population_size, 1200)
  expect_equal(dat$fit_seed, 102L)
  secondary <- pmtp_paper_scenario("c9")
  expect_equal(proximtp:::simulation_data(80, secondary, "gaussian", FALSE,
                                         101L)$fit_seed, 30101L)
})

test_that("replication subsets reproduce the same fitted estimates", {
  set.seed(871)
  before <- .Random.seed
  joint <- suppressWarnings(pmtp_simulate(40, replications = c(2, 4), seed = 123))
  single <- suppressWarnings(pmtp_simulate(40, replications = 4, seed = 123))
  expect_equal(.Random.seed, before)
  row <- joint$estimates[joint$estimates$replication == 4, ]
  rownames(row) <- NULL
  expect_equal(row, single$estimates)
  expect_true(all(is.finite(joint$estimates$estimate)))
  expect_equal(joint$estimates$seed, c(125, 127))
  expect_equal(joint$summary$successful, 2)
  expect_equal(joint$summary$mean_width, mean(joint$estimates$width))
})

test_that("simulation rejects invalid inputs and reports failed fits", {
  expect_error(pmtp_simulate(39), "n")
  expect_error(pmtp_simulate(100, replications = integer()), "replications")
  expect_error(pmtp_simulate(100, replications = c(1, 1)), "duplicates")
  expect_error(pmtp_simulate(100, seed = Inf), "seed")
  expect_error(pmtp_simulate(100, seed = .Machine$integer.max), "seed range")
  testthat::local_mocked_bindings(pmtp = function(...) stop("test failure"))
  expect_warning(out <- pmtp_simulate(40, replications = 1, seed = 1),
                 "Some fits failed")
  expect_equal(out$estimates$error, "test failure")
  expect_equal(out$summary$failed, 1)
  expect_equal(out$summary$successful, 0)
  expect_true(is.na(out$summary$coverage))
})
