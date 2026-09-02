#' Reproduce one simulation setting
#'
#' Runs the paper's estimators on independently generated datasets. Gaussian
#' fits use three outer folds and twice-repeated three-fold tuning; the Sobolev
#' sensitivity uses two outer folds and twice-repeated two-fold tuning.
#'
#' @param n Sample size, or target phase-two sample size when `weighted = TRUE`.
#' @param beta_z,beta_w Proxy coefficients in [pmtp_dgp_spec()].
#' @param scenario Data-generating scenario; see [pmtp_paper_scenario()].
#' @param estimator Estimator family. `nonproximal` returns AIPW and TMLE;
#'   `parametric` returns the parametric bridge estimators.
#' @param weighted Use the two-phase sampling design with cohort size `15 * n`.
#' @param replications Positive integer replication IDs. Subsets can be run
#'   separately without changing the random seeds for any replication.
#' @param seed Base seed. Defaults are 20280722 for primary/weighted Gaussian
#'   simulations, 20260722 for secondary Gaussian simulations, 20290722 for
#'   Sobolev, 20260750 for non-proximal, and 20300723 for parametric simulations.
#'   Replication `r` uses `seed + r`, with fixed offsets for sampling and fitting.
#'
#' @return A list with `estimates` (including failed fits), `summary`,
#'   `settings` (including the seed), and `session_info`. Coverage and width
#'   summaries use finite intervals; failure counts are reported separately.
#' @seealso [pmtp_simulation_grid()]
#' @export
#' @examples
#' \dontrun{
#' result <- pmtp_simulate(n = 750, beta_z = 0.5, beta_w = -0.5,
#'                         replications = 1:500, seed = 20280722)
#' result$summary
#' }
pmtp_simulate <- function(
    n, beta_z = 2, beta_w = -2,
    scenario = c("main", "c6", "c7", "c8", "c9"),
    estimator = c("gaussian", "sobolev", "nonproximal", "parametric"),
    weighted = FALSE, replications = seq_len(500), seed = NULL) {
  scenario <- match.arg(scenario)
  estimator <- match.arg(estimator)
  assert_flag(weighted, "weighted")
  simulation_integer(n, "n", minimum = 40)
  simulation_integer(replications, "replications", vector = TRUE)
  if (anyDuplicated(replications)) {
    stop("`replications` must not contain duplicates.", call. = FALSE)
  }
  if (estimator == "parametric" && !scenario %in% c("main", "c8")) {
    stop("Parametric simulations support scenarios main and c8.", call. = FALSE)
  }
  if (is.null(seed)) {
    seed <- switch(estimator, gaussian = if (scenario == "main") {
      20280722L
    } else 20260722L, sobolev = 20290722L,
    nonproximal = 20260750L, parametric = 20300723L)
  }
  simulation_integer(seed, "seed", minimum = 0)
  if (as.double(seed) + max(replications) + 200000 > .Machine$integer.max) {
    stop("The seed and replication IDs exceed the supported seed range.",
         call. = FALSE)
  }
  scientific <- pmtp_paper_scenario(scenario, beta_z, beta_w)
  settings <- list(n = n, beta_z = beta_z, beta_w = beta_w,
                   scenario = scenario, estimator = estimator,
                   weighted = weighted, replications = replications, seed = seed)
  rows <- lapply(replications, function(replication) {
    replication_seed <- as.integer(seed + replication)
    result <- tryCatch({
      dat <- simulation_data(n, scientific, estimator, weighted, replication_seed)
      if (estimator %in% c("gaussian", "sobolev")) {
        fit <- pmtp(
          dat$data, policy = scientific$policy,
          policy_support = scientific$policy_support, target = scientific$target,
          weights = dat$weights, population_size = dat$population_size,
          control = simulation_control(estimator, dat$fit_seed, weighted)
        )
        out <- fit$estimates[c("estimate", "std_error", "conf_low", "conf_high")]
        out$method <- "proximal_DR"
      } else if (estimator == "nonproximal") {
        fit <- pmtp_nonproximal(
          dat$data, policy = scientific$policy, target = scientific$target,
          weights = dat$weights, population_size = dat$population_size,
          engine = "weighted_point", folds = 5, learner_folds = 5,
          learners_outcome = pmtp_paper_learners(
            weighted = weighted || !is.null(scientific$target)
          ), estimators = c("sdr", "tmle"), seed = dat$fit_seed
        )
        out <- fit$estimates[c("estimate", "std_error", "conf_low", "conf_high")]
        out$method <- c("AIPW", "TMLE")
      } else {
        fit <- if (scenario == "main") {
          pmtp_parametric_suite(dat$data, scientific$spec, weights = dat$weights,
                                population_size = dat$population_size)
        } else {
          pmtp_parametric(dat$data, scientific$spec, weights = dat$weights,
                          population_size = dat$population_size)
        }
        valid <- if (scenario == "main") {
          cv <- fit$converged
          c(cv[["h_correct"]], cv[["h_misspecified"]],
            cv[["g_correct"]], cv[["g_misspecified"]],
            cv[["h_correct"]] && cv[["g_correct"]],
            cv[["h_correct"]] && cv[["g_misspecified"]],
            cv[["h_misspecified"]] && cv[["g_correct"]],
            cv[["h_misspecified"]] && cv[["g_misspecified"]])
        } else c(fit$converged[["h"]], fit$converged[["g"]], all(fit$converged))
        out <- data.frame(method = names(fit$estimates),
                          estimate = unname(fit$estimates),
                          std_error = unname(fit$standard_error))
        out$estimate[!valid] <- out$std_error[!valid] <- NA_real_
        out$conf_low <- out$estimate - stats::qnorm(0.975) * out$std_error
        out$conf_high <- out$estimate + stats::qnorm(0.975) * out$std_error
      }
      out$error <- ifelse(is.finite(out$estimate) & is.finite(out$std_error),
                           NA_character_, "Fit did not converge")
      out
    }, error = function(e) {
      data.frame(method = simulation_methods(estimator, scenario),
                 estimate = NA_real_, std_error = NA_real_, conf_low = NA_real_,
                 conf_high = NA_real_, error = conditionMessage(e))
    })
    result$replication <- replication
    result$seed <- replication_seed
    result$truth <- scientific$truth
    result$covered <- result$conf_low <= scientific$truth &
      result$conf_high >= scientific$truth
    result$width <- result$conf_high - result$conf_low
    result
  })
  estimates <- do.call(rbind, rows)
  rownames(estimates) <- NULL
  if (any(!is.na(estimates$error))) {
    warning("Some fits failed; see `estimates$error` and summary failure counts.",
            call. = FALSE)
  }
  list(estimates = estimates, summary = simulation_summary(estimates),
       settings = settings, session_info = utils::sessionInfo())
}

simulation_integer <- function(x, name, minimum = 1, vector = FALSE) {
  if (!is.numeric(x) || !length(x) || (!vector && length(x) != 1L) ||
      anyNA(x) || any(!is.finite(x)) || any(x < minimum | x != floor(x)) ||
      any(x > .Machine$integer.max)) {
    stop("`", name, "` must contain integers >= ", minimum, ".", call. = FALSE)
  }
}

simulation_control <- function(estimator, seed, weighted) {
  sobolev <- estimator == "sobolev"
  pmtp_control(
    outer_folds = if (sobolev) 2L else 3L,
    inner_folds = if (sobolev) 2L else 3L, inner_repeats = 2L,
    lambda_h = 10^(-5:-1), lambda_g = 10^(-5:-1),
    lambda_gp = if (sobolev) 10^(-3:4) else 10^(-1:2),
    lambda_hp = if (sobolev) 10^(-3:4) else 10^(-1:2),
    bandwidth_h = if (sobolev) 2^(-2:8) else 2^(-2:2),
    bandwidth_g = if (sobolev) 2^(-2:8) else 2^(-2:2),
    bandwidth_gp = 1 / 4, bandwidth_hp = 1 / 4,
    critical_radius_rule = if (sobolev) "matern_sobolev" else "gaussian_dimension",
    kernel_family = if (sobolev) "matern_sobolev" else "gaussian",
    sobolev_l = 4, matern_smoothness = 2,
    kernel_approximation = "nystrom",
    nystrom_rank = pmtp_nystrom_rank(exponent = 2 / 3, multiplier = 3, min_rank = 30),
    nystrom_landmarks = "uniform", cache_kernel_features = TRUE,
    weighted_loss_normalization = if (weighted) "horvitz_thompson" else "hajek",
    max_norm_h = 50, max_norm_g = 50, selection_rule = "minimum",
    keep_cv = FALSE, progress = FALSE, seed = seed
  )
}

simulation_data <- function(n, scientific, estimator, weighted, seed) {
  secondary <- scientific$scenario != "main" &&
    estimator %in% c("gaussian", "sobolev")
  nonproximal <- estimator == "nonproximal"
  population_size <- if (weighted) 15 * n else n
  data_seed <- seed + if (!weighted || nonproximal) 0L else if (secondary) {
    10000L
  } else 20000L
  data <- simulate_pmtp_dgp(population_size, scientific$spec, seed = data_seed)
  if (weighted) {
    sampling_seed <- seed + if (nonproximal) 100000L else if (secondary) {
      20000L
    } else 1L
    data <- sample_pmtp_two_phase(data, target_sample_size = n,
                                  seed = sampling_seed)$phase_two
  }
  fit_seed <- seed + if (nonproximal) 200000L else if (secondary) {
    30000L
  } else if (weighted) 1L else 0L
  list(data = data, weights = if (weighted) "ipw" else NULL,
       population_size = population_size, fit_seed = fit_seed)
}

simulation_methods <- function(estimator, scenario) {
  if (estimator == "nonproximal") return(c("AIPW", "TMLE"))
  if (estimator != "parametric") return("proximal_DR")
  if (scenario != "main") return(c("OR", "DQW", "DR"))
  c("OR_h_correct", "OR_h_misspecified", "DQW_g_correct", "DQW_g_misspecified",
    "DR_h_correct_g_correct", "DR_h_correct_g_misspecified",
    "DR_h_misspecified_g_correct", "DR_h_misspecified_g_misspecified")
}

simulation_summary <- function(estimates) {
  out <- lapply(split(estimates, estimates$method), function(x) {
    valid <- is.na(x$error) & is.finite(x$estimate) & is.finite(x$width)
    good <- x[valid, , drop = FALSE]
    data.frame(method = x$method[1], attempted = nrow(x), successful = sum(valid),
      failed = sum(!valid), bias = mean(good$estimate - good$truth),
      empirical_sd = stats::sd(good$estimate),
      mean_estimated_variance = mean(good$std_error^2),
      coverage = mean(good$covered), mean_width = mean(good$width))
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Simulation settings and seeds used in the paper
#'
#' @param suite Simulation suite.
#' @return A data frame. Each row supplies the arguments for one call to
#'   [pmtp_simulate()], excluding `replications` (500 by default).
#' @export
#' @examples
#' pmtp_simulation_grid("sobolev")
pmtp_simulation_grid <- function(suite = c(
    "main", "secondary", "weighted", "sobolev", "nonproximal",
    "nonproximal_weighted", "parametric", "parametric_misspecified")) {
  suite <- match.arg(suite)
  sizes <- c(750L, 1500L, 3000L, 6000L)
  strengths <- c(2, 1, 0.5)
  grid <- function(n, bz, bw, scenario = "main", diagonal = FALSE) {
    d <- expand.grid(n = n, beta_z = bz, beta_w = bw, KEEP.OUT.ATTRS = FALSE)
    if (diagonal) d <- d[d$beta_z == -d$beta_w, , drop = FALSE]
    d$scenario <- scenario
    d
  }
  if (suite == "secondary") {
    result <- rbind(
      grid(sizes[1:3], 2, -2, "c6"),
      grid(sizes, strengths, -strengths, "c7", TRUE),
      grid(sizes, c(3, 2, 1.5, 1, 0.75, 0.5),
           -c(3, 2, 1.5, 1, 0.75, 0.5), "c8", TRUE),
      grid(sizes, strengths, -strengths, "c9", TRUE)
    )
  } else if (suite %in% c("parametric", "parametric_misspecified")) {
    if (suite == "parametric") strengths <- c(strengths, 0.25)
    result <- grid(c(750L, 3000L, 12000L, 48000L), strengths, -strengths,
                    if (suite == "parametric") "main" else "c8")
  } else {
    result <- grid(sizes, strengths, -strengths,
                    diagonal = suite %in% c("weighted", "sobolev",
                                             "nonproximal_weighted"))
  }
  result$estimator <- switch(suite, sobolev = "sobolev",
    nonproximal = "nonproximal", nonproximal_weighted = "nonproximal",
    parametric = "parametric", parametric_misspecified = "parametric", "gaussian")
  result$weighted <- suite %in% c("weighted", "nonproximal_weighted")
  result$seed <- switch(suite, secondary = 20260722L, sobolev = 20290722L,
    nonproximal = 20260750L, nonproximal_weighted = 20260750L,
    parametric = 20300723L, parametric_misspecified = 20300723L, 20280722L)
  result <- result[order(result$n, result$scenario, -result$beta_z,
                          -abs(result$beta_w)), , drop = FALSE]
  rownames(result) <- NULL
  result
}
