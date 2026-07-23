#' SuperLearner library used for the non-proximal paper comparators
#'
#' Returns the exact ten-learner library used in the original simulation code.
#' When `weighted = TRUE`, the two stepwise wrappers are removed because the
#' installed SuperLearner implementations do not accept observation weights.
#' XGBoost can be added as an explicitly documented extension.
#'
#' @param weighted Whether to return only learners whose wrappers accept
#'   observation weights.
#' @param xgboost Whether to append `"SL.xgboost"`.
#'
#' @return A character vector of SuperLearner wrapper names.
#' @export
pmtp_paper_learners <- function(weighted = FALSE, xgboost = FALSE) {
  assert_flag(weighted, "weighted")
  assert_flag(xgboost, "xgboost")
  learners <- c(
    "SL.glm", "SL.glm.interaction", "SL.step.interaction",
    "SL.bayesglm", "SL.earth", "SL.mean", "SL.gam", "SL.glmnet",
    "SL.step", "SL.ranger"
  )
  if (weighted) {
    learners <- setdiff(learners, c("SL.step.interaction", "SL.step"))
  }
  if (xgboost) learners <- c(learners, "SL.xgboost")
  learners
}

validate_nonproximal_learners <- function(learners, name) {
  if (!is.character(learners) || !length(learners) || anyNA(learners) ||
      any(!nzchar(learners))) {
    stop("`", name, "` must be a nonempty character vector.", call. = FALSE)
  }
  unique(learners)
}

evaluate_nonproximal_policy <- function(policy, data, treatment) {
  if (!is.function(policy)) stop("`policy` must be a function.", call. = FALSE)
  arguments <- formals(policy)
  value <- if (length(arguments) == 1L && !"..." %in% names(arguments)) {
    policy(data[[treatment]])
  } else {
    policy(data, treatment)
  }
  if (!is.numeric(value) || length(value) != nrow(data) || anyNA(value) ||
      any(!is.finite(value))) {
    stop("`policy` must return one finite numeric value per analysis row.",
         call. = FALSE)
  }
  as.numeric(value)
}

nonproximal_analysis_data <- function(data, treatment, outcome, covariates,
                                      weights, target) {
  assert_columns(data, c(treatment, outcome, covariates), "analysis arguments")
  weight <- resolve_vector_argument(data, weights, "weights", default = 1)
  target_value <- resolve_vector_argument(data, target, "target", default = 1)
  if (!is.numeric(weight)) stop("`weights` must be numeric.", call. = FALSE)
  target_value <- normalize_indicator(target_value, "target")
  complete <- stats::complete.cases(
    data[c(treatment, outcome, covariates)], weight, target_value
  )
  keep <- complete & target_value == 1
  if (sum(keep) < 20L) {
    stop("At least 20 complete target-population rows are required.",
         call. = FALSE)
  }
  analysis <- data[keep, c(treatment, outcome, covariates), drop = FALSE]
  if (!is.numeric(analysis[[treatment]]) || !is.numeric(analysis[[outcome]])) {
    stop("Treatment and outcome must be numeric.", call. = FALSE)
  }
  if (!all(analysis[[outcome]] %in% c(0, 1))) {
    stop("The current non-proximal implementation requires a binary outcome.",
         call. = FALSE)
  }
  weight <- as.numeric(weight[keep])
  assert_positive(weight, "weights")
  list(data = analysis, weights = weight, complete = complete, keep = keep)
}

fit_pmtp_superlearner <- function(y, x, weights, ids, learners,
                                  learner_folds, seed) {
  if (length(unique(y)) < 2L) {
    stop("A SuperLearner training outcome had only one observed class.",
         call. = FALSE)
  }
  learner_folds <- min(as.integer(learner_folds), length(unique(ids)))
  if (learner_folds < 2L) {
    stop("At least two internal SuperLearner folds are required.", call. = FALSE)
  }
  cv_control <- SuperLearner::SuperLearner.CV.control(
    V = learner_folds, stratifyCV = FALSE, shuffle = TRUE
  )
  withCallingHandlers(
    withr::with_seed(seed, {
      SuperLearner::SuperLearner(
        Y = y,
        X = x,
        family = stats::binomial(),
        SL.library = learners,
        method = "method.NNLS",
        id = ids,
        obsWeights = weights / mean(weights),
        cvControl = cv_control,
        env = environment(SuperLearner::SuperLearner)
      )
    }),
    warning = function(warning) {
      if (identical(
        conditionMessage(warning),
        "non-integer #successes in a binomial glm!"
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

predict_pmtp_superlearner <- function(fit, newdata) {
  prediction <- stats::predict(fit, newdata = newdata, onlySL = TRUE)$pred
  as.numeric(prediction)
}

extract_pmtp_superlearner_diagnostics <- function(fit) {
  learners <- fit$libraryNames
  coefficient <- stats::setNames(as.numeric(fit$coef), names(fit$coef))[learners]
  cv_risk <- stats::setNames(as.numeric(fit$cvRisk), names(fit$cvRisk))[learners]
  data.frame(
    learner = learners,
    coefficient = as.numeric(coefficient),
    cv_risk = as.numeric(cv_risk),
    error_in_cv = as.logical(fit$errorsInCVLibrary),
    error_in_full_fit = as.logical(fit$errorsInLibrary),
    row.names = NULL
  )
}

make_nonproximal_features <- function(data, treatment, covariates, shifted = NULL) {
  out <- data[c(treatment, covariates)]
  if (!is.null(shifted)) out[[treatment]] <- shifted
  out
}

fit_weighted_point_nonproximal <- function(
    data, treatment, outcome, covariates, policy, weights,
    estimators, folds, learner_folds, learners_outcome,
    learners_treatment, probability_bounds, population_size, seed,
    return_fits) {
  n <- nrow(data)
  fold_id <- make_folds(data[[treatment]], data[[outcome]], folds, seed)
  q_all <- evaluate_nonproximal_policy(policy, data, treatment)
  q_natural <- q_shifted <- density_ratio <- rep(NA_real_, n)
  q_natural_tmle <- q_shifted_tmle <- rep(NA_real_, n)
  tmle_epsilon <- NA_real_
  learner_diagnostics <- vector("list", folds)
  fit_objects <- if (return_fits) vector("list", folds) else NULL

  for (fold in seq_len(folds)) {
    training <- which(fold_id != fold)
    validation <- which(fold_id == fold)
    train_data <- data[training, , drop = FALSE]
    valid_data <- data[validation, , drop = FALSE]
    x_train <- make_nonproximal_features(train_data, treatment, covariates)
    x_valid <- make_nonproximal_features(valid_data, treatment, covariates)
    x_train_shifted <- make_nonproximal_features(
      train_data, treatment, covariates, q_all[training]
    )
    x_valid_shifted <- make_nonproximal_features(
      valid_data, treatment, covariates, q_all[validation]
    )

    outcome_fit <- fit_pmtp_superlearner(
      y = train_data[[outcome]], x = x_train,
      weights = weights[training], ids = training,
      learners = learners_outcome, learner_folds = learner_folds,
      seed = seed + 1000L * fold + 1L
    )
    q0_valid <- predict_pmtp_superlearner(outcome_fit, x_valid)
    q1_valid <- predict_pmtp_superlearner(outcome_fit, x_valid_shifted)

    stacked_x <- rbind(x_train, x_train_shifted)
    stacked_y <- rep(c(0, 1), each = length(training))
    stacked_weights <- rep(weights[training], 2L)
    stacked_ids <- rep(training, 2L)
    treatment_fit <- fit_pmtp_superlearner(
      y = stacked_y, x = stacked_x, weights = stacked_weights,
      ids = stacked_ids, learners = learners_treatment,
      learner_folds = learner_folds,
      seed = seed + 1000L * fold + 2L
    )
    probability_valid <- predict_pmtp_superlearner(treatment_fit, x_valid)
    probability_valid <- pmin(pmax(
      probability_valid, probability_bounds[1L]
    ), probability_bounds[2L])
    ratio_valid <- probability_valid / (1 - probability_valid)

    learner_diagnostics[[fold]] <- list(
      outcome = extract_pmtp_superlearner_diagnostics(outcome_fit),
      treatment = extract_pmtp_superlearner_diagnostics(treatment_fit)
    )

    q_natural[validation] <- q0_valid
    q_shifted[validation] <- q1_valid
    density_ratio[validation] <- ratio_valid

    if (return_fits) {
      fit_objects[[fold]] <- list(
        outcome = outcome_fit, treatment = treatment_fit
      )
    }
  }

  if ("tmle" %in% estimators) {
    # Pooled cross-validated targeting: every initial Q and density-ratio value
    # was predicted without its row. Estimating one fluctuation on those
    # predictions makes the weighted empirical EIF residual score zero on the
    # evaluation sample and avoids using in-sample clever-covariate values.
    q0 <- pmin(pmax(q_natural, 1e-5), 1 - 1e-5)
    q1 <- pmin(pmax(q_shifted, 1e-5), 1 - 1e-5)
    fluctuation <- suppressWarnings(stats::glm(
      data[[outcome]] ~ 1,
      offset = stats::qlogis(q0),
      weights = weights * density_ratio,
      family = stats::binomial()
    ))
    tmle_epsilon <- unname(stats::coef(fluctuation)[1L])
    if (!is.finite(tmle_epsilon)) tmle_epsilon <- 0
    q_natural_tmle <- stats::plogis(stats::qlogis(q0) + tmle_epsilon)
    q_shifted_tmle <- stats::plogis(stats::qlogis(q1) + tmle_epsilon)
  }

  results <- list()
  influence <- list()
  y <- data[[outcome]]
  target_probability <- sum(weights) / population_size
  if ("sdr" %in% estimators) {
    pseudo <- q_shifted + density_ratio * (y - q_natural)
    estimate <- weighted_mean(pseudo, weights)
    centered <- (pseudo - estimate) / target_probability
    results[["sdr"]] <- c(
      estimate = estimate,
      std_error = sqrt(sum((weights * centered)^2)) / population_size
    )
    influence[["sdr"]] <- centered
  }
  if ("tmle" %in% estimators) {
    estimate <- weighted_mean(q_shifted_tmle, weights)
    pseudo <- q_shifted_tmle +
      density_ratio * (y - q_natural_tmle)
    centered <- (pseudo - estimate) / target_probability
    results[["tmle"]] <- c(
      estimate = estimate,
      std_error = sqrt(sum((weights * centered)^2)) / population_size
    )
    influence[["tmle"]] <- centered
  }
  table <- do.call(rbind, results)
  table <- data.frame(
    estimator = c(sdr = "AIPW (SDR)", tmle = "TMLE")[rownames(table)],
    estimate = table[, "estimate"],
    std_error = table[, "std_error"],
    row.names = NULL
  )
  table$conf_low <- table$estimate - stats::qnorm(0.975) * table$std_error
  table$conf_high <- table$estimate + stats::qnorm(0.975) * table$std_error
  list(
    estimates = table,
    influence_function = influence,
    predictions = list(
      outcome_natural = q_natural,
      outcome_shifted = q_shifted,
      density_ratio = density_ratio,
      outcome_natural_tmle = q_natural_tmle,
      outcome_shifted_tmle = q_shifted_tmle,
      tmle_epsilon = tmle_epsilon
    ),
    fold_id = fold_id,
    learner_diagnostics = learner_diagnostics,
    fits = fit_objects
  )
}

fit_lmtp_nonproximal <- function(
    data, treatment, outcome, covariates, policy, weights,
    estimators, folds, learner_folds, learners_outcome,
    learners_treatment, seed) {
  shift <- function(data, trt) evaluate_nonproximal_policy(policy, data, trt)
  control <- lmtp::lmtp_control(
    .learners_outcome_folds = learner_folds,
    .learners_trt_folds = learner_folds,
    .return_full_fits = FALSE
  )
  fits <- list()
  call_one <- function(estimator, seed_offset) {
    fun <- if (estimator == "sdr") lmtp::lmtp_sdr else lmtp::lmtp_tmle
    withr::with_seed(seed + seed_offset, fun(
      data = data, trt = treatment, outcome = outcome,
      baseline = covariates, k = 0, shift = shift, mtp = TRUE,
      outcome_type = "binomial", learners_outcome = learners_outcome,
      learners_trt = learners_treatment, folds = folds, weights = weights,
      control = control
    ))
  }
  for (i in seq_along(estimators)) {
    estimator <- estimators[i]
    fits[[estimator]] <- call_one(estimator, 1000L * i)
  }
  table <- do.call(rbind, lapply(names(fits), function(estimator) {
    estimate <- fits[[estimator]]$estimate
    data.frame(
      estimator = if (estimator == "sdr") "AIPW (SDR)" else "TMLE",
      estimate = estimate@x,
      std_error = estimate@std_error,
      conf_low = estimate@conf_int[1L],
      conf_high = estimate@conf_int[2L]
    )
  }))
  rownames(table) <- NULL
  list(estimates = table, fits = fits)
}

#' Estimate a non-proximal modified-treatment-policy mean
#'
#' Provides the AIPW/SDR and TMLE comparators used in the paper. Unweighted
#' analyses can call `lmtp` directly. For inverse-probability weighted
#' two-phase analyses, `engine = "weighted_point"` uses a point-treatment
#' implementation that passes sampling weights to both the outcome and stacked
#' density-ratio SuperLearners, the TMLE fluctuation, the final empirical mean,
#' and the influence-function variance.
#'
#' @param data A data frame.
#' @param treatment,outcome Single column names.
#' @param covariates Adjustment-variable column names.
#' @param policy A vectorized function of treatment, or an `lmtp`-style
#'   function of `data` and the treatment column name.
#' @param weights Optional sampling-weight column name or numeric vector.
#' @param target Optional target-population column name or 0/1 vector. Rows
#'   outside the target population are excluded before fitting.
#' @param population_size Target-population phase-one size. Defaults to the sum
#'   of the analysis weights.
#' @param estimators Any subset of `c("sdr", "tmle")`.
#' @param engine `"auto"`, `"lmtp"`, or `"weighted_point"`. Auto uses `lmtp`
#'   without sampling weights and the audited point-treatment engine otherwise.
#' @param folds Number of outer cross-fitting folds.
#' @param learner_folds Number of internal SuperLearner folds.
#' @param learners_outcome,learners_treatment SuperLearner libraries. With
#'   `NULL`, the exact paper library is used for `lmtp` and the weight-aware
#'   paper library is used for `weighted_point`.
#' @param probability_bounds Bounds applied to the stacked-classification
#'   probabilities before conversion to density ratios.
#' @param seed Random seed.
#' @param return_fits Whether to retain potentially large SuperLearner fits for
#'   the weighted engine.
#'
#' @return An object of class `pmtp_nonproximal_fit`.
#' @export
pmtp_nonproximal <- function(
    data,
    treatment = "A",
    outcome = "Y",
    covariates = c("L", "Z", "W"),
    policy,
    weights = NULL,
    target = NULL,
    population_size = NULL,
    estimators = c("sdr", "tmle"),
    engine = c("auto", "lmtp", "weighted_point"),
    folds = 5L,
    learner_folds = 10L,
    learners_outcome = NULL,
    learners_treatment = NULL,
    probability_bounds = c(0.001, 0.999),
    seed = 20260722L,
    return_fits = FALSE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  engine <- match.arg(engine)
  estimators <- match.arg(estimators, c("sdr", "tmle"), several.ok = TRUE)
  estimators <- unique(estimators)
  if (engine == "auto") {
    engine <- if (is.null(weights)) "lmtp" else "weighted_point"
  }
  weighted_engine <- identical(engine, "weighted_point")
  if (is.null(learners_outcome)) {
    learners_outcome <- pmtp_paper_learners(weighted = weighted_engine)
  }
  if (is.null(learners_treatment)) {
    learners_treatment <- learners_outcome
  }
  learners_outcome <- validate_nonproximal_learners(
    learners_outcome, "learners_outcome"
  )
  learners_treatment <- validate_nonproximal_learners(
    learners_treatment, "learners_treatment"
  )
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package `SuperLearner` is required for non-proximal estimation.",
         call. = FALSE)
  }
  if (!is.numeric(probability_bounds) || length(probability_bounds) != 2L ||
      anyNA(probability_bounds) || probability_bounds[1L] <= 0 ||
      probability_bounds[2L] >= 1 ||
      probability_bounds[1L] >= probability_bounds[2L]) {
    stop("`probability_bounds` must be increasing values inside (0, 1).",
         call. = FALSE)
  }
  for (argument in c("folds", "learner_folds")) {
    value <- get(argument)
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        value < 2L || value != as.integer(value)) {
      stop("`", argument, "` must be an integer of at least two.",
           call. = FALSE)
    }
  }
  assert_flag(return_fits, "return_fits")

  analysis <- nonproximal_analysis_data(
    data, treatment, outcome, covariates, weights, target
  )
  if (folds >= nrow(analysis$data)) {
    stop("`folds` must be smaller than the analysis sample.", call. = FALSE)
  }
  if (is.null(population_size)) population_size <- sum(analysis$weights)
  assert_positive(population_size, "population_size")
  if (length(population_size) != 1L) {
    stop("`population_size` must be a single positive number.", call. = FALSE)
  }
  if (engine == "lmtp") {
    if (!requireNamespace("lmtp", quietly = TRUE)) {
      stop("Package `lmtp` is required for `engine = \"lmtp\"`.",
           call. = FALSE)
    }
    if (!is.null(weights)) {
      warning(
        "lmtp 1.5.4 does not pass sampling weights to its SuperLearner ",
        "outcome or density-ratio fits; use `engine = \"weighted_point\"` ",
        "for the audited weighted implementation.",
        call. = FALSE
      )
    }
    result <- fit_lmtp_nonproximal(
      analysis$data, treatment, outcome, covariates, policy,
      if (is.null(weights)) NULL else analysis$weights,
      estimators, folds, learner_folds, learners_outcome,
      learners_treatment, seed
    )
    weighting <- if (is.null(weights)) "none" else "final_and_targeting_only"
  } else {
    result <- fit_weighted_point_nonproximal(
      analysis$data, treatment, outcome, covariates, policy,
      analysis$weights, estimators, folds, learner_folds,
      learners_outcome, learners_treatment, probability_bounds,
      population_size, seed, return_fits
    )
    weighting <- "nuisance_targeting_empirical_and_variance"
  }

  out <- c(result, list(
    engine = engine,
    weighting = weighting,
    n_sample = nrow(analysis$data),
    population_size = population_size,
    weights = analysis$weights,
    estimators = estimators,
    learners_outcome = learners_outcome,
    learners_treatment = learners_treatment,
    call = match.call()
  ))
  class(out) <- "pmtp_nonproximal_fit"
  out
}

#' @export
print.pmtp_nonproximal_fit <- function(x, ...) {
  cat("Non-proximal modified-treatment-policy estimates\n")
  cat("Engine:", x$engine, "\n")
  cat("Weighting:", x$weighting, "\n")
  cat("Analysis rows:", x$n_sample, "\n\n")
  print(x$estimates, row.names = FALSE)
  invisible(x)
}
