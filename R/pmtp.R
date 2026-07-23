resolve_vector_argument <- function(data, value, name, default = NULL) {
  if (is.null(value)) {
    if (is.null(default)) return(NULL)
    return(rep(default, nrow(data)))
  }
  if (is.character(value) && length(value) == 1L) {
    assert_columns(data, value, name)
    return(data[[value]])
  }
  if (length(value) != nrow(data)) {
    stop(
      "`", name, "` must be a column name or have one value per data row.",
      call. = FALSE
    )
  }
  value
}

normalize_indicator <- function(x, name) {
  if (is.logical(x)) x <- as.numeric(x)
  if (!is.numeric(x) || anyNA(x) || any(!x %in% c(0, 1))) {
    stop("`", name, "` must contain only 0/1 or FALSE/TRUE.", call. = FALSE)
  }
  as.numeric(x)
}

normalize_policies <- function(policy) {
  if (is.function(policy)) policy <- list(policy)
  if (!is.list(policy) || !length(policy) ||
      !all(vapply(policy, is.function, logical(1)))) {
    stop("`policy` must be a function or a nonempty list of functions.", call. = FALSE)
  }
  bad <- vapply(policy, function(fun) {
    arguments <- formals(fun)
    !length(arguments) %in% c(1L, 2L) || "..." %in% names(arguments)
  }, logical(1))
  if (any(bad)) {
    stop(
      "Each policy must be a vectorized function of treatment, or a function ",
      "of the analysis data and treatment column name.",
      call. = FALSE
    )
  }
  if (is.null(names(policy))) names(policy) <- rep("", length(policy))
  blank <- !nzchar(names(policy))
  names(policy)[blank] <- paste0("policy", which(blank))
  names(policy) <- make.unique(names(policy))
  policy
}

evaluate_pmtp_policy <- function(policy, data, treatment) {
  value <- if (length(formals(policy)) == 1L) {
    policy(data[[treatment]])
  } else {
    policy(data, treatment)
  }
  if (!is.numeric(value) || length(value) != nrow(data) ||
      anyNA(value) || any(!is.finite(value))) {
    stop("Policy must return one finite numeric value per treatment.",
         call. = FALSE)
  }
  as.numeric(value)
}

infer_policy_support <- function(a, target, policy, policy_name) {
  if (length(formals(policy)) != 1L) {
    stop(
      "Policy `", policy_name, "` depends on the analysis data. Supply its ",
      "`policy_support` indicator explicitly.",
      call. = FALSE
    )
  }
  target_a <- a[target == 1]
  if (length(target_a) < 2L) {
    stop("At least two target-population observations are required.", call. = FALSE)
  }
  grid <- seq(min(target_a), max(target_a), length.out = 1001L)
  q_grid <- policy(grid)
  if (!is.numeric(q_grid) || length(q_grid) != length(grid) ||
      anyNA(q_grid) || any(!is.finite(q_grid))) {
    stop("Policy `", policy_name, "` did not return finite numeric values.", call. = FALSE)
  }
  differences <- diff(q_grid)
  tolerance <- sqrt(.Machine$double.eps) * max(1, diff(range(q_grid)))
  monotone <- all(differences >= -tolerance) || all(differences <= tolerance)
  if (!monotone) {
    stop(
      "Policy `", policy_name, "` is not monotone on the empirical target support. ",
      "Supply its `policy_support` indicator explicitly.",
      call. = FALSE
    )
  }
  as.numeric(a >= min(q_grid) - tolerance & a <= max(q_grid) + tolerance)
}

normalize_policy_support <- function(policy_support, policies, a, target,
                                     analysis_data, treatment,
                                     complete_rows, original_n) {
  if (is.null(policy_support)) {
    return(Map(
      function(fun, name) infer_policy_support(a, target, fun, name),
      policies, names(policies)
    ))
  }
  if (!is.list(policy_support)) policy_support <- list(policy_support)
  if (length(policy_support) == 1L && length(policies) > 1L) {
    policy_support <- rep(policy_support, length(policies))
  }
  if (length(policy_support) != length(policies)) {
    stop("`policy_support` must have one element per policy.", call. = FALSE)
  }

  lapply(seq_along(policy_support), function(i) {
    support <- policy_support[[i]]
    if (is.function(support)) {
      support <- if (length(formals(support)) == 1L) {
        support(a)
      } else if (length(formals(support)) == 2L) {
        support(analysis_data, treatment)
      } else {
        stop(
          "Policy-support functions must take treatment, or data and the ",
          "treatment column name.",
          call. = FALSE
        )
      }
    } else if (length(support) == original_n) {
      support <- support[complete_rows]
    }
    if (length(support) != length(a)) {
      stop(
        "Policy-support element ", i,
        " must be a function of treatment or one value per data row.",
        call. = FALSE
      )
    }
    normalize_indicator(support, paste0("policy_support[[", i, "]]"))
  })
}

build_analysis_data <- function(data, treatment, outcome, covariates,
                                negative_control_treatment,
                                negative_control_outcome, weights, target,
                                policies, policy_support) {
  original_n <- nrow(data)
  a_all <- data[[treatment]]
  y_all <- data[[outcome]]
  l_all <- as_numeric_matrix(data, covariates, "covariates")
  z_all <- as_numeric_matrix(
    data, negative_control_treatment, "negative_control_treatment"
  )
  w_proxy_all <- as_numeric_matrix(
    data, negative_control_outcome, "negative_control_outcome"
  )
  weight_all <- resolve_vector_argument(data, weights, "weights", default = 1)
  target_all <- resolve_vector_argument(data, target, "target", default = 1)

  if (!is.numeric(a_all) || !is.numeric(y_all) || !is.numeric(weight_all)) {
    stop("Treatment, outcome, and weights must be numeric.", call. = FALSE)
  }
  complete <- stats::complete.cases(
    a_all, y_all, l_all, z_all, w_proxy_all, weight_all, target_all
  )
  if (!all(complete)) {
    message("Removing ", sum(!complete), " observations with incomplete analysis data.")
  }

  a <- as.numeric(a_all[complete])
  y <- as.numeric(y_all[complete])
  l <- l_all[complete, , drop = FALSE]
  z <- z_all[complete, , drop = FALSE]
  w_proxy <- w_proxy_all[complete, , drop = FALSE]
  weight <- as.numeric(weight_all[complete])
  target_value <- normalize_indicator(target_all[complete], "target")
  analysis_data <- data[complete, , drop = FALSE]

  assert_positive(weight, "weights")
  if (length(a) < 8L) {
    stop("At least eight complete observations are required.", call. = FALSE)
  }

  q <- lapply(seq_along(policies), function(i) {
    tryCatch(
      evaluate_pmtp_policy(policies[[i]], analysis_data, treatment),
      error = function(error) {
        stop("Policy `", names(policies)[i], "`: ", conditionMessage(error),
             call. = FALSE)
      }
    )
  })
  names(q) <- names(policies)
  support <- normalize_policy_support(
    policy_support, policies, a, target_value, analysis_data, treatment,
    complete, original_n
  )
  names(support) <- names(policies)

  colnames(l) <- paste0("L", seq_len(ncol(l)))
  colnames(z) <- paste0("Z", seq_len(ncol(z)))
  colnames(w_proxy) <- paste0("W", seq_len(ncol(w_proxy)))
  core_names <- c("A", colnames(l), colnames(z), colnames(w_proxy))
  a_index <- 1L
  l_index <- seq.int(2L, length.out = ncol(l))
  z_index <- seq.int(max(l_index) + 1L, length.out = ncol(z))
  w_index <- seq.int(max(z_index) + 1L, length.out = ncol(w_proxy))

  list(
    a = a,
    y = y,
    l = l,
    z = z,
    w_proxy = w_proxy,
    weight = weight,
    target = target_value,
    q = q,
    policy_support = support,
    policy_names = names(policies),
    complete_rows = complete,
    original_n = original_n,
    core_names = core_names,
    index = list(a = a_index, l = l_index, z = z_index, w = w_index)
  )
}

#' Estimate a proximal modified-treatment-policy mean
#'
#' Learns outcome and treatment bridge functions in Gaussian RKHSs. Tuning is
#' performed by inner cross-validation entirely within each outer-training
#' sample; selected bridges are refitted on that complete outer-training sample
#' and evaluated on its held-out fold.
#'
#' When inverse two-phase sampling weights are supplied, the same weighted
#' empirical minimax objectives are used for bridge fitting and validation. The
#' regularization rates depend on the number of observed rows in each fold, not
#' on the sum of inverse-probability weights. This keeps the nuisance fits
#' invariant to a common rescaling of all weights. The
#' final estimator is
#' `sum(weights * phi) / sum(weights * target)`. Its variance uses squared
#' inverse sampling weights, as required by the observed-data influence
#' function.
#'
#' @param data A data frame.
#' @param treatment,outcome Single column names.
#' @param covariates Column names for measured covariates.
#' @param negative_control_treatment Column name(s) for negative control
#'   treatment variables.
#' @param negative_control_outcome Column name(s) for negative control outcome
#'   variables.
#' @param policy A vectorized function of treatment; an `lmtp`-style function
#'   of the analysis data and treatment column name; or a named list of such
#'   functions. Data-dependent policies require an explicit `policy_support`.
#' @param weights Optional column name or numeric vector of inverse inclusion
#'   probabilities. Do not normalize these weights to mean one.
#' @param target Optional column name or 0/1 vector indicating the target
#'   population. Defaults to all observations.
#' @param policy_support Optional function, vector, or list defining
#'   `I_q(A,L)`. A function may take treatment alone or the analysis data and
#'   treatment column name. When omitted, support is inferred only for
#'   treatment-only policies that are monotone over the empirical target
#'   support.
#' @param population_size Full phase-one cohort size. By default, uses the sum
#'   of analysis weights. Supplying this explicitly is recommended for a
#'   two-phase application.
#' @param control A control object created by [pmtp_control()].
#'
#' @return An object of class `pmtp_fit`.
#' @export
pmtp <- function(
    data,
    treatment = "A",
    outcome = "Y",
    covariates = "L",
    negative_control_treatment = "Z",
    negative_control_outcome = "W",
    policy,
    weights = NULL,
    target = NULL,
    policy_support = NULL,
    population_size = NULL,
    control = pmtp_control()) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!inherits(control, "pmtp_control")) {
    stop("`control` must be created by `pmtp_control()`.", call. = FALSE)
  }
  required <- c(
    treatment, outcome, covariates,
    negative_control_treatment, negative_control_outcome
  )
  assert_columns(data, required, "analysis arguments")
  policies <- normalize_policies(policy)
  dat <- build_analysis_data(
    data, treatment, outcome, covariates,
    negative_control_treatment, negative_control_outcome,
    weights, target, policies, policy_support
  )

  n_sample <- length(dat$a)
  if (control$outer_folds >= n_sample) {
    stop("`outer_folds` must be smaller than the complete-case sample.", call. = FALSE)
  }
  smallest_outer_training <- n_sample - ceiling(n_sample / control$outer_folds)
  if (control$tune && control$inner_folds >= smallest_outer_training) {
    stop("Too many inner folds for the outer-training sample size.", call. = FALSE)
  }
  if (is.null(population_size)) population_size <- sum(dat$weight)
  assert_positive(population_size, "population_size")
  if (length(population_size) != 1L) {
    stop("`population_size` must be a single positive value.", call. = FALSE)
  }

  outer_fold <- make_folds(
    dat$a, dat$y, control$outer_folds, control$seed
  )
  n_policy <- length(policies)
  h0 <- rep(NA_real_, n_sample)
  hq <- matrix(NA_real_, n_sample, n_policy, dimnames = list(NULL, names(policies)))
  g0 <- matrix(NA_real_, n_sample, n_policy, dimnames = list(NULL, names(policies)))
  phi <- matrix(NA_real_, n_sample, n_policy, dimnames = list(NULL, names(policies)))
  tuning <- vector("list", control$outer_folds)

  for (fold in seq_len(control$outer_folds)) {
    training <- which(outer_fold != fold)
    validation <- which(outer_fold == fold)
    if (control$progress) message("Outer fold ", fold, " of ", control$outer_folds)

    h_cv <- if (control$tune) {
      tune_outcome_bridge(dat, training, control, seed_offset = 1000L * fold)
    } else {
      fixed_outcome_tuning(control)
    }
    h_final <- fit_selected_outcome(
      dat, training, validation, h_cv$selected, control
    )
    h0[validation] <- h_final$h0
    tuning[[fold]] <- list(outcome = h_cv, treatment = vector("list", n_policy))

    for (policy_index in seq_len(n_policy)) {
      g_cv <- if (control$tune) {
        tune_treatment_bridge(
          dat, training, policy_index, control,
          seed_offset = 1000L * fold + 10L * policy_index
        )
      } else {
        fixed_treatment_tuning(control)
      }
      prepared_policy <- prepare_fold_data(
        dat, training, validation, policy_index
      )
      g_final <- fit_selected_treatment(
        dat, training, validation, policy_index,
        g_cv$selected, prepared_policy, control
      )
      hq[validation, policy_index] <- predict_outcome_bridge(
        h_final$fit, prepared_policy$validation$hq
      )
      g0[validation, policy_index] <- g_final$g0
      phi[validation, policy_index] <-
        dat$target[validation] * hq[validation, policy_index] +
        dat$policy_support[[policy_index]][validation] *
          g0[validation, policy_index] *
          (dat$y[validation] - h0[validation])
      g_cv$final <- list(
        norm = g_final$fit$norm,
        constraint = g_final$fit$constraint,
        bandwidth = g_final$base_bandwidth,
        critical_radius_dimension = g_final$critical_radius_dimension,
        approximation = g_final$fit$approximation
      )
      tuning[[fold]]$treatment[[policy_index]] <- g_cv
      names(tuning[[fold]]$treatment)[policy_index] <- names(policies)[policy_index]
    }
    tuning[[fold]]$outcome$final <- list(
      norm = h_final$fit$norm,
      constraint = h_final$fit$constraint,
      bandwidth = h_final$base_bandwidth,
      critical_radius_dimension = h_final$critical_radius_dimension,
      approximation = h_final$fit$approximation
    )
  }

  denominator <- sum(dat$weight * dat$target)
  if (!is.finite(denominator) || denominator <= 0) {
    stop("The weighted target population is empty.", call. = FALSE)
  }
  estimate <- colSums(dat$weight * phi) / denominator
  target_probability <- denominator / population_size
  influence <- (phi - outer(dat$target, estimate)) / target_probability

  # The observed-data IF is Delta * S * influence. Because only Delta = 1
  # rows are present here, E[(Delta S influence)^2] is estimated with S^2.
  asymptotic_variance <- colSums((dat$weight * influence)^2) / population_size
  standard_error <- sqrt(asymptotic_variance / population_size)
  critical <- stats::qnorm(0.975)
  estimates <- data.frame(
    policy = names(policies),
    estimate = as.numeric(estimate),
    std_error = as.numeric(standard_error),
    conf_low = as.numeric(estimate - critical * standard_error),
    conf_high = as.numeric(estimate + critical * standard_error),
    row.names = NULL
  )

  if (!control$keep_cv) {
    tuning <- lapply(tuning, function(fold) {
      fold$outcome$results <- NULL
      fold$outcome$fold_id <- NULL
      fold$outcome$risk_folds <- NULL
      fold$treatment <- lapply(fold$treatment, function(x) {
        x$results <- NULL
        x$fold_id <- NULL
        x$risk_folds <- NULL
        x
      })
      fold
    })
  }

  structure(list(
    estimates = estimates,
    asymptotic_variance = stats::setNames(asymptotic_variance, names(policies)),
    influence_function = influence,
    contributions = phi,
    nuisance = list(h0 = h0, hq = hq, g0 = g0),
    tuning = tuning,
    outer_fold = outer_fold,
    weights = dat$weight,
    target = dat$target,
    policy_support = dat$policy_support,
    complete_rows = dat$complete_rows,
    n_sample = n_sample,
    population_size = population_size,
    target_probability = target_probability,
    weighted = !is.null(weights),
    control = control,
    call = match.call()
  ), class = "pmtp_fit")
}
