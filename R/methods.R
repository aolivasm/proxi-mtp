#' @export
print.pmtp_fit <- function(x, ...) {
  cat("Proximal modified-treatment-policy estimates\n")
  cat("Complete observations:", x$n_sample, "\n")
  cat("Estimated phase-one population size:", round(x$population_size, 2), "\n")
  cat("Weighted target probability:", round(x$target_probability, 4), "\n\n")
  print(x$estimates, row.names = FALSE)
  invisible(x)
}

#' @export
summary.pmtp_fit <- function(object, conf_level = 0.95, ...) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must lie strictly between zero and one.", call. = FALSE)
  }
  out <- object$estimates
  support <- do.call(cbind, object$policy_support)
  weighted_target <- object$weights * object$target
  out$image_proportion <- as.numeric(
    colSums(support * weighted_target) / sum(weighted_target)
  )
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  out$conf_low <- out$estimate - critical * out$std_error
  out$conf_high <- out$estimate + critical * out$std_error
  out <- out[c(
    "policy", "image_proportion", "estimate", "std_error",
    "conf_low", "conf_high"
  )]
  attr(out, "conf_level") <- conf_level
  class(out) <- c("summary_pmtp_fit", class(out))
  out
}

#' @export
print.summary_pmtp_fit <- function(x, ...) {
  level <- 100 * attr(x, "conf_level")
  cat("Proximal MTP estimates with ", level, "% confidence intervals\n\n", sep = "")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}

#' @export
coef.pmtp_fit <- function(object, ...) {
  stats::setNames(object$estimates$estimate, object$estimates$policy)
}

#' @export
vcov.pmtp_fit <- function(object, ...) {
  covariance <- crossprod(object$weights * object$influence_function) /
    object$population_size^2
  dimnames(covariance) <- list(object$estimates$policy, object$estimates$policy)
  covariance
}
