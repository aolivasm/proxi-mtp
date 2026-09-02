# Worked example using the original, entirely simulated teaching dataset.
# No real participant data are included.
library(proximtp)

dat <- read.csv(system.file(
  "extdata", "sim_trial_data.csv", package = "proximtp", mustWork = TRUE
))

# Exposure endpoints in the mock-data construction.
mock_lower <- 0.39
mock_upper <- 3.50

# Original policies: a fixed increase that tapers near the upper endpoint.
taper <- function(a, delta) {
  ifelse(
    a + delta <= mock_upper - 1,
    a + delta,
    (delta * mock_upper + a) / (delta + 1)
  )
}
policy_q1 <- function(a) taper(a, 0.4)
policy_q2 <- function(a) taper(a, 0.8)

# Small tuning grid and low rank for demonstration, not a final analysis.
quick_control <- pmtp_control_fast(
  seed = 1234,
  critical_radius_rule = "gaussian_dimension",
  kernel_approximation = "nystrom",
  nystrom_rank = 60
)

# Weighted estimation for the two treatment-only policies.
# No known phase-one size accompanies the mock file; the default uses sum(wt).
fit <- pmtp(
  data = dat,
  treatment = "A",
  outcome = "Y",
  covariates = c("L1", "L2", "L3"),
  negative_control_treatment = "Z",
  negative_control_outcome = "W",
  weights = "wt",
  policy = list(q1 = policy_q1, q2 = policy_q2),
  control = quick_control
)
print(summary(fit))

# An untapered shift of 0.5, restricted to individuals with A <= 3.
# Original target values: [0.39, 3]; policy-assigned values: [0.89, 3.5].
shift_image <- function(a) as.numeric(a >= 0.89 & a <= 3.5)
restricted_fit <- pmtp(
  data = dat,
  treatment = "A",
  outcome = "Y",
  covariates = c("L1", "L2", "L3"),
  negative_control_treatment = "Z",
  negative_control_outcome = "W",
  weights = "wt",
  policy = list(shift_05 = function(a) a + 0.5),
  target = as.numeric(dat$A <= 3),
  policy_support = list(shift_05 = shift_image),
  control = quick_control
)
print(summary(restricted_fit))

# Original covariate-dependent policies, updated to the two-argument API.
policy_q3 <- function(data, treatment) {
  delta <- 0.4 + 0.2 * (data$L1 > 50)
  taper(data[[treatment]], delta)
}
policy_q4 <- function(data, treatment) {
  delta <- 0.8 + 0.2 * (data$L1 > 50)
  taper(data[[treatment]], delta)
}
image_q3 <- function(data, treatment) {
  delta <- 0.4 + 0.2 * (data$L1 > 50)
  as.numeric(data[[treatment]] >= mock_lower + delta &
               data[[treatment]] <= mock_upper)
}
image_q4 <- function(data, treatment) {
  delta <- 0.8 + 0.2 * (data$L1 > 50)
  as.numeric(data[[treatment]] >= mock_lower + delta &
               data[[treatment]] <= mock_upper)
}

covariate_fit <- pmtp(
  data = dat,
  treatment = "A",
  outcome = "Y",
  covariates = c("L1", "L2", "L3"),
  negative_control_treatment = "Z",
  negative_control_outcome = "W",
  weights = "wt",
  policy = list(q3 = policy_q3, q4 = policy_q4),
  policy_support = list(q3 = image_q3, q4 = image_q4),
  control = quick_control
)
print(summary(covariate_fit))

# Use weights = NULL for an unweighted fit when sampling weights are not needed.
# When the phase-one cohort size is known, also supply population_size.
