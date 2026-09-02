# proximtp

An R package for proximal causal inference for modified treatment policies
(MTPs). The package estimates counterfactual outcome means using negative
control variables to address unmeasured confounding. It implements
doubly robust cross-fitted estimation with Gaussian or Matérn–Sobolev
reproducing kernel Hilbert space (RKHS) bridge functions, including
inverse-probability-weighted estimation under two-phase sampling.

## Modified treatment policies

A modified treatment policy specifies a rule $q(A,L)$ that defines a
hypothetical treatment assignment for each individual, depending on the
observed treatment $A$ and potentially on covariates $L$. The goal of an
MTP analysis is to estimate the expected outcome in a hypothetical
scenario where each individual receives the treatment specified by this
policy instead of their observed treatment.

For example, a simple policy could be a shift intervention that
increases the treatment $A$ by a fixed amount:

$$
q(a,l) = a + 0.5.
$$

This corresponds to a policy that hypothetically shifts each
individual's treatment by 0.5 units. For some individuals, the assigned
treatment under this policy may fall outside the support of the
observed treatment conditional on their covariates, preventing
identification of the counterfactual MTP mean for the entire population
because of a positivity violation. This issue is distinct from the
presence of unmeasured confounding and can arise even in its absence.
One strategy to address it is to restrict the target population to
individuals for whom the counterfactual mean under the policy is
identifiable. Another strategy is to modify the policy itself so that
the counterfactual mean can be identified for the entire population.
Both strategies are described in
[Olivas-Martinez et al. (2025)](https://arxiv.org/abs/2512.12038).

Let $Y(q)$ denote the counterfactual outcome that would have been
observed had, instead of the actual treatment $A$, the individual
received the treatment level $q(A,L)$. Under a user-specified policy,
the function `pmtp()` estimates the mean of these counterfactual
outcomes within the population defined by $(A,L)\in\mathcal S$:

$$
\psi = \mathbb{E}\left[Y(q)\mid (A,L)\in\mathcal S\right].
$$

For the shift policy above, this is the mean outcome that would be
observed in the target population if each individual's treatment were
increased by 0.5 units. If the outcome is binary, the mean is the
probability of that outcome under the policy. The target population is
defined using the observed treatment and covariates, before applying
the policy. By default, it includes the entire population; the `target`
argument allows the user to restrict it. The estimand is a
counterfactual mean, not a contrast between two policies.

## Methodological context

The method uses an outcome bridge $h(A,L,W)$ and a treatment bridge
$g(A,L,Z)$. These are solutions to the observed bridge equations, not
ordinary outcome and treatment regressions. You specify the observed
variables and the policy; the package estimates and tunes the bridges.
The doubly robust estimating equation identifies the target when either
bridge is valid, under the identification assumptions. Valid negative
controls, bridge existence, and the relevant support conditions are
essential; Wald inference additionally requires the convergence and
regularity conditions in the [methodological paper](https://arxiv.org/abs/2512.12038).

## Installation

```r
install.packages("remotes")
remotes::install_github("aolivasm/proxi-mtp", ref = "v0.1.1")
library(proximtp)
```

## Specifying the data and estimator

The main function is `pmtp()`. Column names below are defaults: replace
them with the names in your own data.

| Argument | What to provide |
| --- | --- |
| `data` | A data frame, with one row per observation. |
| `treatment` | The numeric exposure column, e.g. `"A"`. |
| `outcome` | The numeric outcome column, e.g. `"Y"`; binary outcomes should be coded 0/1. |
| `covariates` | Names of measured adjustment covariates, e.g. `c("L1", "L2")`. |
| `negative_control_treatment` | Names of treatment-inducing proxy variables, e.g. `"Z"`. |
| `negative_control_outcome` | Names of outcome-inducing proxy variables, e.g. `"W"`. |
| `policy` | A vectorized function or a named list of such functions; see examples below. |
| `target` | A column name or a 0/1 vector indicating membership in $\mathcal S$. Default: all observations. |
| `policy_support` | Indicator of membership in the policy's image, supplied as a function or 0/1 vector, or a list in the same order as `policy`. See below. |
| `weights` | A column name or vector of positive inverse inclusion probabilities. Default: equal weights. |
| `population_size` | The phase-one cohort size when using two-phase sampling weights. |
| `control` | A `pmtp_control()` object specifying cross-fitting, tuning, and kernel options. |

Covariates and proxies must also be numeric; encode categorical variables
as appropriate numeric indicators before fitting. Rows with missing
analysis variables are removed, with a message. This complete-case
handling is not a general adjustment for missing data. The estimator
accepts binary or continuous outcomes, but does not implement
censoring adjustments for time-to-event outcomes.

Policies are evaluated on the original variable scales. Each function
must return one finite numeric value per row and must have either one
argument, `function(a)`, or two arguments, `function(data, treatment)`.
The latter receives the data frame and the exposure column name.

## A self-contained example

Generate an artificial dataset. Although the generator also returns
latent variables for simulation purposes, only the observed variables
are supplied to the estimator.

```r
dat <- simulate_pmtp_dgp(n = 200, seed = 20260902)
dat <- dat[c("A", "Y", "L", "Z", "W")]
```

In this data-generating mechanism, the exposure has support $[-2,2]$.
A constant positive shift would move some values beyond that support.
Instead, the following policy increases exposure by 0.4 units at lower
values and tapers the increase near the upper boundary. Its image is
$[-1.6,2]$.

```r
shift <- function(a) {
  ifelse(a <= 0.6, a + 0.4, a + (0.4 / 1.4) * (2 - a))
}
shift_image <- function(a) as.numeric(a >= -1.6 & a <= 2)

# A small tuning grid for learning the interface, not the paper simulations.
quick_control <- pmtp_control_fast(
  seed = 1234,
  critical_radius_rule = "gaussian_dimension"
)

fit <- pmtp(
  data = dat,
  treatment = "A",
  outcome = "Y",
  covariates = "L",
  negative_control_treatment = "Z",
  negative_control_outcome = "W",
  policy = list(shift_04 = shift),
  policy_support = list(shift_04 = shift_image),
  control = quick_control
)
summary(fit)
```

The deliberately small grid can produce warnings that the selected
configuration lies on a grid boundary. Inspect `fit$tuning` and consider
expanding the relevant grid for a substantive analysis.

### Target population versus policy image

These arguments have different roles:

- `target` identifies the observations whose counterfactual outcomes are
  averaged: $I\{(A,L)\in\mathcal S\}$.
- `policy_support` evaluates whether the **observed** exposure–covariate
  pair belongs to $\mathcal S_q=\{(q(a,l),l):(a,l)\in\mathcal S\}$.
  It is used in the treatment bridge equation and estimating function.

For example, a constant shift of 0.4 can instead be considered among
individuals with $A\leq1.6$. In this example its image is again $[-1.6,2]$:

```r
restricted_fit <- pmtp(
  data = dat,
  policy = list(constant_shift = function(a) a + 0.4),
  target = as.numeric(dat$A <= 1.6),
  policy_support = list(constant_shift = shift_image),
  control = quick_control
)
```

Keep non-target observations in `data`: they can still contribute to
bridge estimation. Restricting `target` changes the estimand.

For monotone, treatment-only policies, omitting `policy_support` asks the
package to approximate the image using the empirical exposure range in
the target population. This is not a check of the causal support
assumptions. Supply the image explicitly when its definition is known;
this is required for covariate-dependent policies.

### A covariate-dependent policy

Use the two-argument interface to access covariates. Here the increase
depends on `L`, and both the policy and its image use that same rule:

```r
covariate_shift <- function(data, treatment) {
  a <- data[[treatment]]
  delta <- ifelse(data$L < 0, 0.6, 0.4)
  ifelse(a + delta <= 1, a + delta, a + delta / (delta + 1) * (2 - a))
}
covariate_image <- function(data, treatment) {
  delta <- ifelse(data$L < 0, 0.6, 0.4)
  as.numeric(data[[treatment]] >= -2 + delta & data[[treatment]] <= 2)
}

covariate_fit <- pmtp(
  data = dat,
  policy = list(covariate_shift = covariate_shift),
  policy_support = list(covariate_shift = covariate_image),
  control = quick_control
)
```

For multiple policies, provide a named list of policy functions and a
corresponding list of image indicators in the same order.

## Bridge estimation, cross-fitting, and tuning

For each outer cross-fitting fold, the package uses all observations
outside that fold to estimate the bridge functions. It selects
hyperparameters by inner cross-validation within those training data,
refits the selected bridges on the complete outer training sample, and
then evaluates them on the omitted outer fold. Repeated inner partitions
can be requested with `inner_repeats`.

Standardization is computed within each fitting sample and applied to
its validation or evaluation observations. Nyström approximations and
reuse of kernel computations reduce the cost of tuning. Nyström
landmarks define a low-rank kernel representation; they do not replace
the fitting sample with a subsample.

For example, the following specifies the Gaussian configuration used
in the primary simulations:

```r
gaussian_control <- pmtp_control(
  outer_folds = 3,
  inner_folds = 3,
  inner_repeats = 2,
  kernel_family = "gaussian",
  critical_radius_rule = "gaussian_dimension",
  kernel_approximation = "nystrom",
  nystrom_rank = pmtp_nystrom_rank(
    exponent = 2 / 3, multiplier = 3, min_rank = 30
  ),
  selection_rule = "minimum",
  seed = 1234
)
# Use control = gaussian_control in a pmtp() call.
```

The `lambda_h`, `lambda_g`, `lambda_gp`, and `lambda_hp` arguments
specify candidate penalty multipliers; `bandwidth_h`, `bandwidth_g`,
`bandwidth_gp`, and `bandwidth_hp` specify median-heuristic bandwidth
multipliers. The package computes the fold-specific penalties and
bandwidths. See `?pmtp_control` for the grids, scaling rules, norm
constraints, and numerical options.

For Matérn–Sobolev bridges, specify
`kernel_family = "matern_sobolev"` and
`critical_radius_rule = "matern_sobolev"`; the smoothness settings
are `matern_smoothness` and `sobolev_l`, with
`sobolev_l = 2 * matern_smoothness`. The simulation interface below
automatically supplies the configuration for the paper's Sobolev
sensitivity analysis, including its two-fold cross-fitting.

Set controls explicitly when reproducing an analysis: not every
`pmtp_control()` default matches the primary simulation configuration.
The reduced grids in `pmtp_control_fast()` are intended for examples
and software checks.

## Two-phase sampling

When exposure is measured in a subsample selected with known
probabilities, supply the phase-two observations, their original
inverse inclusion probabilities, and the phase-one cohort size.
For example, using artificial data:

```r
phase_one <- simulate_pmtp_dgp(n = 3000, seed = 20260903)
phase_two <- sample_pmtp_two_phase(
  phase_one, target_sample_size = 200, seed = 20260904
)$phase_two

weighted_control <- pmtp_control_fast(
  seed = 1234,
  critical_radius_rule = "gaussian_dimension",
  weighted_loss_normalization = "horvitz_thompson"
)
weighted_fit <- pmtp(
  data = phase_two[c("A", "Y", "L", "Z", "W", "ipw")],
  policy = list(shift_04 = shift),
  policy_support = list(shift_04 = shift_image),
  weights = "ipw",
  population_size = nrow(phase_one),
  control = weighted_control
)
```

Do not rescale the inverse inclusion probabilities to sum to one.
Weights enter bridge fitting, validation, preprocessing, and final
estimation. The small tuning grid above is again for demonstration;
use the weighted simulation interface to reproduce the paper settings.

## Returned results

`pmtp()` returns a `pmtp_fit` object. Useful components and methods are:

```r
summary(fit)                 # Estimates, standard errors, and 95% Wald intervals
summary(fit, conf_level = 0.90)
coef(fit)                    # One estimated mean per policy
vcov(fit)                    # Covariance matrix across policy estimates
fit$estimates                # Results as a data frame
fit$nuisance                 # Cross-fitted bridge evaluations
fit$tuning                   # Selected tuning configurations
```

In `fit$nuisance`, `h0` evaluates the outcome bridge at observed
exposures, `hq` evaluates it at policy-modified exposures, and `g0`
contains the treatment bridge evaluations. Candidate-level validation
results are retained when `keep_cv = TRUE`.

The `image_proportion` in `summary(fit)` is the weighted proportion of
target-population observations whose observed values lie in the policy
image. It is a descriptive quantity, not a positivity diagnostic.

## Reproducing the simulations

`pmtp()` analyzes a supplied dataset; `pmtp_simulate()` generates
datasets and repeats a specified simulation experiment. The proxy
parameters below define the artificial data-generating mechanism and
are not inputs required to analyze an observed dataset.

```r
result <- pmtp_simulate(
  n = 750, beta_z = 0.5, beta_w = -0.5,
  scenario = "main", estimator = "gaussian",
  replications = 1:500, seed = 20280722
)
result$summary
```

Use `replications = 1` for a single replication. Replication identifiers
and the base `seed` determine the random seeds, so a subset of
replications can be reproduced without running the preceding ones.
The returned object includes estimates, summary statistics, settings,
seeds, and session information.

`pmtp_simulation_grid()` supplies the settings and seeds for each
experiment. For example, to reproduce all primary settings:

```r
settings <- pmtp_simulation_grid("main")
settings
results <- lapply(seq_len(nrow(settings)), function(i) {
  do.call(pmtp_simulate, as.list(settings[i, ]))
})
```

Other suites are `"secondary"`, `"weighted"`, `"sobolev"`,
`"nonproximal"`, `"nonproximal_weighted"`, `"parametric"`, and
`"parametric_misspecified"`. Full grids with 500 replications per
setting can require substantial computing time. See `?pmtp_simulate`
and `?pmtp_simulation_grid` for the scenarios and reproducibility
specifications.

The package also provides `pmtp_nonproximal()` for non-proximal AIPW
and TMLE comparisons and `pmtp_parametric_general()` for parametric
bridge estimation. Non-proximal simulations additionally require
`SuperLearner`, `arm`, `earth`, `gam`, `glmnet`, and `ranger`.

## Reference

Olivas-Martinez, A., Gilbert, P. B., and Rotnitzky, A. (2025).
[Proximal Causal Inference for Modified Treatment Policies](https://arxiv.org/abs/2512.12038).
arXiv:2512.12038.
