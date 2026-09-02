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
remotes::install_github("aolivasm/proxi-mtp", ref = "main")
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

## Worked example with mock data

The package includes the simulated dataset from the original worked
example, [sim_trial_data.csv](inst/extdata/sim_trial_data.csv).
It contains 1,000 simulated observations and no real participant data.
The variables are the binary outcome `Y`, treatment `A`, covariates
`L1`, `L2`, and `L3`, negative controls `Z` and `W`, and sampling
weights `wt`.

The complete [example_pmtp.R](inst/examples/example_pmtp.R) script uses
the original treatment-only and covariate-dependent policies, with
argument names updated for the current package. The main steps are
explained below.

```r
dat <- read.csv(system.file(
  "extdata", "sim_trial_data.csv", package = "proximtp", mustWork = TRUE
))
```

### Modifying the policy for the entire population

As in the original example, treatment values in the mock data lie
between 0.39 and 3.50. A constant increase could assign values above
3.50. We therefore use a policy that adds a specified amount at lower
treatment values and gradually reduces the increase near the upper
endpoint. This illustrates the policy-modification strategy described
above, with the entire population as the target.

```r
mock_lower <- 0.39
mock_upper <- 3.50

taper <- function(a, delta) {
  ifelse(
    a + delta <= mock_upper - 1,
    a + delta,
    (delta * mock_upper + a) / (delta + 1)
  )
}
policy_q1 <- function(a) taper(a, 0.4)
policy_q2 <- function(a) taper(a, 0.8)

# A small tuning grid and low rank for demonstration, not a final analysis.
quick_control <- pmtp_control_fast(
  seed = 1234,
  critical_radius_rule = "gaussian_dimension",
  kernel_approximation = "nystrom",
  nystrom_rank = 60
)

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
summary(fit)
```

The two estimates are the counterfactual mean outcomes under `q1` and
`q2`. Sampling weights are supplied through `weights = "wt"`. The mock
file does not include a known phase-one cohort size, so the function
uses the sum of those weights as its default `population_size`.

For these monotone, treatment-only policies, the package can calculate
the policy image from the empirical treatment range; no
`policy_support` argument is needed in this first example.

The deliberately small grid can produce warnings that the selected
configuration lies on a grid boundary. Inspect `fit$tuning` and consider
expanding the relevant grid for a substantive analysis.

### Restricting the target population

The other strategy is to retain a constant shift and restrict the
population for which the counterfactual mean is estimated. For example,
consider a constant increase of 0.5 units and take as the target population
individuals whose observed treatment is at most 3.00. For these
individuals, the assigned value does not exceed 3.50.

Their original treatment values range from 0.39 to 3.00. After adding
0.5, the values that the policy can assign range from 0.89 to 3.50.
This latter range is called the **image of the policy**. It describes
possible assigned treatment values; the target population remains the
individuals with observed treatment at most 3.00.

The `target` argument specifies that population. The estimator also
uses the policy image when estimating the treatment bridge. If supplied,
`policy_support` is a function returning 1 when an observed treatment
value lies in that image and 0 otherwise. For this example:

```r
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
summary(restricted_fit)
```

This estimates the mean outcome that would be observed among individuals
with treatment at most 3.00 if each received 0.5 additional units.
Keep the full dataset in `data`: observations outside the target
population can still contribute to bridge estimation.

Here the image indicator is supplied explicitly using the known
mock-data endpoints. For a monotone treatment-only policy it can instead
be omitted, letting the package approximate the image from the empirical
treatment range in the target population. Neither approach establishes
the causal support assumptions; those require scientific justification.

### A covariate-dependent policy

The original example also allows the increase to depend on `L1`: the
base increase is 0.4 or 0.8, with an additional 0.2 when `L1 > 50`.
These policies again taper near 3.50 and target the entire population.
Use `function(data, treatment)` to access both the exposure and covariates:

```r
policy_q3 <- function(data, treatment) {
  delta <- 0.4 + 0.2 * (data$L1 > 50)
  taper(data[[treatment]], delta)
}
policy_q4 <- function(data, treatment) {
  delta <- 0.8 + 0.2 * (data$L1 > 50)
  taper(data[[treatment]], delta)
}
```

For a given value of `L1`, a policy with increase `delta` can assign
treatment values from `0.39 + delta` through `3.50`. Because this range
now depends on a covariate, provide a corresponding image indicator
for each policy:

```r
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
summary(covariate_fit)
```

For multiple policies, supply image indicators in the same order as
the policy functions. To run the complete worked example:

```r
source(system.file(
  "examples", "example_pmtp.R", package = "proximtp", mustWork = TRUE
))
```

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
probabilities, supply the phase-two observations and their original
inverse inclusion probabilities, as illustrated by `weights = "wt"`
in the mock example. When the phase-one cohort size is known, also
supply it through `population_size`.

To use the weighted empirical objectives from the paper, set
`weighted_loss_normalization = "horvitz_thompson"` in
`pmtp_control()`. Do not rescale inverse inclusion probabilities to
sum to one. Weights enter bridge fitting, validation, preprocessing,
and final estimation. Set `weights = NULL` for unweighted estimation
when sampling weights are not needed.

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
