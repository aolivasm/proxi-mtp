# proximtp

An R package for proximal estimation of modified-treatment-policy means,
with Gaussian and Matérn–Sobolev RKHS bridges, cross-fitting, and support
for two-phase sampling.

## Installation

```r
install.packages("remotes")
remotes::install_github("aolivasm/proxi-mtp", ref = "v0.1.1")
```

## Simulations

```r
library(proximtp)
result <- pmtp_simulate(
  n = 750, beta_z = 0.5, beta_w = -0.5,
  replications = 1:500, seed = 20280722
)
result$summary
```

`pmtp_simulation_grid()` supplies the settings and seeds for the primary,
secondary, weighted, Sobolev, non-proximal, and parametric experiments.
For example, to reproduce every primary setting:

```r
settings <- pmtp_simulation_grid("main")
results <- lapply(seq_len(nrow(settings)), function(i) {
  do.call(pmtp_simulate, as.list(settings[i, ]))
})
```

See `?pmtp_simulate`, `?pmtp_simulation_grid`, and `?pmtp` for arguments.
Non-proximal simulations additionally require `SuperLearner`, `arm`,
`earth`, `gam`, `glmnet`, and `ranger`.
