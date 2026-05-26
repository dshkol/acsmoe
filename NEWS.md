# acsmoe 0.1.0

Initial CRAN release.

* Covariance-aware derived-estimate helpers: `acs_sum()`, `acs_diff()`,
  `acs_ratio()`, `acs_prop()`, `acs_product()`, `acs_linear()`. With no
  covariance supplied, results match the zero-covariance ACS approximation
  formulas used by `tidycensus`.
* Grouped paired estimate/MOE aggregation via `acs_aggregate()`, with
  `cov_strategy = "zero"` (default) and `cov_strategy = "constant"` for
  sensitivity analysis.
* Simulation helpers `acs_simulate()` and `acs_simulate_fn()` for derived
  quantities from estimate/MOE pairs.
* MOE/SE/CI conversion utilities: `moe_to_se()`, `se_to_moe()`, `moe_ci()`.
* Reliability diagnostics: `acs_cv()` and `acs_reliability()`.
