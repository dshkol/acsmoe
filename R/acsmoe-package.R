#' acsmoe: Propagate Uncertainty for ACS Tabular Estimates
#'
#' `acsmoe` provides utilities for American Community Survey workflows that
#' already have published estimates and margins of error. It complements
#' `tidycensus` by matching the standard zero-covariance ACS MOE formulas by
#' default while allowing users to supply covariance structures when available.
#'
#' The package does not download ACS data, estimate variance from microdata, or
#' implement regionalization. Use `tidycensus` for data access and `survey` or
#' `srvyr` for microdata/replicate-weight workflows.
#'
#' Core references include the U.S. Census Bureau ACS handbook guidance on
#' derived estimates, Walker and Herman's `tidycensus` package, and the ACS
#' uncertainty literature by Spielman, Folch, Nagle, Arribas-Bel, and
#' Koschinsky. See the README and vignettes for full references.
#'
#' @section Output shape:
#' Propagation helpers (`acs_sum`, `acs_diff`, `acs_ratio`, `acs_prop`,
#' `acs_product`, `acs_linear`) return a data frame with columns `estimate`,
#' `moe`, and `se`. `moe_ci()` returns `estimate`, `lower`, `upper`, `moe`.
#' `acs_simulate()` returns a numeric matrix of draws. `acs_simulate_fn()`
#' returns either `(estimate, se)` for `summary = "mean"`/`"median"`, or
#' `(estimate, lower, upper)` for `summary = "ci"`. `acs_cv()` and
#' `acs_reliability()` return atomic vectors.
#'
#' @keywords internal
"_PACKAGE"
