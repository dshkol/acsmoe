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
#' @keywords internal
"_PACKAGE"
