# acsmoe

`acsmoe` propagates uncertainty in American Community Survey tabular
workflows. It is scoped to estimate + MOE workflows, not ACS data fetching or
microdata variance estimation.

It provides:

- covariance-aware sums, differences, ratios, proportions, products, and linear
  combinations;
- simulation helpers for derived quantities;
- grouped geographic-style aggregation of paired estimate/MOE columns;
- MOE/SE/CI conversion utilities;
- coefficient-of-variation reliability diagnostics.

With no covariance supplied, results match the standard zero-covariance ACS
approximation formulas used by `tidycensus`, so it drops into existing
estimate/MOE pipelines without changing the baseline numbers.

## Installation

```r
# install.packages("acsmoe")  # once on CRAN

# Development version
# install.packages("pak")
pak::pak("dshkol/acsmoe")
```

## References and attribution

The zero-covariance formulas match the U.S. Census Bureau's ACS guidance for
derived estimates, especially Chapter 8 of *Understanding and Using American
Community Survey Data: What All Data Users Need to Know*.

The standard R workflow and baseline MOE helpers are provided by `tidycensus`;
see Walker and Herman's package and Walker's MOE vignette:
<https://walker-data.com/tidycensus/articles/margins-of-error.html>.

The broader motivation for ACS uncertainty handling and regionalization comes
from Spielman, Folch, Nagle, Arribas-Bel, and Koschinsky's ACS uncertainty
papers. The old `geoss/censumander` repository is cited as historical reference
material only; regionalization is not part of this package.
