# acsmoe

`acsmoe` is an early R package skeleton for propagating uncertainty in American
Community Survey tabular workflows. It is scoped to estimate + MOE workflows,
not ACS data fetching or microdata variance estimation.

The current implementation includes:

- covariance-aware sums, differences, ratios, proportions, products, and linear
  combinations;
- simulation helpers for derived quantities;
- grouped geographic-style aggregation of paired estimate/MOE columns;
- MOE/SE/CI conversion utilities;
- coefficient-of-variation reliability diagnostics;
- tests against `tidycensus` formulas.

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

Run `citation("acsmoe")` after installation for full reference details.
