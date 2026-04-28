rmvnorm_eigen <- function(n, mean, cov) {
  validate_covariance_matrix(cov)
  eig <- eigen(cov, symmetric = TRUE)
  values <- pmax(eig$values, 0)
  z <- matrix(stats::rnorm(n * length(mean)), nrow = n)
  sweep(z %*% (eig$vectors %*% diag(sqrt(values), nrow = length(values))), 2,
        mean, "+")
}

simulate_draws <- function(estimates, moes, cov = NULL, n_sims = 1000,
                           dist = c("normal", "censored_normal"),
                           conf = 0.90) {
  dist <- match.arg(dist)
  validate_numeric(estimates, "estimates")
  validate_numeric(moes, "moes")
  if (length(estimates) != length(moes)) {
    stop("`estimates` and `moes` must have the same length.", call. = FALSE)
  }
  if (!is.numeric(n_sims) || length(n_sims) != 1L || n_sims < 1) {
    stop("`n_sims` must be a positive integer.", call. = FALSE)
  }
  v <- covariance_from_inputs(moes, cov, conf = conf)
  draws <- rmvnorm_eigen(as.integer(n_sims), estimates, v)
  if (dist == "censored_normal") {
    warning(
      "`censored_normal` replaces draws below zero with zero. This matches the convention used in Napierala & Denton (2017) and the ACS handbook semantics for zero estimates, but it is censoring rather than statistical truncation: it does not preserve the input mean or variance for cells whose estimate sits close to zero relative to its MOE.",
      call. = FALSE
    )
    draws <- pmax(draws, 0)
  }
  draws
}

#' Simulate ACS estimates from published estimates and MOEs.
#'
#' @param estimates Numeric estimates.
#' @param moes Numeric MOEs corresponding to `estimates`.
#' @param cov Optional covariance matrix on the standard-error scale.
#' @param n_sims Number of Monte Carlo simulations.
#' @param dist Distribution assumption: `"normal"` or `"censored_normal"`. The
#'   censored variant replaces below-zero draws with zero, matching the
#'   convention used in Napierala & Denton (2017) for ACS counts.
#' @param conf Confidence level associated with input MOEs.
#' @return A numeric matrix of simulated draws.
#' @examples
#' set.seed(1)
#' acs_simulate(c(x = 100, y = 50), c(10, 5), n_sims = 5)
#' @export
acs_simulate <- function(estimates, moes, cov = NULL, n_sims = 1000,
                         dist = c("normal", "censored_normal"), conf = 0.90) {
  draws <- simulate_draws(estimates, moes, cov, n_sims, dist, conf)
  colnames(draws) <- names(estimates)
  draws
}

#' Simulate a derived ACS statistic.
#'
#' @param estimates Numeric estimates.
#' @param moes Numeric MOEs corresponding to `estimates`.
#' @param fn Function applied to each simulated row.
#' @param cov Optional covariance matrix on the standard-error scale.
#' @param n_sims Number of Monte Carlo simulations.
#' @param dist Distribution assumption: `"normal"` or `"censored_normal"`. The
#'   censored variant replaces below-zero draws with zero, matching the
#'   convention used in Napierala & Denton (2017) for ACS counts.
#' @param conf Confidence level associated with input MOEs.
#' @param summary Summary to return: `"mean"`, `"median"`, or `"ci"`.
#' @param point For `summary = "ci"`, point estimate to report alongside the
#'   percentile interval: `"mean"` (default) or `"median"`.
#' @return A data frame summarizing the simulated derived statistic.
#' @details For `summary = "ci"`, the returned interval is the central
#'   `conf`-level percentile interval of the simulated derived values. The
#'   reported `estimate` is the chosen `point` summary of those values.
#' @examples
#' set.seed(1)
#' acs_simulate_fn(c(100, 50), c(10, 5), fn = sum, n_sims = 500)
#' set.seed(1)
#' acs_simulate_fn(c(100, 50), c(10, 5), fn = sum, n_sims = 500,
#'                 summary = "ci", conf = 0.90)
#' @export
acs_simulate_fn <- function(estimates, moes, fn, cov = NULL, n_sims = 1000,
                            dist = c("normal", "censored_normal"),
                            conf = 0.90,
                            summary = c("mean", "median", "ci"),
                            point = c("mean", "median")) {
  if (!is.function(fn)) {
    stop("`fn` must be a function.", call. = FALSE)
  }
  dist <- match.arg(dist)
  summary <- match.arg(summary)
  point <- match.arg(point)
  draws <- simulate_draws(estimates, moes, cov = cov, n_sims = n_sims,
                          dist = dist, conf = conf)
  values <- apply(draws, 1, fn)
  if (summary == "mean") {
    return(data.frame(estimate = mean(values), se = stats::sd(values)))
  }
  if (summary == "median") {
    return(data.frame(estimate = stats::median(values), se = stats::sd(values)))
  }
  validate_conf(conf)
  alpha <- (1 - conf) / 2
  qs <- stats::quantile(values, probs = c(alpha, 1 - alpha), names = FALSE)
  centre <- if (point == "median") stats::median(values) else mean(values)
  data.frame(estimate = centre, lower = qs[1], upper = qs[2])
}
