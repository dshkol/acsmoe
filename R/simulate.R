rmvnorm_eigen <- function(n, mean, cov) {
  validate_covariance_matrix(cov)
  eig <- eigen(cov, symmetric = TRUE)
  values <- pmax(eig$values, 0)
  z <- matrix(stats::rnorm(n * length(mean)), nrow = n)
  sweep(z %*% (eig$vectors %*% diag(sqrt(values), nrow = length(values))), 2,
        mean, "+")
}

simulate_draws <- function(estimates, moes, cov = NULL, n_sims = 1000,
                           dist = c("normal", "truncated_normal"),
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
  if (dist == "truncated_normal") {
    warning(
      "`truncated_normal` censors simulated draws at zero; it prevents negative values but does not preserve the input means or variances.",
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
#' @param dist Distribution assumption: `"normal"` or `"truncated_normal"`.
#' @param conf Confidence level associated with input MOEs.
#' @return A numeric matrix of simulated draws.
acs_simulate <- function(estimates, moes, cov = NULL, n_sims = 1000,
                         dist = c("normal", "truncated_normal"), conf = 0.90) {
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
#' @param dist Distribution assumption: `"normal"` or `"truncated_normal"`.
#' @param conf Confidence level associated with input MOEs.
#' @param summary Summary to return: `"mean"`, `"median"`, or `"ci"`.
#' @return A data frame summarizing the simulated derived statistic.
acs_simulate_fn <- function(estimates, moes, fn, cov = NULL, n_sims = 1000,
                            dist = c("normal", "truncated_normal"),
                            conf = 0.90,
                            summary = c("mean", "median", "ci")) {
  if (!is.function(fn)) {
    stop("`fn` must be a function.", call. = FALSE)
  }
  dist <- match.arg(dist)
  summary <- match.arg(summary)
  draws <- simulate_draws(estimates, moes, cov = cov, n_sims = n_sims,
                          dist = dist, conf = conf)
  values <- apply(draws, 1, fn)
  if (summary == "mean") {
    return(data.frame(estimate = mean(values), se = stats::sd(values)))
  }
  if (summary == "median") {
    return(data.frame(estimate = stats::median(values), se = stats::sd(values)))
  }
  qs <- stats::quantile(values, probs = c(0.05, 0.95), names = FALSE)
  data.frame(estimate = mean(values), lower = qs[1], upper = qs[2])
}
