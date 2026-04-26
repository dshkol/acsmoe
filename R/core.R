result_df <- function(estimate, moe, se = NULL) {
  out <- data.frame(estimate = estimate, moe = moe)
  if (!is.null(se)) {
    out$se <- se
  }
  out
}

#' Propagate ACS uncertainty for a sum.
#'
#' @param estimates Numeric estimates to sum.
#' @param moes Numeric MOEs corresponding to `estimates`.
#' @param cov Optional covariance matrix on the standard-error scale. A scalar
#'   is interpreted as a common off-diagonal covariance.
#' @param conf Confidence level associated with input and output MOEs.
#' @details A scalar `cov` is a covariance, not a correlation. In contrast,
#'   `acs_aggregate(cov_strategy = "constant")` accepts a scalar correlation
#'   because that interface derives covariances from group-specific MOEs.
#' @return A one-row data frame with `estimate`, `moe`, and `se`.
acs_sum <- function(estimates, moes, cov = NULL, conf = 0.90) {
  validate_numeric(estimates, "estimates")
  validate_numeric(moes, "moes")
  if (length(estimates) != length(moes)) {
    stop("`estimates` and `moes` must have the same length.", call. = FALSE)
  }

  keep <- !is.na(estimates) & !is.na(moes)
  estimates <- estimates[keep]
  moes <- moes[keep]
  if (is_zero_cov(cov)) {
    # Match tidycensus::moe_sum(): multiple zero estimates sharing an MOE
    # contribute that MOE once rather than inflating the derived MOE.
    zeros <- estimates == 0
    moes_for_calc <- c(unique(moes[zeros]), moes[!zeros])
    se <- sqrt(sum(moe_to_se(moes_for_calc, conf = conf)^2))
    return(result_df(sum(estimates), se_to_moe(se, conf = conf), se))
  }
  v <- covariance_from_inputs(moes, cov, conf = conf)
  w <- rep(1, length(estimates))
  se <- sqrt(drop(t(w) %*% v %*% w))
  result_df(sum(estimates), se_to_moe(se, conf = conf), se)
}

#' Propagate ACS uncertainty for a difference.
#'
#' @param estimate1 First estimate.
#' @param moe1 MOE for `estimate1`.
#' @param estimate2 Second estimate.
#' @param moe2 MOE for `estimate2`.
#' @param cov Covariance between the two estimates on the standard-error scale.
#' @param conf Confidence level associated with input and output MOEs.
#' @return A data frame with `estimate`, `moe`, and `se`.
acs_diff <- function(estimate1, moe1, estimate2, moe2, cov = 0, conf = 0.90) {
  args <- recycle_common(
    estimate1, moe1, estimate2, moe2, cov,
    names = c("estimate1", "moe1", "estimate2", "moe2", "cov")
  )
  se1 <- moe_to_se(args$moe1, conf = conf)
  se2 <- moe_to_se(args$moe2, conf = conf)
  var <- se1^2 + se2^2 - 2 * args$cov
  if (any(var < -1e-10, na.rm = TRUE)) {
    stop("Covariance implies a negative variance.", call. = FALSE)
  }
  se <- sqrt(pmax(var, 0))
  result_df(args$estimate1 - args$estimate2, se_to_moe(se, conf = conf), se)
}

#' Propagate ACS uncertainty for a ratio.
#'
#' @param num Numerator estimate.
#' @param num_moe Numerator MOE.
#' @param denom Denominator estimate.
#' @param denom_moe Denominator MOE.
#' @param cov Numerator-denominator covariance on the standard-error scale.
#' @param conf Confidence level associated with input and output MOEs.
#' @return A data frame with `estimate`, `moe`, and `se`.
acs_ratio <- function(num, num_moe, denom, denom_moe, cov = 0, conf = 0.90) {
  args <- recycle_common(
    num, num_moe, denom, denom_moe, cov,
    names = c("num", "num_moe", "denom", "denom_moe", "cov")
  )
  se_num <- moe_to_se(args$num_moe, conf = conf)
  se_den <- moe_to_se(args$denom_moe, conf = conf)
  estimate <- args$num / args$denom
  var <- (se_num^2 / args$denom^2) +
    ((args$num^2 * se_den^2) / args$denom^4) -
    ((2 * args$num * args$cov) / args$denom^3)
  if (any(var < -1e-10, na.rm = TRUE)) {
    stop("Covariance implies a negative variance.", call. = FALSE)
  }
  se <- sqrt(pmax(var, 0))
  result_df(estimate, se_to_moe(se, conf = conf), se)
}

#' Propagate ACS uncertainty for a proportion.
#'
#' @param num Numerator estimate.
#' @param num_moe Numerator MOE.
#' @param denom Denominator estimate.
#' @param denom_moe Denominator MOE.
#' @param cov Numerator-denominator covariance on the standard-error scale.
#' @param conf Confidence level associated with input and output MOEs.
#' @details This function follows the Census approximation used for proportions
#'   where the numerator is a subset of the denominator. It does not validate
#'   that `num <= denom`; behavior for non-nested ratios is formulaic rather
#'   than a claim that the inputs define a valid proportion.
#' @return A data frame with `estimate`, `moe`, and `se`.
acs_prop <- function(num, num_moe, denom, denom_moe, cov = 0, conf = 0.90) {
  if (is_zero_cov(cov)) {
    args <- recycle_common(
      num, num_moe, denom, denom_moe,
      names = c("num", "num_moe", "denom", "denom_moe")
    )
    estimate <- args$num / args$denom
    ratio_moe <- acs_ratio(
      args$num, args$num_moe, args$denom, args$denom_moe,
      cov = 0, conf = conf
    )$moe
    x <- args$num_moe^2 - (estimate^2 * args$denom_moe^2)
    moe <- ratio_moe
    pos <- x > 0 & !is.na(x)
    moe[pos] <- sqrt(x[pos]) / args$denom[pos]
    return(result_df(estimate, moe, moe_to_se(moe, conf = conf)))
  }
  acs_ratio(num, num_moe, denom, denom_moe, cov = cov, conf = conf)
}

#' Propagate ACS uncertainty for a product.
#'
#' @param estimate1 First estimate.
#' @param moe1 MOE for `estimate1`.
#' @param estimate2 Second estimate.
#' @param moe2 MOE for `estimate2`.
#' @param cov Covariance between the two estimates on the standard-error scale.
#' @param conf Confidence level associated with input and output MOEs.
#' @return A data frame with `estimate`, `moe`, and `se`.
acs_product <- function(estimate1, moe1, estimate2, moe2, cov = 0, conf = 0.90) {
  args <- recycle_common(
    estimate1, moe1, estimate2, moe2, cov,
    names = c("estimate1", "moe1", "estimate2", "moe2", "cov")
  )
  se1 <- moe_to_se(args$moe1, conf = conf)
  se2 <- moe_to_se(args$moe2, conf = conf)
  estimate <- args$estimate1 * args$estimate2
  var <- (args$estimate2^2 * se1^2) +
    (args$estimate1^2 * se2^2) +
    (2 * args$estimate1 * args$estimate2 * args$cov)
  if (any(var < -1e-10, na.rm = TRUE)) {
    stop("Covariance implies a negative variance.", call. = FALSE)
  }
  se <- sqrt(pmax(var, 0))
  result_df(estimate, se_to_moe(se, conf = conf), se)
}

#' Propagate ACS uncertainty for a linear combination.
#'
#' @param estimates Numeric estimates.
#' @param moes Numeric MOEs corresponding to `estimates`.
#' @param weights Numeric weights for the linear combination.
#' @param cov Optional covariance matrix on the standard-error scale. A scalar
#'   is interpreted as a common off-diagonal covariance.
#' @param conf Confidence level associated with input and output MOEs.
#' @details A scalar `cov` is a covariance, not a correlation. In contrast,
#'   `acs_aggregate(cov_strategy = "constant")` accepts a scalar correlation
#'   because that interface derives covariances from group-specific MOEs.
#' @return A one-row data frame with `estimate`, `moe`, and `se`.
acs_linear <- function(estimates, moes, weights, cov = NULL, conf = 0.90) {
  validate_numeric(estimates, "estimates")
  validate_numeric(moes, "moes")
  validate_numeric(weights, "weights")
  if (length(estimates) != length(moes) || length(estimates) != length(weights)) {
    stop("`estimates`, `moes`, and `weights` must have the same length.",
         call. = FALSE)
  }
  v <- covariance_from_inputs(moes, cov, conf = conf)
  se <- sqrt(drop(t(weights) %*% v %*% weights))
  result_df(sum(weights * estimates), se_to_moe(se, conf = conf), se)
}
