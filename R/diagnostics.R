#' Calculate coefficient of variation from an ACS estimate and MOE.
#'
#' @param estimate Numeric estimate.
#' @param moe Numeric MOE.
#' @param conf Confidence level associated with `moe`.
#' @return Numeric coefficient of variation, using standard error divided by
#'   absolute estimate.
#' @examples
#' acs_cv(estimate = 1000, moe = 80)
#' acs_cv(estimate = c(1000, 100, 0), moe = c(80, 60, 5))
#' @export
acs_cv <- function(estimate, moe, conf = 0.90) {
  validate_numeric(estimate, "estimate")
  validate_numeric(moe, "moe")
  args <- recycle_common(estimate, moe, names = c("estimate", "moe"))
  se <- moe_to_se(args$moe, conf = conf)
  cv <- se / abs(args$estimate)
  cv[args$estimate == 0] <- Inf
  cv
}

#' Categorize ACS estimate reliability by CV thresholds.
#'
#' @param estimate Numeric estimate.
#' @param moe Numeric MOE.
#' @param conf Confidence level associated with `moe`.
#' @param thresholds Named numeric vector with `reliable` and `caveat` CV
#'   thresholds. Defaults of 0.12 and 0.40 reflect commonly used applied
#'   conventions (e.g., reliable below 12% CV, unreliable above 40% CV); these
#'   are not a Census Bureau standard and should be adjusted to your domain,
#'   agency guidance, or project-specific quality standard.
#' @return Ordered factor with levels `reliable`, `caveat`, `unreliable`.
#' @examples
#' acs_reliability(estimate = 1000, moe = 80)
#' acs_reliability(estimate = c(1000, 200, 50), moe = c(80, 60, 50))
#' @export
acs_reliability <- function(estimate, moe,
                            conf = 0.90,
                            thresholds = c(reliable = 0.12, caveat = 0.40)) {
  if (is.null(names(thresholds)) ||
      !all(c("reliable", "caveat") %in% names(thresholds))) {
    stop("`thresholds` must include named values `reliable` and `caveat`.",
         call. = FALSE)
  }
  cv <- acs_cv(estimate, moe, conf = conf)
  out <- rep("unreliable", length(cv))
  out[cv < thresholds[["caveat"]]] <- "caveat"
  out[cv < thresholds[["reliable"]]] <- "reliable"
  factor(out, levels = c("reliable", "caveat", "unreliable"))
}
