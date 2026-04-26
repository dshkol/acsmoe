#' Calculate coefficient of variation from an ACS estimate and MOE.
#'
#' @param estimate Numeric estimate.
#' @param moe Numeric MOE.
#' @param conf Confidence level associated with `moe`.
#' @return Numeric coefficient of variation, using standard error divided by
#'   absolute estimate.
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
#'   thresholds. Defaults follow common ACS applied practice; adjust these for
#'   your domain, agency guidance, or project-specific quality standard.
#' @return Ordered reliability category.
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
