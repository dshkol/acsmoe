acs_z <- function(conf) {
  validate_conf(conf)
  stats::qnorm((1 + conf) / 2)
}

validate_conf <- function(conf) {
  if (!is.numeric(conf) || length(conf) != 1L || is.na(conf) ||
      conf <= 0 || conf >= 1) {
    stop("`conf` must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(conf)
}

validate_numeric <- function(x, name) {
  if (!is.numeric(x)) {
    stop("`", name, "` must be numeric.", call. = FALSE)
  }
  invisible(x)
}

recycle_common <- function(..., names) {
  args <- list(...)
  lengths <- vapply(args, length, integer(1))
  n <- max(lengths)
  bad <- lengths != 1L & lengths != n
  if (any(bad)) {
    stop("Inputs must have length 1 or a common length.", call. = FALSE)
  }
  args <- Map(rep, args, length.out = n)
  stats::setNames(args, names)
}

#' Convert an ACS margin of error to a standard error.
#'
#' @param moe Numeric margin of error.
#' @param conf Confidence level associated with `moe`.
#' @return Numeric standard error.
moe_to_se <- function(moe, conf = 0.90) {
  validate_numeric(moe, "moe")
  moe / acs_z(conf)
}

#' Convert a standard error to an ACS-style margin of error.
#'
#' @param se Numeric standard error.
#' @param conf Desired confidence level.
#' @return Numeric margin of error.
se_to_moe <- function(se, conf = 0.90) {
  validate_numeric(se, "se")
  se * acs_z(conf)
}

#' Convert an ACS MOE to lower and upper confidence bounds.
#'
#' @param estimate Numeric estimate.
#' @param moe Numeric margin of error.
#' @param conf_in Confidence level associated with `moe`.
#' @param conf_out Desired confidence level for the returned interval.
#' @return A data frame with `estimate`, `lower`, `upper`, and `moe`.
moe_ci <- function(estimate, moe, conf_in = 0.90, conf_out = 0.95) {
  validate_numeric(estimate, "estimate")
  validate_numeric(moe, "moe")
  args <- recycle_common(estimate, moe, names = c("estimate", "moe"))
  se <- moe_to_se(args$moe, conf = conf_in)
  out_moe <- se_to_moe(se, conf = conf_out)
  data.frame(
    estimate = args$estimate,
    lower = args$estimate - out_moe,
    upper = args$estimate + out_moe,
    moe = out_moe
  )
}
