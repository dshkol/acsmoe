resolve_col_vector <- function(cols, data, name) {
  if (!is.character(cols)) {
    stop("`", name, "` must be a character vector of column names.",
         call. = FALSE)
  }
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("Unknown columns in `", name, "`: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  cols
}

aggregate_covariance <- function(row_ids, value_col, moes, cov_strategy,
                                 cov_value, conf) {
  if (cov_strategy == "zero") {
    return(NULL)
  }
  if (cov_strategy == "constant") {
    return(covariance_from_correlation(moes, cov_value, conf = conf))
  }
  if (cov_strategy == "supplied") {
    if (!is.list(cov_value) || is.null(cov_value[[value_col]])) {
      stop("For `cov_strategy = 'supplied'`, `cov_value` must be a named list of full covariance matrices.",
           call. = FALSE)
    }
    full <- cov_value[[value_col]]
    if (!is.matrix(full)) {
      stop("Supplied covariance for `", value_col, "` must be a matrix.",
           call. = FALSE)
    }
    return(full[row_ids, row_ids, drop = FALSE])
  }
  stop("Unsupported covariance strategy.", call. = FALSE)
}

#' Aggregate paired ACS estimate and MOE columns by group.
#'
#' @param data A data frame containing ACS estimate and MOE columns.
#' @param group_var Name of the grouping column, supplied as a single string.
#' @param value_cols Character vector of estimate column names to aggregate.
#' @param moe_cols Character vector of MOE column names paired with
#'   `value_cols`.
#' @param cov_strategy Covariance strategy: `"zero"`, `"supplied"`, or
#'   `"constant"`.
#' @param cov_value For `"constant"`, a scalar correlation. For `"supplied"`,
#'   a named list of full covariance matrices keyed by estimate column.
#' @param conf Confidence level associated with input and output MOEs.
#' @details `cov_strategy = "constant"` interprets `cov_value` as a correlation,
#'   not a covariance. This differs from scalar `cov` arguments in core
#'   propagation functions, where a scalar means an off-diagonal covariance on
#'   the standard-error scale.
#'
#'   Output rows are ordered by first appearance of each group level in `data`,
#'   not alphabetically.
#' @return A data frame with one row per group and aggregated estimate/MOE
#'   columns.
#' @examples
#' tracts <- data.frame(
#'   region = c("north", "north", "south", "south"),
#'   pop = c(1000, 1200, 900, 1100),
#'   pop_moe = c(120, 140, 100, 130)
#' )
#' acs_aggregate(tracts, "region", "pop", "pop_moe")
#' acs_aggregate(tracts, "region", "pop", "pop_moe",
#'               cov_strategy = "constant", cov_value = 0.25)
#' @export
acs_aggregate <- function(data, group_var, value_cols, moe_cols,
                          cov_strategy = c("zero", "supplied", "constant"),
                          cov_value = 0, conf = 0.90) {
  cov_strategy <- match.arg(cov_strategy)
  if (!is.character(group_var) || length(group_var) != 1L) {
    stop("`group_var` must be a single column name as a string.",
         call. = FALSE)
  }
  if (!group_var %in% names(data)) {
    stop("`group_var` must identify an existing column: '", group_var,
         "' not found.", call. = FALSE)
  }
  value_cols <- resolve_col_vector(value_cols, data, "value_cols")
  moe_cols <- resolve_col_vector(moe_cols, data, "moe_cols")
  if (length(value_cols) != length(moe_cols)) {
    stop("`value_cols` and `moe_cols` must have the same length.", call. = FALSE)
  }
  groups <- data[[group_var]]
  level_order <- unique(groups)
  level_order <- level_order[!is.na(level_order)]
  split_ids <- split(seq_len(nrow(data)), factor(groups, levels = level_order),
                     drop = TRUE)
  out <- data.frame(stats::setNames(list(names(split_ids)), group_var),
                    stringsAsFactors = FALSE)

  for (i in seq_along(value_cols)) {
    value_col <- value_cols[[i]]
    moe_col <- moe_cols[[i]]
    estimates <- numeric(length(split_ids))
    moes <- numeric(length(split_ids))
    for (g in seq_along(split_ids)) {
      ids <- split_ids[[g]]
      cov <- aggregate_covariance(
        ids, value_col, data[[moe_col]][ids],
        cov_strategy, cov_value, conf
      )
      res <- acs_sum(data[[value_col]][ids], data[[moe_col]][ids],
                     cov = cov, conf = conf)
      estimates[[g]] <- res$estimate
      moes[[g]] <- res$moe
    }
    out[[value_col]] <- estimates
    out[[moe_col]] <- moes
  }

  out
}
