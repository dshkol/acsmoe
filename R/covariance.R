validate_covariance_matrix <- function(cov, tol = 1e-10) {
  if (!is.matrix(cov) || !is.numeric(cov)) {
    stop("Covariance input must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(cov) != ncol(cov)) {
    stop("Covariance matrix must be square.", call. = FALSE)
  }
  if (any(!is.finite(cov))) {
    stop("Covariance matrix must contain only finite values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(cov, t(cov), tolerance = tol))) {
    stop("Covariance matrix must be symmetric.", call. = FALSE)
  }
  eig <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) < -tol) {
    stop("Covariance matrix must be positive semidefinite.", call. = FALSE)
  }
  invisible(cov)
}

covariance_from_inputs <- function(moes, cov, conf = 0.90) {
  validate_numeric(moes, "moes")
  se <- moe_to_se(moes, conf = conf)
  n <- length(se)
  base <- diag(se^2, nrow = n)

  if (is.null(cov)) {
    return(base)
  }

  if (length(cov) == 1L && is.numeric(cov)) {
    out <- matrix(cov, nrow = n, ncol = n)
    diag(out) <- se^2
    validate_covariance_matrix(out)
    return(out)
  }

  if (!is.matrix(cov)) {
    stop("`cov` must be NULL, a scalar covariance, or a covariance matrix.",
         call. = FALSE)
  }
  if (!identical(dim(cov), c(n, n))) {
    stop("`cov` matrix dimensions must match the number of estimates.",
         call. = FALSE)
  }
  validate_covariance_matrix(cov)
  cov
}

covariance_from_correlation <- function(moes, rho, conf = 0.90) {
  se <- moe_to_se(moes, conf = conf)
  n <- length(se)

  if (length(rho) == 1L) {
    if (!is.numeric(rho) || is.na(rho) || rho < -1 || rho > 1) {
      stop("Scalar correlation must be between -1 and 1.", call. = FALSE)
    }
    cor <- matrix(rho, nrow = n, ncol = n)
    diag(cor) <- 1
  } else {
    cor <- rho
    if (!is.matrix(cor) || !identical(dim(cor), c(n, n))) {
      stop("Correlation matrix dimensions must match the number of estimates.",
           call. = FALSE)
    }
    if (any(abs(cor) > 1)) {
      stop("Correlations must be between -1 and 1.", call. = FALSE)
    }
  }

  cov <- cor * tcrossprod(se)
  validate_covariance_matrix(cov)
  cov
}

is_zero_cov <- function(cov) {
  is.null(cov) || (length(cov) == 1L && is.numeric(cov) && identical(as.numeric(cov), 0))
}
