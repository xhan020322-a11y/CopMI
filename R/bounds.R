# Only a floating-point-sized upper-bound discrepancy may be moved inward.
.round_upper_bound <- function(x, upper, scale = abs(upper)) {
  upper <- rep_len(upper, length(x))
  scale <- rep_len(scale, length(x))
  if (any(!is.finite(x)) || any(!is.finite(upper))) {
    stop("Nonfinite value or upper bound.", call. = FALSE)
  }
  tolerance <- 8 * .Machine$double.eps * pmax(abs(x), abs(upper), scale,
                                             .Machine$double.xmin)
  if (any(x - upper > tolerance)) {
    stop("Value is not below its cutoff within floating-point tolerance.", call. = FALSE)
  }
  at_bound <- x >= upper
  adjustment <- numeric(length(x))
  if (any(at_bound)) {
    inside <- upper[at_bound] - .Machine$double.eps *
      pmax(abs(upper[at_bound]), .Machine$double.xmin)
    adjustment[at_bound] <- x[at_bound] - inside
    x[at_bound] <- inside
  }
  if (any(!is.finite(x)) || any(x >= upper)) {
    stop("No finite value can represent the required upper-bound adjustment.", call. = FALSE)
  }
  attr(x, "boundary") <- c(count = sum(at_bound), max_adjustment = max(c(0, adjustment)))
  x
}

# Validate generated values with a scalar or per-column positive-support rule.
.check_completed <- function(dat, X, positive = FALSE) {
  if (!is.matrix(X) || !identical(dim(X), dim(dat$X_cens))) {
    stop("Completed matrix dimensions do not match the input.", call. = FALSE)
  }
  if (!is.logical(positive) || anyNA(positive) || !length(positive) %in% c(1L, ncol(X))) {
    stop("positive must be a logical flag or one flag per column.", call. = FALSE)
  }
  X[dat$ind == 1L] <- dat$X_cens[dat$ind == 1L]
  censored <- dat$ind == 0L
  positive_cells <- censored & matrix(rep_len(positive, ncol(X)), nrow(X), ncol(X), byrow = TRUE)
  if (any(!is.finite(X[censored])) || any(X[positive_cells] <= 0)) {
    stop("Inverse transformation produced nonfinite or invalid censored values.", call. = FALSE)
  }
  limits <- matrix(dat$cutoffs, nrow(X), ncol(X), byrow = TRUE)
  checked <- .round_upper_bound(X[censored], limits[censored])
  X[censored] <- checked
  if (any(X[positive_cells] <= 0)) {
    stop("No positive value remains below the cutoff after boundary adjustment.", call. = FALSE)
  }
  dimnames(X) <- dimnames(dat$X_cens)
  attr(X, "boundary") <- attr(checked, "boundary")
  X
}
