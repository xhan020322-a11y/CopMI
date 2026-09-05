# Validate generated values; never replace a failed draw with a heuristic value.
.check_completed <- function(dat, X, positive = FALSE) {
  X[dat$ind == 1L] <- dat$X_cens[dat$ind == 1L]
  censored <- dat$ind == 0L
  if (any(!is.finite(X[censored])) || (positive && any(X[censored] <= 0))) {
    stop("Inverse transformation produced nonfinite or invalid censored values.", call. = FALSE)
  }
  limits <- matrix(dat$cutoffs, nrow(X), ncol(X), byrow = TRUE)
  if (any(X[censored] >= limits[censored])) {
    stop("An imputed value is not below its cutoff; inspect the sampler and fitted margins.", call. = FALSE)
  }
  dimnames(X) <- dimnames(dat$X_cens)
  X
}
