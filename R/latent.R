#' Transform fitted margins to latent Gaussian scores
#'
#' Apply `qnorm(F_j(x))` using the selected fitted CDF for each variable.
#' @param object A `copmi_margins` object from [copmi_fit_margins()].
#' @return A `copmi_latent` list with `Z` (numeric matrix of Gaussian scores),
#'   `ind` (observation indicators), `z_lod` (latent cutoffs; `Inf` for absent
#'   cutoffs), and `margins` (the supplied fitted margins). Dimensions and names
#'   match the input data. Censored entries in `Z` are cutoff placeholders,
#'   not observations. Both probability tails are evaluated on the log scale.
#' @seealso [copmi_fit_copula()]
#' @export
#' @examples
#' margins <- copmi_fit_margins(nhanes_pah, margin_mode = "normal")
#' latent <- copmi_transform(margins)
#' head(latent$Z)
#' latent$z_lod
copmi_transform <- function(object) {
  .check_class(object, "copmi_margins", "object")
  dat <- object$work_data
  Z <- matrix(NA_real_, nrow(dat$X_cens), ncol(dat$X_cens), dimnames = dimnames(dat$X_cens))
  z_lod <- rep(Inf, ncol(Z))
  for (j in seq_len(ncol(Z))) {
    fit <- object$fits[[j]]
    Z[, j] <- .gaussian_scores(dat$X_cens[, j], fit)
    if (is.finite(dat$cutoffs[j])) z_lod[j] <- .gaussian_scores(dat$cutoffs[j], fit)
  }
  if (any(!is.finite(Z)) || any(!is.finite(z_lod[colSums(dat$ind == 0L) > 0L]))) {
    stop("Marginal transformation produced invalid Gaussian scores.", call. = FALSE)
  }
  structure(list(Z = Z, ind = dat$ind, z_lod = z_lod, margins = object), class = "copmi_latent")
}

.gaussian_scores <- function(x, fit) {
  pars <- as.list(fit$parameters)
  if (fit$family == "norm") return((x - pars$mean) / pars$sd)
  cdf <- .distribution(fit$family)$p
  lp <- do.call(cdf, c(list(q = x, log.p = TRUE), pars))
  z <- stats::qnorm(lp, log.p = TRUE)
  upper <- which(lp > log(0.5))
  z[upper] <- stats::qnorm(do.call(cdf,
    c(list(q = x[upper], lower.tail = FALSE, log.p = TRUE), pars)),
    lower.tail = FALSE, log.p = TRUE)
  z
}

.inverse_scores <- function(z, fit) {
  if (any(!is.finite(z))) stop("Inverse transform requires finite latent scores.", call. = FALSE)
  pars <- as.list(fit$parameters)
  if (fit$family == "norm") return(pars$mean + pars$sd * z)
  dist <- .distribution(fit$family)
  x <- numeric(length(z))
  fallback_count <- 0L
  for (lower in c(TRUE, FALSE)) {
    idx <- which((z <= 0) == lower)
    if (!length(idx)) next
    lp <- stats::pnorm(z[idx], lower.tail = lower, log.p = TRUE)
    value <- tryCatch(do.call(dist$q, c(list(p = lp, lower.tail = lower,
      log.p = TRUE), pars)), error = function(e) rep(NaN, length(idx)))
    check <- do.call(dist$p, c(list(q = value, lower.tail = lower, log.p = TRUE), pars))
    bad <- !is.finite(value) | !is.finite(check) | abs(check - lp) > 1e-8 * pmax(1, abs(lp))
    if (dist$positive) bad <- bad | value <= 0
    for (k in which(bad)) {
      value[k] <- .solve_marginal_quantile(lp[k], fit, lower)
      fallback_count <- fallback_count + 1L
    }
    x[idx] <- value
  }
  attr(x, "quantile_fallbacks") <- fallback_count
  x
}

# Solve the same positive-support quantile on the representable log-x interval.
.solve_marginal_quantile <- function(log_p, fit, lower.tail) {
  dist <- .distribution(fit$family)
  if (!dist$positive || !is.finite(log_p) || log_p >= 0) {
    stop("No valid positive-support quantile fallback for this input.", call. = FALSE)
  }
  cdf <- dist$p
  residual <- function(log_x) do.call(cdf, c(list(q = exp(log_x),
    lower.tail = lower.tail, log.p = TRUE), as.list(fit$parameters))) - log_p
  interval <- c(log(.Machine$double.xmin * .Machine$double.eps), log(.Machine$double.xmax))
  endpoints <- vapply(interval, residual, numeric(1))
  if (anyNA(endpoints) || prod(sign(endpoints)) > 0) {
    stop("Requested quantile is outside the representable CDF range.", call. = FALSE)
  }
  root <- stats::uniroot(residual, interval, f.lower = endpoints[1], f.upper = endpoints[2],
                         tol = 1e-12, maxiter = 200L, check.conv = TRUE)$root
  x <- exp(root)
  if (!is.finite(x) || x <= 0 || !is.finite(residual(root)) ||
      abs(residual(root)) > 1e-8 * max(1, abs(log_p))) {
    stop("Inverse-CDF root did not pass the probability residual check.", call. = FALSE)
  }
  x
}

.inverse_latent <- function(Z, margins) {
  X <- margins$work_data$X_cens
  inverse_fallbacks <- 0L
  for (j in seq_len(ncol(X))) {
    idx <- which(margins$work_data$ind[, j] == 0L)
    if (!length(idx)) next
    value <- .inverse_scores(Z[idx, j], margins$fits[[j]])
    inverse_fallbacks <- inverse_fallbacks + (attr(value, "quantile_fallbacks") %||% 0L)
    X[idx, j] <- value
  }
  # Support and censoring are checked on the fitted scale before undoing the shift.
  positive <- vapply(margins$fits, function(fit) .distribution(fit$family)$positive, logical(1))
  X <- .check_completed(margins$work_data, X, positive = positive)
  boundary <- attr(X, "boundary")
  if (!is.null(margins$shift)) {
    X <- .check_completed(margins$analysis_data,
      sweep(X, 2, margins$shift$anchor, "+"))
    boundary <- c(sum(boundary[1], attr(X, "boundary")[1]),
                  max(boundary[2], attr(X, "boundary")[2]))
  }
  if (margins$input_scale == "raw") {
    X <- .check_completed(margins$input, exp(X), positive = TRUE)
    boundary <- c(sum(boundary[1], attr(X, "boundary")[1]),
                  max(boundary[2], attr(X, "boundary")[2]))
  }
  attr(X, "boundary") <- NULL
  attr(X, "numerics") <- c(quantile_fallbacks = inverse_fallbacks,
    boundary_adjustments = unname(boundary[1]), max_boundary_adjustment = unname(boundary[2]))
  X
}
