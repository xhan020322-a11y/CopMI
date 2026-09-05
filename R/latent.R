#' Transform fitted margins to latent Gaussian scores
#'
#' Apply `qnorm(F_j(x))` using the selected fitted CDF for each variable.
#' @param object A `copmi_margins` object from [copmi_fit_margins()].
#' @return A `copmi_latent` list with `Z` (numeric matrix of Gaussian scores),
#'   `ind` (observation indicators), `z_lod` (latent cutoffs; `Inf` for absent
#'   cutoffs), and `margins` (the supplied fitted margins). Dimensions and names
#'   match the input data. Censored entries in `Z` are cutoff placeholders,
#'   not observations. Both probability tails are evaluated on the log scale; no probability clipping is used.
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
  pars <- as.list(fit$parameters)
  if (fit$family == "norm") return(pars$mean + pars$sd * z)
  quantile <- .distribution(fit$family)$q
  x <- do.call(quantile, c(list(p = stats::pnorm(z, log.p = TRUE), log.p = TRUE), pars))
  upper <- which(z > 0)
  x[upper] <- do.call(quantile, c(list(p = stats::pnorm(z[upper], lower.tail = FALSE,
    log.p = TRUE), lower.tail = FALSE, log.p = TRUE), pars))
  x
}

.inverse_latent <- function(Z, margins) {
  X <- matrix(NA_real_, nrow(Z), ncol(Z))
  for (j in seq_len(ncol(X))) X[, j] <- .inverse_scores(Z[, j], margins$fits[[j]])
  if (!is.null(margins$shift)) X <- sweep(X, 2, margins$shift$anchor, "+")
  X <- .check_completed(margins$analysis_data, X)
  if (margins$input_scale == "raw") X <- exp(X)
  .check_completed(margins$input, X, positive = margins$input_scale == "raw")
}
