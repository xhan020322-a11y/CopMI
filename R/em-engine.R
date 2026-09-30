# Compute the conditional truncated-normal moments used by the E-step.
.truncated_moments <- function(mean, sigma, lower, upper) {
  sigma <- .checked_covariance(sigma)
  q <- length(mean)
  if (nrow(sigma) != q || length(lower) != q || length(upper) != q ||
      any(!is.finite(mean)) || anyNA(lower) || anyNA(upper) || any(lower >= upper)) {
    stop("Invalid truncated-normal parameters or bounds.", call. = FALSE)
  }
  out <- tmvtnorm::mtmvnorm(mean = mean, sigma = sigma, lower = lower,
    upper = upper, doComputeVariance = TRUE)
  if (length(out$tmean) != q || !identical(dim(out$tvar), c(q, q)) ||
      any(!is.finite(out$tmean)) || any(!is.finite(out$tvar)) ||
      any(out$tmean < lower | out$tmean > upper) || any(diag(out$tvar) < 0)) {
    stop("Truncated-normal moments are invalid.", call. = FALSE)
  }
  # Moment entries use numerical integration, not only floating-point arithmetic.
  moment_scale <- max(max(abs(out$tvar)), .Machine$double.xmin)
  out$asymmetry <- max(abs(out$tvar - t(out$tvar))) / moment_scale
  # The symmetric part preserves every quadratic form.
  out$tvar <- (out$tvar + t(out$tvar)) / 2
  tolerance <- 64 * q * .Machine$double.eps * max(abs(out$tvar))
  if (min(eigen(out$tvar, symmetric = TRUE, only.values = TRUE)$values) < -tolerance) {
    stop("Truncated-normal covariance is not positive semidefinite.", call. = FALSE)
  }
  out
}

# Every E-step uses the requested truncated moments. Numerical failures stop
# with row/iteration context; an unconstrained moment is never substituted.
.fit_copula_engine <- function(Z_full, IND, z_lod, max_iter, tol, verbose, lyles_control) {
  n <- nrow(Z_full)
  p <- ncol(Z_full)
  init <- do.call(estimate_lyles_cor_init,
    c(list(Z_full = Z_full, IND = IND, z_lod = z_lod, verbose = verbose), lyles_control))
  Sigma <- init$R
  converged <- FALSE
  history <- numeric()
  stable <- 0L
  moment_asymmetry <- numeric()
  for (iter in seq_len(max_iter)) {
    S <- matrix(0, p, p)
    moment_asymmetry[iter] <- 0
    # Sigma is fixed within this iteration; discard its cache before the next.
    conditional_cache <- new.env(hash = TRUE, parent = emptyenv())
    for (i in seq_len(n)) {
      missing <- which(IND[i, ] == 0L)
      z <- Z_full[i, ]
      if (!length(missing)) {
        S <- S + tcrossprod(z)
        next
      }
      pattern <- paste(missing, collapse = ",")
      parameters <- conditional_cache[[pattern]]
      if (is.null(parameters)) {
        parameters <- .prepare_conditional_normal(which(IND[i, ] == 1L), missing, Sigma)
        conditional_cache[[pattern]] <- parameters
      }
      conditional_mean <- .conditional_normal_mean(z, parameters)
      moments <- tryCatch(.truncated_moments(mean = conditional_mean,
        sigma = parameters$sigma, lower = rep(-Inf, length(missing)),
        upper = z_lod[missing]),
        error = function(e) stop("EM iteration ", iter, ", row ", i,
          ": ", conditionMessage(e), call. = FALSE))
      moment_asymmetry[iter] <- max(moment_asymmetry[iter], moments$asymmetry)
      z[missing] <- moments$tmean
      EZZ <- tcrossprod(z)
      EZZ[missing, missing] <- EZZ[missing, missing, drop = FALSE] + moments$tvar
      S <- S + EZZ
    }
    updated <- make_pd_cor(S / n)
    history[iter] <- max(abs(updated - Sigma))
    Sigma <- updated
    if (verbose) message("EM iteration ", iter, ": change = ", signif(history[iter], 4))
    stable <- if (history[iter] < tol) stable + 1L else 0L
    if (stable >= 3L) {
      converged <- TRUE
      break
    }
  }
  list(rng_state = get(".Random.seed", envir = .GlobalEnv), Sigma_hat = Sigma,
    Sigma_init = init$R, init_diagnostics = init, em_change_history = history,
    moment_asymmetry_history = moment_asymmetry,
    fallback_records = init$fallback_records, converged = converged, n_iter = length(history))
}
