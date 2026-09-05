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
  for (iter in seq_len(max_iter)) {
    S <- matrix(0, p, p)
    for (i in seq_len(n)) {
      missing <- which(IND[i, ] == 0L)
      z <- Z_full[i, ]
      if (!length(missing)) {
        S <- S + tcrossprod(z)
        next
      }
      conditional <- .conditional_normal(z, which(IND[i, ] == 1L), missing, Sigma)
      moments <- tryCatch(.truncated_moments(mean = conditional$mean,
        sigma = conditional$sigma, lower = rep(-Inf, length(missing)),
        upper = z_lod[missing], doComputeVariance = TRUE),
        error = function(e) stop("EM iteration ", iter, ", row ", i,
          ": ", conditionMessage(e), call. = FALSE))
      z[missing] <- moments$tmean
      EZZ <- tcrossprod(z)
      EZZ[missing, missing] <- EZZ[missing, missing, drop = FALSE] + moments$tvar
      S <- S + EZZ
    }
    updated <- make_pd_cor(S / n)
    history[iter] <- max(abs(updated - Sigma))
    Sigma <- updated
    if (verbose) message("EM iteration ", iter, ": change = ", signif(history[iter], 4))
    if (history[iter] < tol) {
      converged <- TRUE
      break
    }
  }
  list(rng_state = get(".Random.seed", envir = .GlobalEnv), Sigma_hat = Sigma,
    Sigma_init = init$R, init_diagnostics = init, em_change_history = history,
    fallback_records = init$fallback_records, converged = converged, n_iter = length(history))
}
