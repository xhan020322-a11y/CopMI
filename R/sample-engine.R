# Conditional Gaussian algebra is shared by estimation and sampling.
.conditional_normal <- function(z, observed, missing, Sigma) {
  covariance <- Sigma[missing, missing, drop = FALSE]
  mean <- rep(0, length(missing))
  if (length(observed)) {
    cross <- Sigma[missing, observed, drop = FALSE]
    coefficients <- t(solve(Sigma[observed, observed, drop = FALSE], t(cross)))
    mean <- as.vector(coefficients %*% z[observed])
    covariance <- covariance - coefficients %*% t(cross)
  }
  list(mean = mean, sigma = (covariance + t(covariance)) / 2)
}

.draw_copula_engine <- function(Z_full, IND, z_lod, Sigma_hat, M, ag, thinning) {
  imputations <- replicate(M, Z_full, simplify = FALSE)
  for (i in seq_len(nrow(Z_full))) {
    missing <- which(IND[i, ] == 0L)
    if (!length(missing)) next
    conditional <- .conditional_normal(Z_full[i, ], which(IND[i, ] == 1L), missing, Sigma_hat)
    upper <- z_lod[missing]
    start <- pmin(conditional$mean, upper) - 1
    for (m in seq_len(M)) {
      draw <- tryCatch(.truncated_draw(n = 1, mean = conditional$mean,
        sigma = conditional$sigma, lower = rep(-Inf, length(missing)), upper = upper,
        algorithm = ag, start.value = start, thinning = thinning),
        error = function(e) stop("Imputation ", m, ", row ", i,
          ": ", conditionMessage(e), call. = FALSE))
      imputations[[m]][i, missing] <- as.numeric(draw)
    }
  }
  list(Z_imp_list = imputations)
}

.truncated_moments <- function(...) {
  out <- tmvtnorm::mtmvnorm(...)
  if (any(!is.finite(out$tmean)) || any(!is.finite(out$tvar)) || any(diag(out$tvar) < 0)) {
    stop("Truncated-normal moments are invalid.", call. = FALSE)
  }
  out
}

.truncated_draw <- function(...) {
  out <- tmvtnorm::rtmvnorm(...)
  if (any(!is.finite(out))) stop("Truncated-normal draw is nonfinite.", call. = FALSE)
  out
}
