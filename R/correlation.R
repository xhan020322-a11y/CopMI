# Projection is allowed only for a finite pairwise initialization matrix.
# EM covariance updates must already be positive definite.
make_pd_cor <- function(S, project = FALSE) {
  if (!is.matrix(S) || nrow(S) != ncol(S) || any(!is.finite(S)) || any(diag(S) <= 0)) {
    stop("Correlation input must be a finite square matrix with positive diagonal.", call. = FALSE)
  }
  R <- stats::cov2cor((S + t(S)) / 2)
  pd <- tryCatch({ chol(R); TRUE }, error = function(e) FALSE)
  if (pd) return(R)
  if (!project) stop("The EM covariance update is not positive definite.", call. = FALSE)
  projected <- Matrix::nearPD(R, corr = TRUE)
  if (!projected$converged) stop("Initial correlation projection did not converge.", call. = FALSE)
  R <- as.matrix(projected$mat)
  attr(R, "pd_adjustment") <- projected$normF
  R
}

lyles_pair_negloglik <- function(par, z1, z2, obs1, obs2, lod1, lod2, rho_max = 1) {
  mu1 <- par[1]
  sd1 <- exp(par[2])
  mu2 <- par[3]
  sd2 <- exp(par[4])
  rho <- rho_max * tanh(par[5])
  if (any(!is.finite(c(mu1, sd1, mu2, sd2, rho))) || sd1 <= 0 || sd2 <= 0 || abs(rho) >= 1) return(1e100)
  Sigma <- matrix(c(sd1^2, rho * sd1 * sd2, rho * sd1 * sd2, sd2^2), 2, 2)
  ll <- 0
  oo <- obs1 & obs2
  if (any(oo)) ll <- ll + sum(mvtnorm::dmvnorm(cbind(z1[oo], z2[oo]),
    mean = c(mu1, mu2), sigma = Sigma, log = TRUE))
  co <- !obs1 & obs2
  if (any(co)) {
    mu <- mu1 + rho * sd1 / sd2 * (z2[co] - mu2)
    ll <- ll + sum(stats::dnorm(z2[co], mu2, sd2, log = TRUE)) +
      sum(stats::pnorm(lod1, mu, sd1 * sqrt(1 - rho^2), log.p = TRUE))
  }
  oc <- obs1 & !obs2
  if (any(oc)) {
    mu <- mu2 + rho * sd2 / sd1 * (z1[oc] - mu1)
    ll <- ll + sum(stats::dnorm(z1[oc], mu1, sd1, log = TRUE)) +
      sum(stats::pnorm(lod2, mu, sd2 * sqrt(1 - rho^2), log.p = TRUE))
  }
  cc <- sum(!obs1 & !obs2)
  if (cc > 0L) ll <- ll + cc * log(as.numeric(mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf), upper = c(lod1, lod2), mean = c(mu1, mu2), sigma = Sigma)))
  if (is.finite(ll)) -ll else 1e100
}

estimate_lyles_pair_cor <- function(z1, z2, ind1, ind2, lod1, lod2,
                                    min_obs = 2L, maxit = 2000L, rho_max = 1) {
  obs1 <- ind1 == 1L
  obs2 <- ind2 == 1L
  if (sum(obs1) < min_obs || sum(obs2) < min_obs) {
    stop("Too few observed values for pairwise initialization; check min_obs.", call. = FALSE)
  }
  both <- obs1 & obs2
  # Zero is only an optimizer starting value when the observed pairs cannot
  # supply a correlation; it is never substituted for a failed fitted pair.
  r0 <- if (sum(both) >= 2L && stats::sd(z1[both]) > 0 && stats::sd(z2[both]) > 0)
    stats::cor(z1[both], z2[both]) else 0
  if (all(obs1) && all(obs2)) return(list(rho = r0, optimizer = NA_character_,
    attempts = data.frame(optimizer = character(), converged = logical(), message = character(), warnings = character())))
  start <- c(mean(z1[obs1]), log(stats::sd(z1[obs1])), mean(z2[obs2]),
             log(stats::sd(z2[obs2])), atanh(0.95 * r0 / max(rho_max, abs(r0))))
  objective <- function(par) lyles_pair_negloglik(par, z1, z2, obs1, obs2, lod1, lod2, rho_max)
  attempts <- list()
  for (method in c("L-BFGS-B", "Nelder-Mead")) {
    fit <- .run_optimizer(start, objective, method, maxit)
    attempts[[method]] <- data.frame(optimizer = method, converged = fit$valid,
                                     message = fit$message, warnings = fit$warnings)
    if (fit$valid) return(list(rho = rho_max * tanh(fit$par[5]), optimizer = method,
                               attempts = do.call(rbind, attempts)))
  }
  stop("Both optimizers failed during pairwise initialization: ",
       paste(vapply(attempts, function(a) a$message, character(1)), collapse = " | "), call. = FALSE)
}

estimate_lyles_cor_init <- function(Z_full, IND, z_lod, min_obs = 2L,
                                    maxit = 2000L, rho_max = 1, verbose = FALSE) {
  p <- ncol(Z_full)
  R <- diag(1, p)
  pairs <- attempts <- list()
  for (a in seq_len(p)) {
    for (b in seq_len(p)) {
      if (b <= a) next
      pair <- tryCatch(estimate_lyles_pair_cor(Z_full[, a], Z_full[, b], IND[, a], IND[, b],
        z_lod[a], z_lod[b], min_obs, maxit, rho_max),
        error = function(e) stop("Initial correlation, columns ", a, " and ", b,
          ": ", conditionMessage(e), call. = FALSE))
      R[a, b] <- R[b, a] <- pair$rho
      key <- length(pairs) + 1L
      pairs[[key]] <- data.frame(var1 = a, var2 = b, rho_lyles = pair$rho,
        n_obs_pair = sum(IND[, a] == 1L & IND[, b] == 1L), optimizer = pair$optimizer,
        optimizer_retries = max(nrow(pair$attempts) - 1L, 0L))
      attempts[[key]] <- data.frame(var1 = rep(a, nrow(pair$attempts)),
                                    var2 = rep(b, nrow(pair$attempts)), pair$attempts)
      if (verbose) message("Initial correlation (", a, ", ", b, "): ", signif(pair$rho, 4))
    }
  }
  R <- make_pd_cor(R, project = TRUE)
  adjustment <- attr(R, "pd_adjustment") %||% 0
  attr(R, "pd_adjustment") <- NULL
  pair_table <- if (length(pairs)) do.call(rbind, pairs) else data.frame()
  list(R = R, pair_results = pair_table,
    optimization = if (length(attempts)) do.call(rbind, attempts) else data.frame(),
    pd_adjustment = adjustment,
    fallback_records = combine_fallbacks(
      fallback_record("Copula_EM", "initialization_optimizer_retry", sum(pair_table$optimizer_retries),
        "Retried the same pairwise censored-normal likelihood with Nelder-Mead."),
      fallback_record("Copula_EM", "initialization_pd_projection", as.integer(adjustment > 0),
        paste0("Projected the pairwise initial matrix; Frobenius adjustment = ", signif(adjustment, 5), "."))))
}
