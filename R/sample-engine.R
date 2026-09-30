# Prepare conditional Gaussian matrices once for a given censoring pattern.
.prepare_conditional_normal <- function(observed, missing, Sigma) {
  Sigma <- .checked_covariance(Sigma)
  order <- c(observed, missing)
  factor <- chol(Sigma[order, order, drop = FALSE])
  q <- length(missing)
  k <- length(observed)
  block <- k + seq_len(q)
  covariance <- crossprod(factor[block, block, drop = FALSE])
  coefficients <- NULL
  if (k) {
    coefficients <- t(backsolve(factor[seq_len(k), seq_len(k), drop = FALSE],
                               factor[seq_len(k), block, drop = FALSE]))
  }
  list(observed = observed, n_missing = q, coefficients = coefficients,
       sigma = covariance)
}

# Means still depend on each row's observed values.
.conditional_normal_mean <- function(z, parameters) {
  observed <- parameters$observed
  mean <- if (is.matrix(z)) matrix(0, nrow(z), parameters$n_missing) else
    rep(0, parameters$n_missing)
  if (length(observed)) {
    mean <- if (is.matrix(z)) z[, observed, drop = FALSE] %*% t(parameters$coefficients) else
      as.vector(parameters$coefficients %*% z[observed])
  }
  if (any(!is.finite(mean))) stop("Conditional mean is not finite.", call. = FALSE)
  mean
}

# Conditional Gaussian algebra is shared by estimation and sampling.
.conditional_normal <- function(z, observed, missing, Sigma) {
  parameters <- .prepare_conditional_normal(observed, missing, Sigma)
  list(mean = .conditional_normal_mean(z, parameters), sigma = parameters$sigma)
}

# Draw N(mean, sd^2) below upper using log probabilities in the left tail.
.rnorm_upper_trunc <- function(mean, sd, upper) {
  if (any(!is.finite(mean)) || any(!is.finite(sd) | sd <= 0) ||
      any(!is.finite(upper))) {
    stop("Invalid truncated-normal mean, SD, or upper bound.", call. = FALSE)
  }
  log_p <- stats::pnorm((upper - mean) / sd, log.p = TRUE)
  out <- mean + sd * stats::qnorm(log_p + log(stats::runif(length(mean))), log.p = TRUE)
  if (any(!is.finite(out))) {
    stop("Truncated-normal draw must be finite and strictly below its upper bound.", call. = FALSE)
  }
  checked <- .round_upper_bound(out, upper, scale = pmax(abs(mean), sd, abs(upper)))
  tail_scale <- sd / pmax(1, (mean - upper) / sd)
  if (any(abs(checked - out) > 1e-6 * tail_scale)) {
    stop("Floating-point resolution is too coarse for this truncated draw.", call. = FALSE)
  }
  checked
}

# Rows with the same censoring pattern share a conditional precision matrix.
.prepare_gibbs_groups <- function(Z_full, IND, z_lod, Sigma_hat) {
  rows <- which(rowSums(IND == 0L) > 0L)
  if (!length(rows)) return(list())
  patterns <- apply(IND[rows, , drop = FALSE] == 0L, 1L, paste0, collapse = "")
  row_groups <- split(rows, patterns, drop = TRUE)
  lapply(names(row_groups), function(pattern) {
    rows <- row_groups[[pattern]]
    missing <- which(IND[rows[1L], ] == 0L)
    observed <- which(IND[rows[1L], ] == 1L)
    tryCatch({
      conditional <- .conditional_normal(Z_full[rows, , drop = FALSE],
                                          observed, missing, Sigma_hat)
      precision <- chol2inv(chol(conditional$sigma))
      if (any(!is.finite(precision)) || any(diag(precision) <= 0)) {
        stop("Conditional precision must be finite and positive definite.", call. = FALSE)
      }
      state <- matrix(NA_real_, length(rows), length(missing))
      boundary <- c(count = 0, max_adjustment = 0)
      for (k in seq_along(missing)) {
        draw <- .rnorm_upper_trunc(conditional$mean[, k],
          sqrt(conditional$sigma[k, k]), z_lod[missing[k]])
        event <- attr(draw, "boundary")
        boundary <- c(boundary[1] + event[1], max(boundary[2], event[2]))
        state[, k] <- draw
      }
      list(pattern = pattern, rows = rows, missing = missing, upper = z_lod[missing],
           mean = conditional$mean, precision = precision, state = state,
           boundary = unname(boundary))
    }, error = function(e) stop("Gibbs initialization, pattern ", pattern,
                                ": ", conditionMessage(e), call. = FALSE))
  })
}

# Each sweep updates every censored coordinate; states persist between draws.
.draw_copula_engine <- function(Z_full, IND, z_lod, Sigma_hat, M,
                                gibbs_burn, gibbs_thin) {
  groups <- .prepare_gibbs_groups(Z_full, IND, z_lod, Sigma_hat)
  kept <- if (length(groups)) as.integer(gibbs_burn + seq_len(M) * as.double(gibbs_thin)) else integer()
  total <- if (length(kept)) max(kept) else 0L
  n_censored <- sum(IND == 0L)
  trace_mean <- trace_update <- numeric(total)
  imputations <- replicate(M, Z_full, simplify = FALSE)
  keep_id <- 1L
  boundary_count <- sum(vapply(groups, function(g) g$boundary[1], numeric(1)))
  boundary_max <- max(c(0, vapply(groups, function(g) g$boundary[2], numeric(1))))
  for (sweep in seq_len(total)) {
    state_sum <- update_sum <- 0
    for (g in seq_along(groups)) {
      gr <- groups[[g]]
      centered <- gr$state - gr$mean
      for (k in seq_along(gr$missing)) {
        other <- setdiff(seq_along(gr$missing), k)
        qkk <- gr$precision[k, k]
        mean <- gr$mean[, k]
        if (length(other)) mean <- mean -
          as.numeric(centered[, other, drop = FALSE] %*% gr$precision[k, other]) / qkk
        value <- tryCatch(.rnorm_upper_trunc(mean, sqrt(1 / qkk), gr$upper[k]),
          error = function(e) stop("Gibbs sweep ", sweep, ", pattern ", gr$pattern,
            ", variable ", gr$missing[k], ": ", conditionMessage(e), call. = FALSE))
        event <- attr(value, "boundary")
        boundary_count <- boundary_count + event[1]
        boundary_max <- max(boundary_max, event[2])
        update_sum <- update_sum + sum(abs(value - gr$state[, k]))
        gr$state[, k] <- value
        centered[, k] <- value - gr$mean[, k]
      }
      groups[[g]] <- gr
      state_sum <- state_sum + sum(gr$state)
    }
    trace_mean[sweep] <- state_sum / n_censored
    trace_update[sweep] <- update_sum / n_censored
    if (sweep == kept[keep_id]) {
      for (gr in groups) imputations[[keep_id]][gr$rows, gr$missing] <- gr$state
      keep_id <- keep_id + 1L
    }
  }
  list(Z_imp_list = imputations, sampling = list(
    scan = "systematic", burn_in = as.integer(gibbs_burn), thin = as.integer(gibbs_thin),
    completed = TRUE, total_sweeps = total, kept_sweeps = kept,
    n_patterns = length(groups), n_rows_with_censoring = sum(rowSums(IND == 0L) > 0L),
    n_censored_latent = n_censored, n_scalar_updates = as.double(total) * n_censored,
    boundary_adjustments = unname(boundary_count), max_boundary_adjustment = unname(boundary_max),
    trace_mean_censored_z = trace_mean, trace_mean_abs_update = trace_update))
}
