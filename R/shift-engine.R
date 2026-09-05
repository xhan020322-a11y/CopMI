estimate_observed_log_sd <- function(X_log_cens, IND) {
  sd_hat <- vapply(seq_len(ncol(X_log_cens)), function(j)
    stats::sd(X_log_cens[IND[, j] == 1L, j]), numeric(1))
  if (any(!is.finite(sd_hat) | sd_hat <= 0)) {
    stop("SD-shift requires at least two varying observed values in every column.", call. = FALSE)
  }
  list(sd = sd_hat, source = rep("observed_sd", length(sd_hat)))
}

make_log_sd_shift_dat <- function(dat, shift_k = 3) {
  X_log_cens <- as.matrix(dat$X_cens)
  IND <- as.matrix(dat$ind)
  LOD_log <- as.numeric(dat$cutoffs)
  p <- ncol(X_log_cens)
  
  sd_est <- estimate_observed_log_sd(X_log_cens, IND)
  anchor <- rep(NA_real_, p)
  min_shifted_observed_or_cutoff <- rep(NA_real_, p)
  for (j in seq_len(p)) {
    obs_j <- X_log_cens[IND[, j] == 1L, j]
    ref_j <- if (is.finite(LOD_log[j])) LOD_log[j] else min(obs_j)
    if (!is.finite(ref_j)) {
      stop("Cannot define SD-shift anchor for variable ", colnames(X_log_cens)[j] %||% j, ".")
    }
    
    anchor_j <- ref_j - shift_k * sd_est$sd[j]
    shifted_visible <- X_log_cens[, j] - anchor_j
    check_values <- shifted_visible
    lod_shift_j <- if (is.finite(LOD_log[j])) LOD_log[j] - anchor_j else NA_real_
    if (is.finite(lod_shift_j)) {
      check_values <- c(check_values, lod_shift_j)
    }
    min_shift_j <- if (length(check_values)) min(check_values, na.rm = TRUE) else NA_real_
    
    if (any(!is.finite(check_values)) || !is.finite(min_shift_j) || min_shift_j <= 0) {
      stop(
        "Fixed SD-shift does not make visible values/cutoff positive for variable ",
        colnames(X_log_cens)[j] %||% j,
        ": shift_k=", shift_k,
        ", observed_sd=", sd_est$sd[j],
        ", min shifted visible/cutoff=", min_shift_j
      )
    }
    
    anchor[j] <- anchor_j
    min_shifted_observed_or_cutoff[j] <- min_shift_j
  }
  
  X_shift_cens <- sweep(X_log_cens, 2, anchor, "-")
  LOD_shift <- ifelse(is.finite(LOD_log), LOD_log - anchor, NA_real_)
  colnames(X_shift_cens) <- colnames(X_log_cens)
  
  shift_dat <- dat
  shift_dat$X_cens <- X_shift_cens
  shift_dat$cutoffs <- as.numeric(LOD_shift)
  shift_dat$sd_shift <- list(
    anchor = anchor,
    sd = sd_est$sd,
    sd_source = sd_est$source,
    shift_multiplier = shift_k,
    min_shifted_observed_or_cutoff = min_shifted_observed_or_cutoff,
    original_cutoffs = LOD_log,
    shifted_cutoffs = as.numeric(LOD_shift)
  )
  
  list(
    dat = shift_dat,
    anchor = anchor,
    sd = sd_est$sd,
    sd_source = sd_est$source,
    shift_multiplier = shift_k,
    min_shifted_observed_or_cutoff = min_shifted_observed_or_cutoff,
    original_cutoffs = LOD_log,
    shifted_cutoffs = as.numeric(LOD_shift)
  )
}
