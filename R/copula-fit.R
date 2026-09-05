#' Estimate the latent Gaussian copula correlation
#'
#' Estimate pairwise censored-normal correlations for initialization, then
#' update the latent correlation matrix by EM using truncated-normal moments.
#' This step estimates the correlation model without drawing completed data.
#' @param object A `copmi_latent` object from [copmi_transform()].
#' @param max_iter Nonnegative integer maximum number of EM iterations. Zero
#'   returns initialization only and sets `converged = FALSE`.
#' @param tol Positive finite convergence tolerance: the maximum absolute
#'   elementwise change in the correlation matrix must be smaller than `tol`.
#' @param seed Nonnegative integer seed, or `NULL` to use and advance the current
#'   RNG stream. Explicit seeds restore the caller's RNG state on exit.
#' @param verbose Whether to print initialization and EM progress.
#' @param lyles_control Named list with optional `min_obs` (default 2), `maxit`
#'   (2000 per optimizer), and `rho_max` (1) for pairwise censored-normal initialization.
#'   Pairwise likelihoods try L-BFGS-B then Nelder-Mead. Failed pairs
#'   raise errors. `rho_max < 1` explicitly restricts censored-pair fits;
#'   fully observed pairs use their sample correlation.
#' @return A `copmi_copula` list containing `Sigma_hat`, `Sigma_init` (correlation
#'   matrices), `converged`, `n_iter`, `em_change_history`, `init_diagnostics`
#'   (including the per-pair table), `fallback_records`, `latent` (the supplied
#'   object), `rng_state` (state after estimation for reproducible draws),
#'   `seed`, `max_iter`, `tol`, and `lyles_control`.
#' @details Convergence refers only to the correlation-change criterion, not
#'   convergence of the marginal optimizers or proof of a global optimum.
#'   An indefinite finite pairwise initial matrix is projected with
#'   [Matrix::nearPD()]; its adjustment is recorded as `init_diagnostics$pd_adjustment`.
#'   EM updates are not projected. Moment or linear-algebra failures raise errors
#'   rather than replacing moments or adding a ridge. Unconverged fits can be
#'   inspected but require explicit opt-in before imputation.
#' @seealso [copmi_impute()], [diagnostics()]
#' @export
#' @examples
#' # Use the first 100 rows for a quick help example; full data have 1330 rows.
#' dat <- make_lod_data(nhanes_pah$X_cens[1:100, ],
#'                      nhanes_pah$ind[1:100, ], nhanes_pah$cutoffs, scale = "log")
#' margins <- copmi_fit_margins(dat, margin_mode = "normal")
#' model <- copmi_fit_copula(copmi_transform(margins), max_iter = 100, seed = 1)
#' model$Sigma_hat
#' model$em_change_history
copmi_fit_copula <- function(object, max_iter = 100, tol = 1e-4, seed = 1,
                             verbose = FALSE, lyles_control = list()) {
  .check_class(object, "copmi_latent", "object")
  .check_number(max_iter, "max_iter", 0, .Machine$integer.max, TRUE)
  .check_number(tol, "tol", .Machine$double.eps)
  .check_flag(verbose, "verbose")
  lyles_control <- .check_lyles_control(lyles_control)
  res <- .with_rng(seed, code = .fit_copula_engine(
    object$Z, object$ind, object$z_lod, max_iter, tol, verbose, lyles_control))
  dimnames(res$Sigma_hat) <- dimnames(res$Sigma_init) <-
    list(colnames(object$Z), colnames(object$Z))
  res$latent <- object
  res$seed <- seed
  res$max_iter <- max_iter
  res$tol <- tol
  res$lyles_control <- lyles_control
  structure(res, class = "copmi_copula")
}
