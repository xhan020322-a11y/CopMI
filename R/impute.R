#' Draw completed data from a fitted copula model
#'
#' Sample censored latent values conditional on observed values, invert the
#' fitted marginal transforms, and return completed data on the input scale.
#' @param object A `copmi_copula` object from [copmi_fit_copula()].
#' @param m Positive integer number of completed data sets.
#' @param seed Nonnegative integer seed. The default `NULL` continues the RNG
#'   state saved by the fitted model, while restoring the caller's RNG state.
#'   Repeated calls with the same fit and seed give the same result; use a new
#'   seed for additional independent runs of the sampler.
#' @param ag Sampling algorithm: `"gibbs"`, `"gibbsR"`, or `"rejection"`, passed
#'   to [tmvtnorm::rtmvnorm()].
#' @param allow_unconverged Allow exploratory draws from an unconverged or
#'   initialization-only fit. Defaults to `FALSE`; `converged` remains `FALSE`.
#' @param thinning Positive integer number of Gibbs steps per retained draw.
#'   Ignored by the rejection sampler.
#' @return A `copmi_mi` object with `imp_list` (list of `m` numeric matrices),
#'   `M`, `method`, `is_mi`, `input`, `input_scale`, `margin_mode`, `converged`,
#'   `n_iter`, `Sigma_hat`, `family_selected`, `bic_table`, `fallback_records`,
#'   `shift_anchor`, `model`, `call`, and `extra`. `extra` retains latent draws,
#'   selected margin parameters, optimization diagnostics, initialization,
#'   and EM history. Use [complete()] and [diagnostics()] for stable extraction.
#' @details Observed values and dimension names are preserved. Censored values
#'   must be finite and strictly below their cutoff. Moment, sampling, or inverse
#'   transformation failures raise errors; no fixed-value replacement is used. Marginal parameters and
#'   correlations are held fixed across draws; parameter and family-selection
#'   uncertainty is not sampled. The Gibbs sampler starts inside the
#'   truncation region; users should assess mixing and select appropriate thinning.
#' @seealso [copula_em_impute()], [complete()]
#' @export
#' @examples
#' # Use the first 100 rows for a quick help example; full data have 1330 rows.
#' dat <- make_lod_data(nhanes_pah$X_cens[1:100, ],
#'                      nhanes_pah$ind[1:100, ], nhanes_pah$cutoffs, scale = "log")
#' margins <- copmi_fit_margins(dat, margin_mode = "normal")
#' model <- copmi_fit_copula(copmi_transform(margins), max_iter = 100, seed = 1)
#' result <- copmi_impute(model, m = 2)
#' dim(complete(result, 1))
#' summary(result)
copmi_impute <- function(object, m = 5, seed = NULL, ag = "gibbs", thinning = 2,
                         allow_unconverged = FALSE) {
  .check_class(object, "copmi_copula", "object")
  .check_number(m, "m", 1, .Machine$integer.max, TRUE)
  .check_number(thinning, "thinning", 1, .Machine$integer.max, TRUE)
  ag <- match.arg(ag, c("gibbs", "gibbsR", "rejection"))
  .check_flag(allow_unconverged, "allow_unconverged")
  if (!object$converged && !allow_unconverged) {
    stop("Copula EM has not converged. Refit with more iterations; use allow_unconverged = TRUE only for exploratory draws.", call. = FALSE)
  }
  latent <- object$latent
  margins <- latent$margins
  drawn <- .with_rng(seed, state = if (is.null(seed)) object$rng_state else NULL,
    code = .draw_copula_engine(latent$Z, latent$ind, latent$z_lod,
                               object$Sigma_hat, m, ag, thinning))
  imp_list <- lapply(drawn$Z_imp_list, .inverse_latent, margins = margins)
  ct <- margins$candidate_table
  ot <- margins$optimization
  retries <- if (nrow(ot)) sum(duplicated(ot[, c("variable", "candidate", "start")])) else 0L
  fallbacks <- combine_fallbacks(
    object$fallback_records,
    fallback_record("Copula_EM", "margin_optimizer_retry", retries,
      "Retried the same censored marginal likelihood with the second optimizer."),
    fallback_record("Copula_EM", "margin_family_fit_failed", sum(!ct$converged),
      summarize_messages(ct$candidate[!ct$converged]))
  )
  bic_table <- lapply(margins$fits, `[[`, "candidate_table")
  shift <- margins$shift
  extra <- list(Z_imp_list = drawn$Z_imp_list, Sigma_hat = object$Sigma_hat,
    Sigma_init = object$Sigma_init, init_diagnostics = object$init_diagnostics,
    em_change_history = object$em_change_history, fallback_records = fallbacks,
    n_iter = object$n_iter, converged = object$converged,
    family_hat = margins$family_selected, family_selected = margins$family_selected,
    margin_fits = margins$fits, fit_how = vapply(margins$fits, `[[`, character(1), "how"),
    bic_table = bic_table, candidate_table = ct, optimization = ot, z_lod = latent$z_lod,
    shift_anchor = shift$anchor, shift_sd = shift$sd, shift_sd_source = shift$sd_source,
    shift_multiplier = shift$shift_multiplier, shifted_cutoffs = shift$shifted_cutoffs,
    original_cutoffs = margins$analysis_data$cutoffs, scale = margins$input_scale,
    sampling = list(ag = ag, thinning = thinning, seed = seed))
  structure(list(method = if (margins$margin_mode == "normal") "Copula_EM_normal" else "Copula_EM",
    imp_list = imp_list, M = as.integer(m), is_mi = TRUE, extra = extra,
    call = match.call(), input = margins$input, input_scale = margins$input_scale,
    margin_mode = margins$margin_mode, converged = object$converged, n_iter = object$n_iter,
    Sigma_hat = object$Sigma_hat, family_selected = margins$family_selected,
    bic_table = bic_table, fallback_records = fallbacks, shift_anchor = shift$anchor,
    model = object), class = "copmi_mi")
}

#' Multiple imputation for left-censored exposure data
#'
#' Run the complete pipeline: validate data, optionally log the input, fit
#' margins, transform to Gaussian scores, estimate copula correlation by EM,
#' and draw completed data sets on the supplied scale.
#' @param x A `copmi_lod_data` object, or a numeric matrix/data frame.
#' @param ind Observation indicator matrix, required for matrix/data frame `x`.
#' @param cutoffs Column cutoffs, required for matrix/data frame `x`.
#' @inheritParams copmi_fit_margins
#' @inheritParams copmi_fit_copula
#' @inheritParams copmi_impute
#' @param ... Reserved; unused arguments raise an error.
#' @return A `copmi_mi` object as documented in [copmi_impute()]. Each completed
#'   matrix has the original dimensions, names, observed values, and input
#'   scale. `converged` records EM convergence, not marginal optimizer status.
#' @details Use `input_scale = "raw", margin_mode = "normal"` when positive raw
#'   exposures are assumed log-normal; the package applies natural log and fits
#'   censored normal margins. Use `input_scale = "log"` for data already logged.
#'   Without the normal assumption, use `margin_mode = "sd_shift"` (or
#'   `"select"`) for BIC selection on the shifted log scale. No normality test
#'   silently decides the mode. Both modes use the ordered `optim_methods`.
#'
#'   Separate calls to [copmi_fit_margins()], [copmi_transform()],
#'   [copmi_fit_copula()] with the same seed, and [copmi_impute()] with its
#'   default seed reproduce the one-call result. `max_iter = 0` is available
#'   for initialization-only diagnostics, not a converged EM fit.
#' @seealso [make_lod_data()], [complete()], [diagnostics()]
#' @export
#' @examples
#' # Use the first 100 rows for a quick help example; full data have 1330 rows.
#' dat <- make_lod_data(nhanes_pah$X_cens[1:100, ],
#'                      nhanes_pah$ind[1:100, ], nhanes_pah$cutoffs, scale = "log")
#' # Input is already on the log scale; assume normal log margins.
#' fit <- copula_em_impute(dat, margin_mode = "normal",
#'                          m = 2, max_iter = 100, seed = 1)
#' complete(fit, 1)
#' diagnostics(fit)$optimization
#' # Positive raw-scale data: transform internally and return raw-scale values.
#' raw <- make_lod_data(exp(dat$X_cens), dat$ind,
#'                      exp(dat$cutoffs))
#' fit_raw <- copula_em_impute(raw, input_scale = "raw", margin_mode = "normal",
#'                              m = 1, max_iter = 100, seed = 1)
#' head(complete(fit_raw, 1))
copula_em_impute <- function(x, ind = NULL, cutoffs = NULL, m = 5,
                             margin_mode = c("sd_shift", "normal"),
                             max_iter = 100, tol = 1e-4, seed = 1, shift_k = 3,
                             margin_candidates = c("norm", "logis", "lnorm", "gamma",
                               "weibull", "exp", "invgauss", "gengamma", "llogis", "lomax", "burr"),
                             ag = "gibbs", thinning = 2, verbose = FALSE,
                             lyles_control = list(), input_scale = c("log", "raw"),
                             optim_methods = c("L-BFGS-B", "Nelder-Mead"),
                              allow_unconverged = FALSE, ...) {
  .check_dots(...)
  .check_flag(allow_unconverged, "allow_unconverged")
  .check_number(m, "m", 1, .Machine$integer.max, TRUE)
  .check_number(max_iter, "max_iter", 0, .Machine$integer.max, TRUE)
  .check_number(tol, "tol", .Machine$double.eps)
  .check_number(thinning, "thinning", 1, .Machine$integer.max, TRUE)
  .check_flag(verbose, "verbose")
  if (!is.null(seed)) .check_number(seed, "seed", 0, .Machine$integer.max, TRUE)
  ag <- match.arg(ag, c("gibbs", "gibbsR", "rejection"))
  lyles_control <- .check_lyles_control(lyles_control)
  if (inherits(x, "copmi_lod_data")) {
    if (!is.null(ind) || !is.null(cutoffs)) {
      stop("Do not supply ind/cutoffs again when x is a copmi_lod_data object.", call. = FALSE)
    }
    dat <- .as_lod_data(x)
  } else {
    if (is.null(ind) || is.null(cutoffs)) stop("ind and cutoffs are required for matrix x.", call. = FALSE)
    dat <- make_lod_data(x, ind, cutoffs)
  }
  margins <- copmi_fit_margins(dat, .match_margin_mode(margin_mode), match.arg(input_scale),
                              shift_k, margin_candidates, optim_methods)
  model <- copmi_fit_copula(copmi_transform(margins), max_iter, tol, seed, verbose, lyles_control)
  result <- copmi_impute(model, m, ag = ag, thinning = thinning, allow_unconverged = allow_unconverged)
  result$call <- match.call()
  result
}
