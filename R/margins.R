#' Fit and compare left-censored marginal distributions
#'
#' Fit a single exposure on the supplied scale. This function does not log or
#' shift the data. Pass `candidates = "norm"` for a single censored normal fit,
#' or multiple families for minimum-BIC selection.
#' @param x Numeric vector of exposure values. Censored values may be `NA`.
#' @param ind Numeric or logical vector, with `1` observed and `0` censored.
#' @param cutoff One finite detection limit; `NA` if all values are observed.
#' @param candidates Unique candidate names: `norm`, `logis`, `lnorm`, `gamma`,
#'   `weibull`, `exp`, `invgauss`, `gengamma`, `llogis`, `lomax`, or `burr`.
#' @param optim_methods Ordered optimization methods. The default first tries
#'   `"L-BFGS-B"`, then `"Nelder-Mead"` on error, nonconvergence, or invalid
#'   estimates/log-likelihood/BIC. Both optimize the same left-censored
#'   likelihood on log-transformed positive parameters, without empirical
#'   parameter bounds. Location parameters remain unconstrained. A single
#'   method can be supplied for controlled comparisons.
#' @return A `copmi_margin` list containing `family`, `parameters` (named numeric
#'   vector), `bic`, `loglik`, `optimizer`, `how`, `candidate_table` (one row per
#'   candidate, including failure messages), and `optimization` (one row per
#'   starting-value/optimizer attempt). Failed candidates have `bic = Inf`.
#' @details At least two observed values are required. Nonzero observed variation is
#'   required except for the exponential family. Positive-support families require positive observed values and
#'   a positive cutoff when censoring is present. Failed candidates are retained
#'   in the diagnostic table but cannot win BIC selection. If every candidate
#'   fails, the function raises an error. Family-specific multiple starting values are retained. Optimization warnings are captured in
#'   `optimization$warnings`, rather than printed during fitting.
#' @seealso [copmi_fit_margins()], [copmi_transform()]
#' @export
#' @examples
#' d <- nhanes_pah
#' marginal <- copmi_fit_margin(d$X_cens[, 1], d$ind[, 1], d$cutoffs[1],
#'                              candidates = "norm")
#' marginal$parameters
#' marginal$candidate_table
#' marginal$optimization
copmi_fit_margin <- function(x, ind, cutoff,
                             candidates = c("norm", "logis", "lnorm", "gamma",
                                            "weibull", "exp", "invgauss", "gengamma",
                                            "llogis", "lomax", "burr"),
                             optim_methods = c("L-BFGS-B", "Nelder-Mead")) {
  candidates <- .check_candidates(candidates)
  optim_methods <- .check_optim_methods(optim_methods)
  if (!is.numeric(x) || !is.null(dim(x)) || !length(x) ||
      !(is.numeric(ind) || is.logical(ind)) || !is.null(dim(ind)) ||
      length(ind) != length(x)) {
    stop("x and ind must be vectors of the same positive length.", call. = FALSE)
  }
  dat <- make_lod_data(matrix(x, ncol = 1), matrix(ind, ncol = 1), unname(cutoff))
  x <- dat$X_cens[, 1]
  ind <- dat$ind[, 1]
  if (sum(ind == 1L) < 2L) stop("At least two observed values are required per margin.", call. = FALSE)
  results <- vector("list", length(candidates))
  table_rows <- optimization_rows <- results
  for (i in seq_along(candidates)) {
    family <- candidates[i]
    result <- tryCatch(.fit_margin_family(x, ind == 1L, cutoff, family, optim_methods),
      error = function(e) list(fit = NULL, optimization = .empty_optimization(),
                                message = conditionMessage(e)))
    fit <- result$fit
    results[i] <- list(fit)
    table_rows[[i]] <- data.frame(candidate = family,
      bic = if (is.null(fit)) Inf else fit$bic,
      converged = !is.null(fit),
      optimizer = if (is.null(fit)) NA_character_ else fit$optimizer,
      message = result$message)
    attempts <- result$optimization
    attempts$candidate <- rep(family, nrow(attempts))
    optimization_rows[[i]] <- attempts[, c("candidate", setdiff(names(attempts), "candidate"))]
  }
  candidate_table <- do.call(rbind, table_rows)
  if (!any(is.finite(candidate_table$bic))) {
    messages <- unlist(lapply(optimization_rows, function(z) z$message[nzchar(z$message)]))
    stop("All candidate margin fits failed: ", paste(candidate_table$message, collapse = " | "),
         if (length(messages)) paste0("; ", summarize_messages(messages)), call. = FALSE)
  }
  best <- results[[which.min(candidate_table$bic)]]
  structure(list(family = best$family, parameters = best$parameters, bic = best$bic,
    loglik = best$loglik, optimizer = best$optimizer, how = best$how,
    candidate_table = candidate_table, optimization = do.call(rbind, optimization_rows)),
    class = "copmi_margin")
}

#' Fit all margins under a log-normal assumption or by BIC selection
#'
#' Separately choose the supplied data scale and the marginal assumption.
#' `margin_mode = "normal"` fits censored normal margins on the log scale.
#' `margin_mode = "sd_shift"` shifts log-scale values and selects a margin for
#' each column by BIC. The alias `"select"` is also accepted for the latter.
#' @param x A [make_lod_data()] object.
#' @param margin_mode `"normal"` if the user assumes normality after logging,
#'   or `"sd_shift"`/`"select"` for heterogeneous candidate selection. This is a
#'   user-specified modelling assumption, not an automatic normality test.
#' @param input_scale `"log"` (default) if values and cutoffs are already logged;
#'   `"raw"` to apply natural log internally to strictly positive values and
#'   finite cutoffs. Imputations are returned on this supplied scale.
#' @param shift_k Positive SD-shift multiplier, used only for selection.
#' @param margin_candidates Candidate family names, used only for selection.
#'   See [copmi_fit_margin()]. In normal mode the sole candidate is `norm`.
#' @inheritParams copmi_fit_margin
#' @return A `copmi_margins` list with `fits` (named list of `copmi_margin`
#'   objects), `family_selected` (named character vector), `candidate_table` and
#'   `optimization` (tables with a `variable` column), `input` (original data),
#'   `analysis_data` (log-scale data), `work_data` (data used for fitting),
#'   `shift` (`copmi_shift` or `NULL`), `margin_mode`, `input_scale`, and
#'   `optim_methods`. No copula estimation or random imputation occurs here.
#' @seealso [copmi_shift()], [copmi_transform()], [copula_em_impute()]
#' @export
#' @examples
#' margins <- copmi_fit_margins(nhanes_pah, margin_mode = "normal")
#' margins$family_selected
#' margins$fits[[1]]$parameters
#' selected <- copmi_fit_margins(nhanes_pah, margin_mode = "select",
#'                               margin_candidates = c("norm", "logis"))
#' selected$candidate_table
copmi_fit_margins <- function(x, margin_mode = c("sd_shift", "normal"),
                              input_scale = c("log", "raw"), shift_k = 3,
                              margin_candidates = c("norm", "logis", "lnorm", "gamma",
                                "weibull", "exp", "invgauss", "gengamma", "llogis", "lomax", "burr"),
                              optim_methods = c("L-BFGS-B", "Nelder-Mead")) {
  x <- .as_lod_data(x)
  margin_mode <- .match_margin_mode(margin_mode)
  input_scale <- match.arg(input_scale)
  optim_methods <- .check_optim_methods(optim_methods)
  analysis_data <- x
  if (input_scale == "raw") {
    if (any(x$X_cens <= 0) || any(x$cutoffs[is.finite(x$cutoffs)] <= 0)) {
      stop("input_scale = 'raw' requires positive values and finite cutoffs before log().", call. = FALSE)
    }
    analysis_data$X_cens <- log(x$X_cens)
    analysis_data$cutoffs <- log(x$cutoffs)
    analysis_data$scale <- "log"
  }
  shift <- NULL
  work_data <- analysis_data
  if (margin_mode == "sd_shift") {
    margin_candidates <- .check_candidates(margin_candidates)
    shift <- copmi_shift(analysis_data, shift_k)
    work_data <- shift$dat
  } else margin_candidates <- "norm"
  nm <- colnames(x$X_cens)
  fits <- lapply(seq_along(nm), function(j) {
    tryCatch(copmi_fit_margin(work_data$X_cens[, j], work_data$ind[, j],
      work_data$cutoffs[j], candidates = margin_candidates, optim_methods = optim_methods),
      error = function(e) stop("Margin ", nm[j], ": ", conditionMessage(e), call. = FALSE))
  })
  names(fits) <- nm
  bind_tables <- function(field) {
    out <- do.call(rbind, lapply(seq_along(fits), function(j) {
      tab <- fits[[j]][[field]]
      data.frame(variable = rep(nm[j], nrow(tab)), tab, row.names = NULL,
                 stringsAsFactors = FALSE)
    }))
    rownames(out) <- NULL
    out
  }
  structure(list(fits = fits, family_selected = vapply(fits, `[[`, character(1), "family"),
    candidate_table = bind_tables("candidate_table"), optimization = bind_tables("optimization"),
    input = x, analysis_data = analysis_data, work_data = work_data, shift = shift,
    margin_mode = margin_mode, input_scale = input_scale, optim_methods = optim_methods),
    class = "copmi_margins")
}

.match_margin_mode <- function(mode) {
  if (identical(mode, "select")) return("sd_shift")
  match.arg(mode, c("sd_shift", "normal"))
}
