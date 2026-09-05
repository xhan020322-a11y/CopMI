#' Extract completed data from a CopMI result
#' @param object A `copmi_mi` object from [copula_em_impute()] or [copmi_impute()].
#' @param action One positive integer selecting an imputation, or `"all"`,
#'   `"mean"`, or `"long"`.
#' @param ... Reserved; unused arguments raise an error.
#' @return For an integer, an `n` by `p` matrix; for `"all"`, a list of `m`
#'   matrices; for `"mean"`, an `n` by `p` elementwise average; for `"long"`,
#'   a data frame with `n * m` rows and columns `.imp` (1 through `m`), `.id`
#'   (original row number), then the exposure variables in their original order.
#' @details The mean completion is a descriptive convenience. An analysis of
#'   this average does not propagate between-imputation uncertainty. CopMI does
#'   not pool downstream model estimates or standard errors. Use
#'   `CopMI::complete()` if another package defines a function named `complete`.
#' @export
#' @examples
#' # Use the first 100 rows for a quick help example; full data have 1330 rows.
#' dat <- make_lod_data(nhanes_pah$X_cens[1:100, ],
#'                      nhanes_pah$ind[1:100, ], nhanes_pah$cutoffs, scale = "log")
#' fit <- copula_em_impute(dat, margin_mode = "normal",
#'                          m = 2, max_iter = 100)
#' dim(complete(fit, 1))
#' length(complete(fit, "all"))
#' dim(complete(fit, "mean"))
#' head(complete(fit, "long"))
complete <- function(object, action = "all", ...) UseMethod("complete")

#' @rdname complete
#' @export
complete.copmi_mi <- function(object, action = "all", ...) {
  .check_dots(...)
  if (is.numeric(action)) {
    .check_number(action, "action", 1, length(object$imp_list), TRUE)
    return(object$imp_list[[as.integer(action)]])
  }
  action <- match.arg(action, c("all", "mean", "long"))
  if (action == "all") return(object$imp_list)
  if (action == "mean") {
    out <- Reduce(`+`, object$imp_list) / length(object$imp_list)
    return(out)
  }
  rows <- lapply(seq_along(object$imp_list), function(i) {
    df <- as.data.frame(object$imp_list[[i]], stringsAsFactors = FALSE, optional = TRUE)
    data.frame(.imp = i, .id = seq_len(nrow(df)), df,
               row.names = NULL, check.names = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Extract fitting, optimization, and imputation diagnostics
#' @param object A `copmi_mi`, `copmi_copula`, `copmi_margins`, or `copmi_margin`
#'   object returned by the corresponding fitting function.
#' @param ... Reserved; unused arguments raise an error.
#' @return A named list. For imputation results: `converged`, `n_iter`,
#'   `Sigma_hat`, `Sigma_init`, `em_change_history`, `family_selected`,
#'   `bic_table`, `candidate_table`, `optimization`, `fallback_records`,
#'   `init_diagnostics`, `shift_anchor`, `shift_sd`, `shift_sd_source`,
#'   `input_scale`, and `sampling`. The fallback table has columns `method`,
#'   `stage`, `count`, and `detail`; `count` counts events, not necessarily
#'   distinct rows or variables. Marginal objects return the relevant candidate
#'   and optimization tables; copula objects return correlation and EM details.
#' @details A successful Nelder-Mead retry remains visible in `optimization`
#'   even if the final fit converges. The fallback table records marginal and initial-pair optimizer
#'   retries, failed marginal candidates, and initial-matrix projection. EM or
#'   sampling failures raise errors; they do not produce replacement values. Consult the per-attempt `message` and `warnings`
#'   columns when diagnosing a failed candidate.
#' @export
#' @examples
#' # Use the first 100 rows for a quick help example; full data have 1330 rows.
#' dat <- make_lod_data(nhanes_pah$X_cens[1:100, ],
#'                      nhanes_pah$ind[1:100, ], nhanes_pah$cutoffs, scale = "log")
#' fit <- copula_em_impute(dat, margin_mode = "normal",
#'                          m = 1, max_iter = 100)
#' dx <- diagnostics(fit)
#' dx$converged
#' dx$optimization
#' dx$fallback_records
diagnostics <- function(object, ...) UseMethod("diagnostics")

#' @rdname diagnostics
#' @export
diagnostics.copmi_mi <- function(object, ...) {
  .check_dots(...)
  list(converged = object$converged, n_iter = object$n_iter,
    Sigma_hat = object$Sigma_hat, Sigma_init = object$extra$Sigma_init,
    em_change_history = object$extra$em_change_history,
    family_selected = object$family_selected, bic_table = object$bic_table,
    candidate_table = object$extra$candidate_table, optimization = object$extra$optimization,
    fallback_records = object$fallback_records, init_diagnostics = object$extra$init_diagnostics,
    shift_anchor = object$shift_anchor, shift_sd = object$extra$shift_sd,
    shift_sd_source = object$extra$shift_sd_source, input_scale = object$input_scale,
    sampling = object$extra$sampling)
}

#' @rdname diagnostics
#' @export
diagnostics.copmi_copula <- function(object, ...) {
  .check_dots(...)
  object[c("converged", "n_iter", "Sigma_hat", "Sigma_init", "em_change_history",
           "init_diagnostics", "fallback_records")]
}

#' @rdname diagnostics
#' @export
diagnostics.copmi_margins <- function(object, ...) {
  .check_dots(...)
  object[c("family_selected", "candidate_table", "optimization", "margin_mode", "input_scale")]
}

#' @rdname diagnostics
#' @export
diagnostics.copmi_margin <- function(object, ...) {
  .check_dots(...)
  object[c("family", "parameters", "bic", "loglik", "optimizer", "candidate_table", "optimization")]
}
