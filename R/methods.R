#' Print and summarize CopMI objects
#' @param x A CopMI data, marginal, latent, fitted copula, or imputation object.
#' @param object A `copmi_mi` result to summarize.
#' @param ... Additional arguments, currently unused.
#' @return `print()` displays key dimensions and fit information, then invisibly
#'   returns its input. `summary()` returns a one-row data frame with `method`,
#'   `margin_mode`, `input_scale`, `m`, `converged`, `n_iter`, `n`, `p`, and
#'   `n_fallbacks` (total number of recorded retry, candidate-failure, and initial-projection events).
#' @name copmi-methods
#' @examples
#' # Use the first 100 rows for a quick help example; full data have 1330 rows.
#' dat <- make_lod_data(nhanes_pah$X_cens[1:100, ],
#'                      nhanes_pah$ind[1:100, ], nhanes_pah$cutoffs, scale = "log")
#' print(dat)
#' fit <- copula_em_impute(dat, margin_mode = "normal",
#'                          max_iter = 100, m = 1)
#' print(fit)
#' summary(fit)
NULL

#' @rdname copmi-methods
#' @export
print.copmi_lod_data <- function(x, ...) {
  cat("CopMI left-censored data\n  n:", nrow(x$X_cens), " p:", ncol(x$X_cens),
      "\n  censored cells:", sum(x$ind == 0L), "\n  scale label:", x$scale, "\n")
  invisible(x)
}

#' @rdname copmi-methods
#' @export
print.copmi_margin <- function(x, ...) {
  cat("CopMI fitted margin\n  family:", x$family, " optimizer:", x$optimizer,
      "\n  BIC:", format(x$bic), "\n")
  print(x$parameters)
  invisible(x)
}

#' @rdname copmi-methods
#' @export
print.copmi_margins <- function(x, ...) {
  cat("CopMI fitted margins\n  mode:", x$margin_mode, " input scale:", x$input_scale, "\n")
  print(x$family_selected)
  invisible(x)
}

#' @rdname copmi-methods
#' @export
print.copmi_shift <- function(x, ...) {
  cat("CopMI SD shift\n  multiplier:", x$shift_multiplier, "\n  anchors:\n")
  print(x$anchor)
  invisible(x)
}

#' @rdname copmi-methods
#' @export
print.copmi_latent <- function(x, ...) {
  cat("CopMI latent Gaussian scores\n  n:", nrow(x$Z), " p:", ncol(x$Z),
      "\n  latent cutoffs:\n")
  print(x$z_lod)
  invisible(x)
}

#' @rdname copmi-methods
#' @export
print.copmi_copula <- function(x, ...) {
  cat("CopMI fitted copula\n  converged:", x$converged, " after ", x$n_iter,
      " EM iteration(s)\n", sep = "")
  print(x$Sigma_hat)
  invisible(x)
}

#' @rdname copmi-methods
#' @export
print.copmi_mi <- function(x, ...) {
  cat("CopMI multiple-imputation fit\n  method:", x$method,
      "\n  margin mode:", x$margin_mode, " input scale:", x$input_scale,
      "\n  imputations:", x$M, "\n  converged:", x$converged,
      " after ", x$n_iter, " EM iteration(s)\n", sep = " ")
  invisible(x)
}

#' @rdname copmi-methods
#' @export
summary.copmi_mi <- function(object, ...) {
  data.frame(method = object$method, margin_mode = object$margin_mode,
    input_scale = object$input_scale, m = object$M, converged = object$converged,
    n_iter = object$n_iter, n = nrow(object$input$X_cens), p = ncol(object$input$X_cens),
    n_fallbacks = sum(object$fallback_records$count), stringsAsFactors = FALSE)
}
