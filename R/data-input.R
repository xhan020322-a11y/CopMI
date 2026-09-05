#' Validate and construct left-censored exposure data
#'
#' Construct the common input to all CopMI workflows. No log or exponential
#' transformation is performed: values and cutoffs must already use the same
#' analysis scale. The indicator, not the placeholder value, defines censoring.
#'
#' @param X_cens Numeric matrix or data frame, with samples in rows and exposure
#'   variables in columns. Observed cells must be finite; censored cells may be
#'   finite or `NA` and are replaced by the corresponding cutoff.
#' @param ind Matrix or data frame of the same dimensions, containing exactly
#'   `0`/`FALSE` (left-censored) and `1`/`TRUE` (observed). Fractional values and
#'   missing indicators are rejected before any integer conversion.
#' @param cutoffs Numeric vector, one cutoff per column, in column order. A
#'   named vector must match the column names in order. Use a finite cutoff for
#'   each censored column; `NA` is allowed for a fully observed column.
#' @param scale One nonempty character label, such as `"analysis"` or `"log"`.
#'
#' @return A `copmi_lod_data` list with `X_cens` (double matrix), `ind` (integer
#'   matrix), `cutoffs` (numeric vector), and `scale`.
#'   Matrix dimensions and sample/variable names are preserved. If absent,
#'   column names are generated as `X1`, `X2`, etc.
#' @details Observed values below a finite cutoff are inconsistent with the
#'   left-censoring model and are rejected. Observed values equal to the cutoff
#'   are allowed. Column names must be unique and cannot be `.imp` or `.id`,
#'   which are reserved by [complete()]. General missingness, right censoring,
#'   and cell-specific cutoffs are not represented by this input contract.
#' @seealso [copula_em_impute()], [copmi_fit_margins()]
#' @export
#' @examples
#' dat <- make_lod_data(nhanes_pah$X_cens,
#'                      nhanes_pah$ind, nhanes_pah$cutoffs)
#' print(dat)
#' dim(dat$X_cens)
#' colSums(dat$ind == 0L)
make_lod_data <- function(X_cens, ind, cutoffs, scale = "analysis") {
  X_cens <- .numeric_matrix(X_cens, "X_cens")
  if (!(is.matrix(ind) || is.data.frame(ind))) {
    stop("ind must be a matrix or data frame.", call. = FALSE)
  }
  ind <- as.matrix(ind)
  if (!(is.numeric(ind) || is.logical(ind)) || !identical(dim(ind), dim(X_cens)) ||
      anyNA(ind) || any(!ind %in% c(0, 1))) {
    stop("ind must have the same dimensions as X_cens and contain only 0/1 or FALSE/TRUE.",
         call. = FALSE)
  }
  for (k in 1:2) {
    if (!is.null(dimnames(ind)[[k]]) && !is.null(dimnames(X_cens)[[k]]) &&
        !identical(dimnames(ind)[[k]], dimnames(X_cens)[[k]])) {
      stop("ind dimnames must match X_cens in the same order.", call. = FALSE)
    }
  }
  if (is.null(colnames(X_cens))) colnames(X_cens) <- paste0("X", seq_len(ncol(X_cens)))
  nm <- colnames(X_cens)
  if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) || any(nm %in% c(".imp", ".id"))) {
    stop("Column names must be nonempty, unique, and different from .imp and .id.", call. = FALSE)
  }
  if (!is.numeric(cutoffs) || length(cutoffs) != ncol(X_cens) || any(is.infinite(cutoffs))) {
    stop("cutoffs must contain one finite number or NA per column.", call. = FALSE)
  }
  if (!is.null(names(cutoffs)) && !identical(names(cutoffs), nm)) {
    stop("Named cutoffs must match X_cens columns in the same order.", call. = FALSE)
  }
  if (!is.character(scale) || length(scale) != 1L || is.na(scale) || !nzchar(scale)) {
    stop("scale must be one nonempty character label.", call. = FALSE)
  }
  if (any(!is.finite(X_cens[ind == 1L])) || any(is.infinite(X_cens))) {
    stop("Observed X_cens values must be finite; only censored cells may contain NA.", call. = FALSE)
  }
  for (j in seq_len(ncol(X_cens))) {
    censored <- ind[, j] == 0L
    if (any(censored) && !is.finite(cutoffs[j])) {
      stop("Censored column ", nm[j], " requires a finite cutoff.", call. = FALSE)
    }
    if (is.finite(cutoffs[j]) && any(X_cens[!censored, j] < cutoffs[j])) {
      stop("Observed values in ", nm[j], " are below the cutoff.", call. = FALSE)
    }
    X_cens[censored, j] <- cutoffs[j]
  }
  storage.mode(ind) <- "integer"
  dimnames(ind) <- dimnames(X_cens)
  structure(list(X_cens = X_cens, ind = ind, cutoffs = as.numeric(cutoffs),
                 scale = scale), class = "copmi_lod_data")
}

#' Construct the SD shift used by heterogeneous margins
#'
#' Shift each column by `anchor = reference - shift_k * observed_sd` and return
#' both shifted values and the information needed to undo the shift.
#' @param x A [make_lod_data()] object.
#' @param shift_k Positive finite multiplier of the observed standard deviation.
#' @return A `copmi_shift` list with `dat` (shifted `copmi_lod_data`), `anchor`,
#'   `sd`, `sd_source`, `shift_multiplier`,
#'   `min_shifted_observed_or_cutoff`, `original_cutoffs`, and `shifted_cutoffs`.
#'   Each vector has one element per column. Recover the input scale by adding
#'   `anchor` to the shifted columns with `sweep(..., 2, anchor, "+")`.
#' @details The reference is the finite cutoff, or the minimum observed value
#'   if the cutoff is `NA`. SD is estimated from each column's observed cells only.
#'   An unavailable or zero SD raises an error; scales are never borrowed from
#'   other variables. All visible shifted values must be positive. This is an
#'   additive shift, not a log transformation. Positive-support selected families
#'   imply values above `anchor`; selected normal/logistic families retain their
#'   full real support and do not impose this extra lower bound.
#' @seealso [copmi_fit_margins()]
#' @export
#' @examples
#' shifted <- copmi_shift(nhanes_pah)
#' shifted$anchor
#' restored <- sweep(shifted$dat$X_cens, 2, shifted$anchor, "+")
#' all.equal(restored, nhanes_pah$X_cens)
copmi_shift <- function(x, shift_k = 3) {
  x <- .as_lod_data(x)
  .check_number(shift_k, "shift_k", .Machine$double.eps)
  structure(make_log_sd_shift_dat(x, shift_k), class = "copmi_shift")
}
