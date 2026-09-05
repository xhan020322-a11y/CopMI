.check_number <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if ((!is.numeric(x) || is.complex(x)) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    stop(name, " must be a finite ", if (integer) "integer" else "number",
         " in [", lower, ", ", upper, "].", call. = FALSE)
  }
  invisible(x)
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

.check_dots <- function(...) {
  dots <- list(...)
  if (length(dots)) stop("Unused argument(s) in ...: ",
                         paste(names(dots) %||% "unnamed", collapse = ", "),
                         call. = FALSE)
}

.numeric_matrix <- function(x, name) {
  if (!(is.matrix(x) || is.data.frame(x)) ||
      (is.data.frame(x) && !all(vapply(x, is.numeric, logical(1)))) ||
      (is.matrix(x) && !is.numeric(x))) {
    stop(name, " must be a numeric matrix or numeric data frame.", call. = FALSE)
  }
  x <- as.matrix(x)
  if (is.complex(x)) stop(name, " must contain real numbers.", call. = FALSE)
  if (length(dim(x)) != 2L || any(dim(x) == 0L)) {
    stop(name, " must have at least one row and one column.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.check_class <- function(x, class, name) {
  if (!inherits(x, class)) stop(name, " must be a ", class, " object.", call. = FALSE)
  invisible(x)
}

.as_lod_data <- function(x) {
  .check_class(x, "copmi_lod_data", "x")
  make_lod_data(x$X_cens, x$ind, x$cutoffs, scale = x$scale)
}

.margin_families <- function() {
  c("norm", "logis", "lnorm", "gamma", "weibull", "exp", "invgauss",
    "gengamma", "llogis", "lomax", "burr")
}

.check_candidates <- function(x) {
  if (!is.character(x) || !length(x) || anyNA(x) || anyDuplicated(x) ||
      any(!x %in% .margin_families())) {
    stop("candidates must be unique names from: ",
         paste(.margin_families(), collapse = ", "), ".", call. = FALSE)
  }
  x
}

.check_lyles_control <- function(x) {
  defaults <- list(min_obs = 2L, maxit = 2000L, rho_max = 1)
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      anyDuplicated(names(x)) || any(!names(x) %in% names(defaults))))) {
    stop("lyles_control must be a named list with min_obs, maxit, or rho_max.",
         call. = FALSE)
  }
  for (name in names(x)) if (is.null(x[[name]])) {
    stop("lyles_control entries cannot be NULL.", call. = FALSE)
  }
  x <- utils::modifyList(defaults, x)
  .check_number(x$min_obs, "lyles_control$min_obs", 2, .Machine$integer.max, TRUE)
  .check_number(x$maxit, "lyles_control$maxit", 1, .Machine$integer.max, TRUE)
  .check_number(x$rho_max, "lyles_control$rho_max", .Machine$double.eps, 1)
  x
}

# Restore the caller's RNG state for explicitly seeded operations and for
# sampling continued from a saved fit. NULL without a state uses the caller RNG.
.with_rng <- function(seed = NULL, state = NULL, code) {
  if (!is.null(seed)) .check_number(seed, "seed", 0, .Machine$integer.max, TRUE)
  local_rng <- !is.null(seed) || !is.null(state)
  if (local_rng) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    if (!is.null(seed)) set.seed(seed)
    else assign(".Random.seed", state, envir = .GlobalEnv)
  }
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) stats::runif(1)
  force(code)
}
