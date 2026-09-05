.check_optim_methods <- function(x) {
  allowed <- c("L-BFGS-B", "Nelder-Mead")
  if (!is.character(x) || !length(x) || anyNA(x) ||
      anyDuplicated(x) || any(!x %in% allowed)) {
    stop("optim_methods must be one or both of L-BFGS-B and Nelder-Mead, in order.", call. = FALSE)
  }
  x
}

.distribution <- function(family) {
  switch(family,
    norm = list(d = stats::dnorm, p = stats::pnorm, q = stats::qnorm),
    logis = list(d = stats::dlogis, p = stats::plogis, q = stats::qlogis),
    lnorm = list(d = stats::dlnorm, p = stats::plnorm, q = stats::qlnorm),
    gamma = list(d = stats::dgamma, p = stats::pgamma, q = stats::qgamma),
    weibull = list(d = stats::dweibull, p = stats::pweibull, q = stats::qweibull),
    exp = list(d = stats::dexp, p = stats::pexp, q = stats::qexp),
    invgauss = list(d = statmod::dinvgauss, p = statmod::pinvgauss, q = statmod::qinvgauss),
    gengamma = list(d = dgengamma, p = pgengamma, q = qgengamma),
    llogis = list(d = dllogis, p = pllogis, q = qllogis),
    lomax = list(d = dlomax, p = plomax, q = qlomax),
    burr = list(d = dburr, p = pburr, q = qburr),
    stop("Unsupported distribution: ", family, call. = FALSE)
  )
}

.empty_optimization <- function() {
  data.frame(start = integer(), optimizer = character(), converged = logical(),
             bic = double(), message = character(), warnings = character())
}

# This is the only optimizer error boundary. Invalid trials receive a finite
# penalty for L-BFGS-B; a final penalized value is never accepted as a fit.
.run_optimizer <- function(start, objective, method, maxit = 2000L) {
  messages <- character()
  fit <- tryCatch(withCallingHandlers(
    stats::optim(start, objective, method = method, control = list(maxit = maxit)),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) list(convergence = 99L, message = conditionMessage(e)))
  fit$valid <- isTRUE(fit$convergence == 0L) &&
    length(fit$par) == length(start) && all(is.finite(fit$par)) &&
    isTRUE(is.finite(fit$value) && fit$value < 1e100)
  fit$message <- fit$message %||% if (fit$valid) "" else
    paste0("Optimization failed (code ", fit$convergence, ").")
  fit$warnings <- paste(unique(messages), collapse = " | ")
  fit
}

.censored_loglik <- function(parameters, family, x, observed, cutoff) {
  dist <- .distribution(family)
  pars <- as.list(parameters)
  sum(do.call(dist$d, c(list(x = x[observed], log = TRUE), pars))) +
    if (any(!observed)) sum(!observed) *
      do.call(dist$p, c(list(q = cutoff, log.p = TRUE), pars)) else 0
}

.fit_margin_family <- function(x, observed, cutoff, family, optim_methods) {
  starts <- .margin_starts(x[observed], cutoff, any(!observed), family)
  positive <- !names(starts[[1]]) %in% c("meanlog", "location")
  if (family == "norm") positive <- names(starts[[1]]) != "mean"
  location_center <- starts[[1]][!positive]
  location_scale <- if (family == "norm") starts[[1]]["sd"] else
    if (family == "logis") starts[[1]]["scale"] else 1
  decode <- function(theta) {
    theta[!positive] <- location_center + location_scale * theta[!positive]
    theta[positive] <- exp(theta[positive])
    theta
  }
  objective <- function(theta) {
    parameters <- decode(theta)
    if (any(!is.finite(parameters)) || any(parameters[positive] <= 0)) return(1e100)
    ll <- .censored_loglik(parameters, family, x, observed, cutoff)
    if (is.finite(ll)) -ll / length(x) else 1e100
  }
  attempts <- .empty_optimization()
  best <- NULL
  for (i in seq_along(starts)) {
    start <- starts[[i]]
    start[!positive] <- (start[!positive] - location_center) / location_scale
    start[positive] <- log(start[positive])
    for (method in optim_methods) {
      fit <- .run_optimizer(start, objective, method)
      bic <- if (fit$valid) 2 * fit$value * length(x) + length(start) * log(length(x)) else Inf
      attempts <- rbind(attempts, data.frame(start = i, optimizer = method,
        converged = fit$valid, bic = bic, message = fit$message, warnings = fit$warnings))
      if (fit$valid) {
        if (is.null(best) || bic < best$bic) {
          best <- list(family = family, parameters = decode(fit$par), bic = bic,
            loglik = -fit$value * length(x), optimizer = method,
            how = "censored_maximum_likelihood")
        }
        break
      }
    }
  }
  list(fit = best, optimization = attempts,
       message = if (is.null(best)) summarize_messages(attempts$message) else "")
}
