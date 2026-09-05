# Family-specific initialization only. The likelihood and optimizer are shared.
# A positive family rejects incompatible observations instead of clipping them.
.margin_starts <- function(obs, cutoff, has_censored, family) {
  if (length(obs) < 2L) stop("At least two observed values are required.", call. = FALSE)
  if (!family %in% c("norm", "logis") &&
      (any(obs <= 0) || (has_censored && cutoff <= 0))) {
    stop("Positive-support margins require positive observations and a positive cutoff.", call. = FALSE)
  }
  m <- mean(obs)
  s <- stats::sd(obs)
  if (!is.finite(s) || (s == 0 && family != "exp")) {
    stop("A margin must have finite, nonzero observed variation.", call. = FALSE)
  }
  med <- stats::median(obs)
  switch(family,
    norm = list(c(mean = m, sd = s)),
    logis = list(c(location = med, scale = s * sqrt(3) / pi)),
    lnorm = list(c(meanlog = mean(log(obs)), sdlog = stats::sd(log(obs)))),
    exp = list(c(rate = 1 / m)),
    weibull = lapply(c(0.8, 1.1, 1.6), function(a)
      c(shape = a, scale = m / gamma(1 + 1 / a))),
    gamma = lapply(c(0.6, 1, 1.8), function(a) {
      shape <- a * (m / s)^2
      c(shape = shape, scale = m / shape)
    }),
    invgauss = lapply(c(0.5, 1, 2), function(a)
      c(mean = m, shape = a * m * (m / s)^2)),
    llogis = lapply(c(1.2, 2, 3), function(a) c(shape = a, scale = med)),
    burr = list(c(shape1 = 1.2, shape2 = 2, scale = med),
                c(shape1 = 1.8, shape2 = 1.5, scale = med),
                c(shape1 = 2.5, shape2 = 1.2, scale = m)),
    gengamma = list(c(shape = (m / s)^2, scale = s^2 / m, k = 1),
                   c(shape = 1, scale = m, k = 1.2),
                   c(shape = 1.5, scale = med, k = 0.8)),
    lomax = lapply(c(2.2, 2.8, 3.5), function(a) c(shape = a, scale = m * (a - 1)))
  )
}
