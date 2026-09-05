# Stable distribution functions. Preserve true support, endpoints and log tails.
.softplus <- function(x) -stats::plogis(-x, log.p = TRUE)
.log_expm1 <- function(x) x + log(-expm1(-x))

# Density evaluation on nonnegative support, including the limit at zero.
.positive_density <- function(x, parameters, log_density, at_zero, log) {
  if (any(!is.finite(parameters)) || any(parameters <= 0)) return(rep(NaN, length(x)))
  out <- rep(-Inf, length(x))
  out[is.na(x)] <- NA_real_
  inside <- is.finite(x) & x > 0
  out[inside] <- log_density(x[inside])
  out[which(x == 0)] <- at_zero
  if (log) out else exp(out)
}

# log(softplus(x)) and log(1 - exp(-exp(x))) use their analytic limits
# when the correction is smaller than floating-point precision.
.log_softplus <- function(x) {
  out <- x
  regular <- which(x > log(.Machine$double.eps))
  out[regular] <- log(.softplus(x[regular]))
  out
}
.from_log_hazard <- function(log_hazard, lower.tail, log.p) {
  if (!lower.tail) return(if (log.p) -exp(log_hazard) else exp(-exp(log_hazard)))
  log_cdf <- log_hazard
  regular <- which(log_hazard > log(.Machine$double.eps))
  log_cdf[regular] <- log(-expm1(-exp(log_hazard[regular])))
  if (log.p) log_cdf else exp(log_cdf)
}

dlomax <- function(x, shape, scale, log = FALSE) {
  .positive_density(x, c(shape, scale), function(x)
    log(shape) - log(scale) - (shape + 1) * .softplus(log(x) - log(scale)),
    log(shape) - log(scale), log)
}
plomax <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE) {
  if (any(!is.finite(c(shape, scale))) || shape <= 0 || scale <= 0) return(rep(NaN, length(q)))
  lh <- log(shape) + .log_softplus(log(pmax(q, 0)) - log(scale))
  .from_log_hazard(lh, lower.tail, log.p)
}
qlomax <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE) {
  if (any(!is.finite(c(shape, scale))) || shape <= 0 || scale <= 0) return(rep(NaN, length(p)))
  exp(log(scale) + .log_expm1(stats::qexp(p, rate = shape, lower.tail = lower.tail, log.p = log.p)))
}

dllogis <- function(x, shape, scale, log = FALSE) {
  .positive_density(x, c(shape, scale), function(x)
    stats::dlogis(log(x), location = log(scale), scale = 1 / shape, log = TRUE) - log(x),
    if (shape == 1) -log(scale) else if (shape < 1) Inf else -Inf, log)
}
pllogis <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE) {
  if (any(!is.finite(c(shape, scale))) || shape <= 0 || scale <= 0) return(rep(NaN, length(q)))
  stats::plogis(log(pmax(q, 0)), location = log(scale), scale = 1 / shape,
                lower.tail = lower.tail, log.p = log.p)
}
qllogis <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE) {
  if (any(!is.finite(c(shape, scale))) || shape <= 0 || scale <= 0) return(rep(NaN, length(p)))
  exp(log(scale) + stats::qlogis(p, lower.tail = lower.tail, log.p = log.p) / shape)
}

dburr <- function(x, shape1, shape2, scale, log = FALSE) {
  .positive_density(x, c(shape1, shape2, scale), function(x) {
    z <- log(x) - log(scale)
    log(shape1) + log(shape2) - log(scale) + (shape1 - 1) * z -
      (shape2 + 1) * .softplus(shape1 * z)
  }, if (shape1 == 1) log(shape2) - log(scale) else if (shape1 < 1) Inf else -Inf, log)
}
pburr <- function(q, shape1, shape2, scale, lower.tail = TRUE, log.p = FALSE) {
  if (any(!is.finite(c(shape1, shape2, scale))) || any(c(shape1, shape2, scale) <= 0)) return(rep(NaN, length(q)))
  lh <- log(shape2) + .log_softplus(shape1 * (log(pmax(q, 0)) - log(scale)))
  .from_log_hazard(lh, lower.tail, log.p)
}
qburr <- function(p, shape1, shape2, scale, lower.tail = TRUE, log.p = FALSE) {
  if (any(!is.finite(c(shape1, shape2, scale))) || any(c(shape1, shape2, scale) <= 0)) return(rep(NaN, length(p)))
  exp(log(scale) + .log_expm1(stats::qexp(p, rate = shape2,
      lower.tail = lower.tail, log.p = log.p)) / shape1)
}

dgengamma <- function(x, shape, scale, k, log = FALSE) {
  .positive_density(x, c(shape, scale, k), function(x) {
    z <- log(x) - log(scale)
    log(k) - log(scale) - lgamma(shape) + (k * shape - 1) * z - exp(k * z)
  }, if (k * shape == 1) log(k) - log(scale) - lgamma(shape) else if (k * shape < 1) Inf else -Inf, log)
}
pgengamma <- function(q, shape, scale, k, lower.tail = TRUE, log.p = FALSE) {
  if (any(!is.finite(c(shape, scale, k))) || any(c(shape, scale, k) <= 0)) return(rep(NaN, length(q)))
  stats::pgamma(exp(k * (log(pmax(q, 0)) - log(scale))), shape = shape,
                lower.tail = lower.tail, log.p = log.p)
}
qgengamma <- function(p, shape, scale, k, lower.tail = TRUE, log.p = FALSE) {
  if (any(!is.finite(c(shape, scale, k))) || any(c(shape, scale, k) <= 0)) return(rep(NaN, length(p)))
  exp(log(scale) + log(stats::qgamma(p, shape = shape,
      lower.tail = lower.tail, log.p = log.p)) / k)
}
