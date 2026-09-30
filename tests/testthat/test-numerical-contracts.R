test_that("only roundoff-sized boundary discrepancies may move inward", {
  x <- .round_upper_bound(c(0.5, 1, 1 + .Machine$double.eps), 1)
  expect_identical(x[1], 0.5)
  expect_true(all(x < 1))
  expect_equal(unname(attr(x, "boundary")[1]), 2)
  expect_error(.round_upper_bound(1.001, 1), "floating-point")
  expect_error(.round_upper_bound(Inf, 1), "Nonfinite")
  expect_error(.round_upper_bound(NA_real_, 1), "Nonfinite")
})

test_that("a block Cholesky produces the Gaussian conditional covariance", {
  S <- matrix(c(1, .5, .3, .5, 1, .2, .3, .2, 1), 3)
  out <- .conditional_normal(c(2, 0, 0), 1, 2:3, S)
  expect_equal(out$mean, c(1, .6))
  expect_equal(out$sigma, S[2:3, 2:3] - tcrossprod(S[2:3, 1]))
  bad <- S
  bad[1, 2] <- .6
  expect_error(.conditional_normal(c(2, 0, 0), 1, 2:3, bad), "symmetric")
  expect_error(make_pd_cor(matrix(c(1, 2, 2, 1), 2)), "not positive definite")
})

test_that("quantile root fallback solves the same log probability without RNG", {
  fit <- list(family = "gamma", parameters = c(shape = 2, rate = 1))
  set.seed(410)
  state <- .Random.seed
  for (lower in c(TRUE, FALSE)) {
    x <- .solve_marginal_quantile(-30, fit, lower)
    expect_equal(pgamma(x, 2, lower.tail = lower, log.p = TRUE), -30, tolerance = 1e-8)
  }
  expect_identical(.Random.seed, state)
  expect_error(.solve_marginal_quantile(-1e6, fit, TRUE), "representable")
})

test_that("invalid moment output stops without projecting covariance eigenvalues", {
  testthat::local_mocked_bindings(mtmvnorm = function(...) list(tmean = c(-1, -1),
    tvar = matrix(c(1, 2, 2, 1), 2)), .package = "tmvtnorm")
  expect_error(.truncated_moments(c(0, 0), diag(2), c(-Inf, -Inf), c(0, 0)), "semidefinite")
})

test_that("a failed quantile is recovered and counted without resampling", {
  distribution <- .distribution
  testthat::local_mocked_bindings(.distribution = function(family) {
    out <- distribution(family)
    out$q <- function(p, ...) rep(NA_real_, length(p))
    out
  })
  set.seed(101)
  state <- .Random.seed
  fit <- list(family = "gamma", parameters = c(shape = 2, rate = 1))
  out <- .inverse_scores(c(-2, 2), fit)
  expect_equal(as.numeric(out), qgamma(pnorm(c(-2, 2)), 2), tolerance = 1e-8)
  expect_equal(attr(out, "quantile_fallbacks"), 2)
  expect_identical(.Random.seed, state)
})

test_that("seeded failures restore both present and absent caller RNG state", {
  set.seed(919)
  state <- .Random.seed
  expect_error(.with_rng(10, code = { runif(1); stop("forced") }), "forced")
  expect_identical(.Random.seed, state)
  rm(".Random.seed", envir = .GlobalEnv)
  expect_error(.with_rng(10, code = { runif(1); stop("forced") }), "forced")
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  assign(".Random.seed", state, envir = .GlobalEnv)
})
