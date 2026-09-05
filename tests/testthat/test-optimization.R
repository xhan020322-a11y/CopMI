test_that("a failed primary triggers real Nelder-Mead on the same likelihood", {
  original <- .run_optimizer
  testthat::local_mocked_bindings(.run_optimizer = function(start, objective, method, maxit = 2000L) {
    if (method == "L-BFGS-B") return(list(valid = FALSE, message = "forced primary failure", warnings = ""))
    original(start, objective, method, maxit)
  })
  d <- nhanes_pah
  out <- copmi_fit_margin(d$X_cens[, 1], d$ind[, 1], d$cutoffs[1], "norm")
  expect_identical(out$optimization$optimizer, c("L-BFGS-B", "Nelder-Mead"))
  expect_identical(out$optimizer, "Nelder-Mead")
  expect_identical(out$optimization$converged, c(FALSE, TRUE))
  obs <- d$ind[, 1] == 1L
  mu <- out$parameters["mean"]
  sd <- out$parameters["sd"]
  ll <- sum(stats::dnorm(d$X_cens[obs, 1], mu, sd, log = TRUE)) +
    sum(!obs) * stats::pnorm(d$cutoffs[1], mu, sd, log.p = TRUE)
  expect_equal(out$loglik, as.numeric(ll), tolerance = 1e-8)
  expect_equal(out$bic, as.numeric(-2 * ll + 2 * log(length(obs))), tolerance = 1e-8)
})

test_that("optimizer failures and penalty plateaus cannot become fitted models", {
  penalty <- .run_optimizer(c(a = 1, b = 1), function(x) 1e100, "L-BFGS-B")
  expect_false(penalty$valid)
  failed <- .run_optimizer(c(a = 1, b = 1), function(x) stop("test failure"), "L-BFGS-B")
  expect_false(failed$valid)
  expect_match(failed$message, "test failure")
  testthat::local_mocked_bindings(.run_optimizer = function(...) list(valid = FALSE,
    message = "forced optimizer failure", warnings = ""))
  d <- nhanes_pah
  expect_error(copmi_fit_margin(d$X_cens[, 1], d$ind[, 1], d$cutoffs[1], "norm"), "All candidate")
})

test_that("successful primary skips backup and incompatible candidates are recorded", {
  d <- nhanes_pah
  x <- d$X_cens[, 1] - 10
  out <- copmi_fit_margin(x, d$ind[, 1], d$cutoffs[1] - 10, c("norm", "exp"))
  expect_identical(out$family, "norm")
  expect_identical(out$optimization$optimizer, "L-BFGS-B")
  expect_false(out$candidate_table$converged[2])
  expect_match(out$candidate_table$message[2], "positive")
})

test_that("all families share the correct censored likelihood and inverse CDF", {
  shifted <- copmi_shift(nhanes_pah)$dat
  j <- 1
  obs <- shifted$ind[, j] == 1L
  for (family in .margin_families()) {
    fit <- copmi_fit_margin(shifted$X_cens[, j], shifted$ind[, j], shifted$cutoffs[j], family)
    dist <- .distribution(family)
    pars <- as.list(fit$parameters)
    ll <- sum(do.call(dist$d, c(list(x = shifted$X_cens[obs, j], log = TRUE), pars))) +
      sum(!obs) * do.call(dist$p, c(list(q = shifted$cutoffs[j], log.p = TRUE), pars))
    expect_equal(fit$loglik, as.numeric(ll), tolerance = 1e-7, info = family)
    probabilities <- c(1e-8, 0.1, 0.5, 0.9, 1 - 1e-8)
    quantiles <- do.call(dist$q, c(list(p = probabilities), pars))
    expect_equal(do.call(dist$p, c(list(q = quantiles), pars)), probabilities,
                 tolerance = 1e-6, info = family)
  }
})

test_that("pairwise initialization retries rather than substituting observed correlation", {
  original <- .run_optimizer
  testthat::local_mocked_bindings(.run_optimizer = function(start, objective, method, maxit = 2000L) {
    if (method == "L-BFGS-B") return(list(valid = FALSE, message = "forced primary failure", warnings = ""))
    original(start, objective, method, maxit)
  })
  d <- nhanes_pah
  z <- copmi_transform(copmi_fit_margins(d, margin_mode = "normal"))
  pair <- estimate_lyles_pair_cor(z$Z[, 1], z$Z[, 2], z$ind[, 1], z$ind[, 2], z$z_lod[1], z$z_lod[2])
  expect_identical(pair$optimizer, "Nelder-Mead")
  expect_equal(nrow(pair$attempts), 2)
  expect_error(estimate_lyles_pair_cor(1:6, 6:1, rep(1, 6), rep(1, 6), NA_real_, NA_real_, min_obs = 10), "Too few")
})
