test_that("positive families reject negative observations without changing data", {
  x <- c(-1, 1:9)
  expect_error(copmi_fit_margin(x, rep(1, 10), NA_real_, "exp"), "positive observations")
  expect_identical(x, c(-1, 1:9))
  tiny <- 1e-20 * (1:10)
  fit <- copmi_fit_margin(tiny, rep(1, 10), NA_real_, "exp")
  expect_equal(unname(fit$parameters["rate"]) * mean(tiny), 1, tolerance = 1e-5)
  normal <- copmi_fit_margin(tiny, rep(1, 10), NA_real_, "norm")
  expect_equal(unname(normal$parameters["mean"]) / mean(tiny), 1, tolerance = 1e-4)
  expect_equal(unname(normal$parameters["sd"]) / sqrt(mean((tiny - mean(tiny))^2)), 1, tolerance = 1e-4)
})

test_that("a failed E-step stops instead of substituting untruncated moments", {
  testthat::local_mocked_bindings(.truncated_moments = function(...) stop("forced moments failure"))
  expect_error(copula_em_impute(nhanes_pah, margin_mode = "normal"), "EM iteration.*row.*forced moments failure")
})

test_that("a failed sample stops instead of returning its fixed starting value", {
  model <- copmi_fit_copula(copmi_transform(copmi_fit_margins(nhanes_pah, margin_mode = "normal")))
  testthat::local_mocked_bindings(.truncated_draw = function(...) stop("forced draw failure"))
  expect_error(copmi_impute(model, m = 3), "Imputation.*row.*forced draw failure")
})

test_that("inverse transformation preserves real-support shifted normal tails", {
  m <- copmi_fit_margins(nhanes_pah, margin_mode = "select", margin_candidates = "norm")
  latent <- copmi_transform(m)
  latent$Z[latent$ind == 0] <- -10
  x <- .inverse_latent(latent$Z, m)
  for (j in seq_len(ncol(x))) {
    cens <- latent$ind[, j] == 0L
    expected <- m$fits[[j]]$parameters["mean"] - 10 * m$fits[[j]]$parameters["sd"] + m$shift$anchor[j]
    expect_equal(unname(x[cens, j]), rep(unname(expected), sum(cens)))
    expect_true(all(x[cens, j] < m$shift$anchor[j]))
  }
  fit <- list(family = "norm", parameters = c(mean = 0, sd = 1))
  expect_equal(.gaussian_scores(c(-40, 0, 40), fit), c(-40, 0, 40))
  expect_equal(.inverse_scores(c(-40, 0, 40), fit), c(-40, 0, 40))
})

test_that("invalid inverse values fail rather than becoming LOD replacements", {
  d <- make_lod_data(matrix(c(0, 1, 2), 3, 1), matrix(c(0, 1, 1), 3, 1), 0)
  for (bad in c(NA_real_, Inf, 0, 1)) {
    x <- d$X_cens
    x[1] <- bad
    expect_error(.check_completed(d, x), "nonfinite|below its cutoff")
  }
})

test_that("SD-shift does not borrow a different column's scale", {
  x <- cbind(constant = rep(1, 10), varying = 1:10)
  dat <- make_lod_data(x, matrix(1, 10, 2), c(NA_real_, NA_real_))
  expect_error(copmi_shift(dat), "varying observed values in every column")
  dat <- make_lod_data(x[, 2, drop = FALSE] * 1e-20, matrix(1, 10, 1), NA_real_)
  expect_equal(copmi_shift(dat)$sd, stats::sd((1:10) * 1e-20))
})

test_that("only finite initialization matrices may be projected", {
  bad <- matrix(c(1, NA, NA, 1), 2)
  expect_error(make_pd_cor(bad, project = TRUE), "finite square")
  indefinite <- matrix(c(1, .9, .9, .9, 1, -.9, .9, -.9, 1), 3)
  expect_error(make_pd_cor(indefinite), "not positive definite")
  good <- make_pd_cor(indefinite, project = TRUE)
  expect_gt(attr(good, "pd_adjustment"), 0)
  expect_true(all(eigen(good, symmetric = TRUE)$values > 0))
})

test_that("unconverged models require explicit exploratory use", {
  latent <- copmi_transform(copmi_fit_margins(nhanes_pah, margin_mode = "normal"))
  model <- copmi_fit_copula(latent, max_iter = 0)
  expect_error(copmi_impute(model), "has not converged")
  out <- copmi_impute(model, m = 1, allow_unconverged = TRUE)
  expect_false(out$converged)
  expect_error(copula_em_impute(nhanes_pah, margin_mode = "normal", max_iter = 0), "has not converged")
})
