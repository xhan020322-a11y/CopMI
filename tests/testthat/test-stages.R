test_that("staged and one-call pipelines agree in both modes", {
  for (mode in c("normal", "sd_shift")) {
    margins <- copmi_fit_margins(nhanes_pah, margin_mode = mode,
                                  margin_candidates = c("norm", "logis"))
    latent <- copmi_transform(margins)
    model <- copmi_fit_copula(latent, max_iter = 100, seed = 42)
    staged <- copmi_impute(model, m = 2)
    direct <- copula_em_impute(nhanes_pah, margin_mode = mode,
      margin_candidates = c("norm", "logis"), max_iter = 100, m = 2, seed = 42)
    expect_identical(complete(staged), complete(direct))
    expect_identical(staged$Sigma_hat, direct$Sigma_hat)
    expect_s3_class(margins, "copmi_margins")
    expect_s3_class(latent, "copmi_latent")
    expect_s3_class(model, "copmi_copula")
    expect_equal(length(model$em_change_history), model$n_iter)
    expect_true(all(eigen(model$Sigma_hat, symmetric = TRUE)$values > 0))
    expect_equal(unname(diag(model$Sigma_hat)), rep(1, ncol(nhanes_pah$X_cens)))
  }
})

test_that("raw and log input paths return the documented scales", {
  d <- nhanes_pah
  raw <- make_lod_data(exp(d$X_cens), d$ind, exp(d$cutoffs))
  for (mode in c("normal", "select")) {
    raw_fit <- copula_em_impute(raw, input_scale = "raw", margin_mode = mode,
      margin_candidates = c("norm", "logis"), m = 1, max_iter = 100, seed = 12)
    log_fit <- copula_em_impute(d, input_scale = "log", margin_mode = mode,
      margin_candidates = c("norm", "logis"), m = 1, max_iter = 100, seed = 12)
    x <- complete(raw_fit, 1)
    expect_equal(as.vector(x), as.vector(exp(complete(log_fit, 1))), tolerance = 1e-6)
    expect_identical(x[d$ind == 1], raw$X_cens[d$ind == 1])
    expect_true(all(x > 0))
    cut <- matrix(raw$cutoffs, nrow(d$X_cens), ncol(d$X_cens), byrow = TRUE)
    expect_true(all(x[d$ind == 0] < cut[d$ind == 0]))
  }
})

test_that("normal assumption bypasses candidate selection and shift", {
  m <- copmi_fit_margins(nhanes_pah, margin_mode = "normal",
                         margin_candidates = "not_used", shift_k = NA_real_)
  expect_true(all(m$family_selected == "norm"))
  expect_null(m$shift)
  expect_equal(nrow(m$candidate_table), ncol(nhanes_pah$X_cens))
  expect_equal(m$work_data$X_cens, nhanes_pah$X_cens)
})

test_that("single-column, no-censoring, and all-censored-row cases are supported", {
  d <- nhanes_pah
  one <- make_lod_data(d$X_cens[, 1, drop = FALSE], d$ind[, 1, drop = FALSE], d$cutoffs[1])
  fit <- copula_em_impute(one, margin_mode = "normal", m = 1, max_iter = 100)
  expect_equal(dim(complete(fit, "mean")), c(nrow(d$X_cens), 1L))
  expect_equal(as.numeric(fit$Sigma_hat), 1)
  observed_rows <- rowSums(d$ind == 0L) == 0L
  x <- d$X_cens[observed_rows, , drop = FALSE]
  full <- make_lod_data(x, matrix(1, nrow(x), ncol(x)), rep(NA_real_, ncol(x)))
  fit <- copula_em_impute(full, margin_mode = "normal", m = 2, max_iter = 100)
  expect_equal(as.vector(complete(fit, 1)), as.vector(x))
  cols <- which(colSums(d$ind == 0L) > 0L)
  ind <- d$ind[, cols]
  ind[1, ] <- 0L
  allrow <- make_lod_data(d$X_cens[, cols], ind, d$cutoffs[cols])
  fit <- copula_em_impute(allrow, margin_mode = "normal", m = 1, max_iter = 100)
  expect_true(all(complete(fit, 1)[1, ] < d$cutoffs[cols]))
  init <- copula_em_impute(d, margin_mode = "normal", m = 1, max_iter = 0, allow_unconverged = TRUE)
  expect_false(init$converged)
  expect_equal(init$n_iter, 0L)
  expect_length(init$model$em_change_history, 0)
})

test_that("names, extraction, and RNG state survive the workflow", {
  d <- nhanes_pah
  rownames(d$X_cens) <- rownames(d$ind) <- paste0("sample_", seq_len(nrow(d$X_cens)))
  set.seed(345)
  before <- .Random.seed
  fit <- copula_em_impute(d, margin_mode = "normal", m = 1, max_iter = 100, seed = 7)
  expect_identical(.Random.seed, before)
  expect_identical(dimnames(complete(fit, 1)), dimnames(d$X_cens))
  expect_equal(as.vector(complete(fit, "mean")), as.vector(complete(fit, 1)))
  expect_identical(names(complete(fit, "long")), c(".imp", ".id", colnames(d$X_cens)))
  expect_equal(nrow(complete(fit, "long")), nrow(d$X_cens))
  expect_error(complete(fit, 1.2), "action")
  expect_error(complete(fit, c(1, 2)), "action")
  expect_identical(complete(copmi_impute(fit$model, m = 1)), complete(fit))
  expect_false(identical(complete(copmi_impute(fit$model, m = 1, seed = 987)), complete(fit)))
})
