test_that("make_lod_data standardizes input", {
  dat <- make_lod_data(
    nhanes_pah$X_cens,
    nhanes_pah$ind,
    nhanes_pah$cutoffs
  )

  expect_s3_class(dat, "copmi_lod_data")
  expect_equal(dim(dat$X_cens), dim(nhanes_pah$X_cens))
  expect_true(all(dat$ind %in% c(0L, 1L)))
  expect_true(all(dat$X_cens[dat$ind == 0L] <= rep(dat$cutoffs, each = nrow(dat$X_cens))[dat$ind == 0L]))
})

test_that("copula_em_impute returns stable completed matrices", {
  fit <- copula_em_impute(
    nhanes_pah,
    m = 2,
    margin_mode = "normal",
    max_iter = 100,
    seed = 42
  )

  expect_s3_class(fit, "copmi_mi")
  expect_length(fit$imp_list, 2)
  expect_equal(dim(complete(fit, 1)), dim(nhanes_pah$X_cens))
  expect_equal(dim(complete(fit, "mean")), dim(nhanes_pah$X_cens))

  observed <- nhanes_pah$ind == 1L
  censored <- nhanes_pah$ind == 0L
  imp1 <- complete(fit, 1)
  expect_equal(imp1[observed], nhanes_pah$X_cens[observed])
  cutoff_mat <- matrix(nhanes_pah$cutoffs, nrow = nrow(imp1), ncol = ncol(imp1), byrow = TRUE)
  expect_true(all(imp1[censored] < cutoff_mat[censored]))

  dx <- diagnostics(fit)
  expect_true(is.matrix(dx$Sigma_hat))
  expect_equal(dim(dx$Sigma_hat), c(ncol(imp1), ncol(imp1)))
  expect_true(all(is.finite(dx$Sigma_hat)))
  expect_true(is.data.frame(dx$fallback_records))
})

test_that("copula_em_impute is reproducible with a fixed seed", {
  fit1 <- copula_em_impute(nhanes_pah, m = 1, margin_mode = "normal", max_iter = 100, seed = 99)
  fit2 <- copula_em_impute(nhanes_pah, m = 1, margin_mode = "normal", max_iter = 100, seed = 99)

  expect_equal(complete(fit1, 1), complete(fit2, 1))
})

test_that("default sd_shift mode runs on example data", {
  fit <- copula_em_impute(
    nhanes_pah,
    m = 1,
    margin_mode = "sd_shift",
    max_iter = 100,
    seed = 7,
    margin_candidates = c("norm", "logis")
  )

  expect_s3_class(fit, "copmi_mi")
  expect_equal(dim(complete(fit, 1)), dim(nhanes_pah$X_cens))
  expect_length(fit$family_selected, ncol(nhanes_pah$X_cens))
})
