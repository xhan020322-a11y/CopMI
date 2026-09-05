test_that("indicators are validated before coercion and labels cannot misalign", {
  d <- nhanes_pah
  bad <- d$ind
  bad[1, 1] <- 0.5
  expect_error(make_lod_data(d$X_cens, bad, d$cutoffs), "only 0/1")
  bad[1, 1] <- NA
  expect_error(make_lod_data(d$X_cens, bad, d$cutoffs), "only 0/1")
  expect_error(make_lod_data(d$X_cens, d$ind[, ncol(d$ind):1], d$cutoffs), "dimnames")
  expect_error(make_lod_data(d$X_cens, d$ind, stats::setNames(d$cutoffs, rev(colnames(d$X_cens)))), "Named cutoffs")
  expect_error(make_lod_data(matrix(character(6), 2), matrix(1, 2, 3), rep(NA_real_, 3)), "numeric")
  expect_error(make_lod_data(matrix(numeric(), 0, 3), matrix(1, 0, 3), rep(0, 3)), "at least")
})

test_that("only censored cells may be missing and observed cutoffs are consistent", {
  d <- nhanes_pah
  x <- d$X_cens
  x[d$ind == 0L] <- NA_real_
  expect_equal(make_lod_data(x, d$ind, d$cutoffs)$X_cens, d$X_cens)
  x[which(d$ind == 1L)[1]] <- NA_real_
  expect_error(make_lod_data(x, d$ind, d$cutoffs), "Observed")
  expect_error(make_lod_data(d$X_cens, d$ind, rep(NA_real_, ncol(d$ind))), "finite cutoff")
  x <- d$X_cens
  x[which(d$ind[, 1] == 1L)[1], 1] <- d$cutoffs[1] - 1
  expect_error(make_lod_data(x, d$ind, d$cutoffs), "below the cutoff")
})

test_that("unsupported parameters and unidentifiable margins fail early", {
  d <- nhanes_pah
  for (m in list(0, -1, 1.5, c(1, 2), NA_real_, Inf)) {
    expect_error(copula_em_impute(d, m = m), "m must")
  }
  expect_error(copula_em_impute(d, max_iter = -1), "max_iter")
  expect_error(copula_em_impute(d, tol = 0), "tol")
  expect_error(copula_em_impute(d, thinning = 0), "thinning")
  expect_error(copula_em_impute(d, seeds = 10), "Unused")
  expect_error(copula_em_impute(d, ind = d$ind), "again")
  expect_error(copula_em_impute(d, lyles_control = list(max_it = 5)), "lyles_control")
  expect_error(copmi_fit_margins(d, margin_candidates = "unknown"), "candidates")
  expect_error(copmi_fit_margins(d, optim_methods = "BFGS"), "optim_methods")
  expect_error(copmi_fit_margin(1, 1, NA_real_, "norm"), "two")
  expect_error(copmi_fit_margin(rep(1, 8), rep(1, 8), NA_real_, "norm"), "variation")
  d$X_cens[which(d$ind == 1L)[1]] <- 0
  d$cutoffs[] <- NA_real_
  d$ind[] <- 1L
  expect_error(copmi_fit_margins(d, input_scale = "raw"), "positive")
})
