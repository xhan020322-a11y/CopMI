gibbs_case <- function() {
  Z <- rbind(c(-0.5, -1, 0.4), c(-0.5, 0.2, -1), c(0.2, 0.3, 0.4),
             c(-0.5, -1, -1), c(-0.5, -1, 0.6))
  ind <- rbind(c(0L, 0L, 1L), c(0L, 1L, 0L), c(1L, 1L, 1L),
               c(0L, 0L, 0L), c(0L, 0L, 1L))
  dimnames(Z) <- dimnames(ind) <- list(paste0("row", 1:5), paste0("X", 1:3))
  list(Z_full = Z, IND = ind, z_lod = c(-0.5, -1, -1),
       Sigma_hat = diag(0.35, 3) + 0.65)
}

test_that("Gibbs retains the requested sweeps and preserves observed cells", {
  dat <- gibbs_case()
  set.seed(321)
  out <- do.call(.draw_copula_engine, c(dat, list(M = 5, gibbs_burn = 200, gibbs_thin = 50)))
  dx <- out$sampling
  expect_identical(dx$kept_sweeps, seq.int(250L, 450L, 50L))
  expect_identical(dx$total_sweeps, 450L)
  expect_true(dx$completed)
  expect_equal(dx$n_patterns, 3)
  expect_equal(dx$n_rows_with_censoring, 4)
  expect_equal(dx$n_scalar_updates, 450 * sum(dat$IND == 0L))
  expect_length(dx$trace_mean_censored_z, 450)
  expect_true(all(is.finite(dx$trace_mean_abs_update)))
  expect_true(all(dx$trace_mean_abs_update >= 0))
  upper <- matrix(dat$z_lod, 5, 3, byrow = TRUE)
  expect_length(out$Z_imp_list, 5)
  for (Z in out$Z_imp_list) {
    expect_identical(dimnames(Z), dimnames(dat$Z_full))
    expect_identical(Z[dat$IND == 1L], dat$Z_full[dat$IND == 1L])
    expect_true(all(is.finite(Z)))
    expect_true(all(Z[dat$IND == 0L] < upper[dat$IND == 0L]))
  }
  expect_false(identical(out$Z_imp_list[[1]], out$Z_imp_list[[2]]))
  expect_equal(dx$trace_mean_censored_z[dx$kept_sweeps],
    vapply(out$Z_imp_list, function(Z) mean(Z[dat$IND == 0L]), numeric(1)))
})

test_that("retention uses one continuous chain", {
  dat <- gibbs_case()
  set.seed(47)
  many <- do.call(.draw_copula_engine, c(dat, list(M = 3, gibbs_burn = 2, gibbs_thin = 4)))
  for (k in 1:3) {
    set.seed(47)
    one <- do.call(.draw_copula_engine,
      c(dat, list(M = 1, gibbs_burn = 2 + (k - 1) * 4, gibbs_thin = 4)))
    expect_identical(one$Z_imp_list[[1]], many$Z_imp_list[[k]])
  }
  set.seed(47)
  zero <- do.call(.draw_copula_engine, c(dat, list(M = 2, gibbs_burn = 0, gibbs_thin = 1)))
  expect_identical(zero$sampling$kept_sweeps, 1:2)
})

test_that("truncated draws have the expected moments and handle extreme tails", {
  set.seed(91)
  x <- .rnorm_upper_trunc(rep(0, 50000), 1, 0)
  expect_equal(mean(x), -sqrt(2 / pi), tolerance = 0.02)
  expect_equal(stats::var(x), 1 - 2 / pi, tolerance = 0.02)
  for (upper in c(-40, -8, 0, 8)) {
    x <- .rnorm_upper_trunc(rep(0, 1000), 1, upper)
    expect_true(all(is.finite(x) & x < upper))
  }
  expect_error(.rnorm_upper_trunc(0, 0, 1), "Invalid")
  expect_error(.rnorm_upper_trunc(NA_real_, 1, 1), "Invalid")
  expect_error(.rnorm_upper_trunc(0, 1, -Inf), "Invalid")
})

test_that("correlated Gibbs draws agree with truncated-normal moments", {
  sigma <- matrix(c(1, 0.8, 0.8, 1), 2)
  target <- tmvtnorm::mtmvnorm(mean = c(0, 0), sigma = sigma,
    lower = c(-Inf, -Inf), upper = c(0, 0), doComputeVariance = TRUE)
  set.seed(113)
  out <- .draw_copula_engine(matrix(0, 5000, 2), matrix(0L, 5000, 2),
                             c(0, 0), sigma, 1, 200, 100)$Z_imp_list[[1]]
  expect_equal(colMeans(out), as.vector(target$tmean), tolerance = 0.035)
  expect_equal(unname(stats::cov(out)), unname(target$tvar), tolerance = 0.035)
})

test_that("no censoring needs no sweeps and a single variable retains its shape", {
  dat <- gibbs_case()
  dat$IND[] <- 1L
  set.seed(89)
  before <- .Random.seed
  out <- do.call(.draw_copula_engine, c(dat, list(M = 2, gibbs_burn = 200, gibbs_thin = 50)))
  expect_identical(.Random.seed, before)
  expect_identical(out$Z_imp_list, rep(list(dat$Z_full), 2))
  expect_equal(out$sampling$total_sweeps, 0)
  expect_length(out$sampling$kept_sweeps, 0)
  expect_length(out$sampling$trace_mean_censored_z, 0)
  expect_length(out$sampling$trace_mean_abs_update, 0)
  out <- .draw_copula_engine(matrix(c(0, 1), 2, 1), matrix(c(0L, 1L), 2, 1),
                             0, matrix(1), 2, 0, 1)
  expect_equal(dim(out$Z_imp_list[[1]]), c(2, 1))
  expect_equal(out$Z_imp_list[[1]][2, 1], 1)
})

test_that("invalid schedules and conditional covariances fail explicitly", {
  for (bad in list(-1, 0.5, NA_real_, Inf, c(1, 2))) {
    expect_error(copula_em_impute(nhanes_pah, gibbs_burn = bad), "gibbs_burn")
    expect_error(copula_em_impute(nhanes_pah, gibbs_thin = bad), "gibbs_thin")
  }
  expect_error(copula_em_impute(nhanes_pah, gibbs_thin = 0), "gibbs_thin")
  expect_error(copula_em_impute(nhanes_pah, m = 2, gibbs_thin = .Machine$integer.max),
               "integer range")
  expect_error(copula_em_impute(nhanes_pah, ag = "gibbs"), "Unused argument")
  expect_error(copula_em_impute(nhanes_pah, thinning = 2), "Unused argument")
  dat <- gibbs_case()
  for (sigma in list(matrix(1, 3, 3), diag(-1, 3), matrix(NA_real_, 3, 3))) {
    dat$Sigma_hat <- sigma
    expect_error(do.call(.draw_copula_engine,
      c(dat, list(M = 1, gibbs_burn = 0, gibbs_thin = 1))), "Gibbs initialization")
  }
})

test_that("a sweep failure reports its coordinate without returning partial draws", {
  original <- .rnorm_upper_trunc
  calls <- 0L
  testthat::local_mocked_bindings(.rnorm_upper_trunc = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("forced update failure")
    original(...)
  })
  expect_error(.draw_copula_engine(matrix(0), matrix(0L), 0, matrix(1), 1, 0, 1),
               "Gibbs sweep 1, pattern TRUE, variable 1.*forced update failure")
})


test_that("both public imputation APIs default to burn 200 and lag 50", {
  dat <- make_lod_data(nhanes_pah$X_cens[1:100, ],
                       nhanes_pah$ind[1:100, ], nhanes_pah$cutoffs, scale = "log")
  fit <- copula_em_impute(dat, margin_mode = "normal", seed = 1)
  staged <- copmi_impute(fit$model)
  for (result in list(fit, staged)) {
    dx <- result$extra$sampling
    expect_identical(dx$burn_in, 200L)
    expect_identical(dx$thin, 50L)
    expect_identical(dx$kept_sweeps, seq.int(250L, 450L, 50L))
    expect_identical(dx$total_sweeps, 450L)
    expect_length(result$imp_list, 5L)
  }
  expect_identical(staged$imp_list, fit$imp_list)
  custom <- copmi_impute(fit$model, gibbs_thin = 100L)
  expect_identical(custom$extra$sampling$kept_sweeps, seq.int(300L, 700L, 100L))
})
