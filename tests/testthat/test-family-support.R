support_case <- function() {
  parameters <- list(norm = c(mean = 0, sd = 1), logis = c(location = 0, scale = 1),
    lnorm = c(meanlog = 0, sdlog = 1), gamma = c(shape = 2, rate = 1),
    weibull = c(shape = 2, scale = 1), exp = c(rate = 1),
    invgauss = c(mean = 1, shape = 1), gengamma = c(shape = 2, scale = 1, k = 2),
    llogis = c(shape = 2, scale = 1), lomax = c(shape = 2, scale = 1),
    burr = c(shape1 = 2, shape2 = 2, scale = 1))
  p <- length(parameters)
  X <- matrix(rep(c(2, 3, 5), p), 3, p,
    dimnames = list(c("censored", "observed1", "observed2"), names(parameters)))
  dat <- make_lod_data(X, matrix(rep(c(0, 1, 1), p), 3, p), rep(2, p))
  analysis <- dat
  analysis$X_cens <- dat$X_cens - 5
  analysis$cutoffs <- dat$cutoffs - 5
  list(work_data = dat, analysis_data = analysis, input = analysis,
    input_scale = "log", shift = list(anchor = rep(-5, p)),
    fits = setNames(lapply(names(parameters), function(family)
      list(family = family, parameters = parameters[[family]])), names(parameters)))
}

test_that("the distribution registry separates real and positive support", {
  positive <- vapply(.margin_families(), function(family) .distribution(family)$positive, logical(1))
  expect_identical(names(positive)[!positive], c("norm", "logis"))
  expect_equal(sum(positive), 9)
  expect_error(.distribution("unknown"), "Unsupported")
})

test_that("inverse transformation preserves both support groups without substitution", {
  margins <- support_case()
  Z <- matrix(-3, 3, length(margins$fits))
  X <- .inverse_latent(Z, margins)
  shifted <- X[1, ] - margins$shift$anchor
  expected <- vapply(margins$fits, function(fit) as.numeric(.inverse_scores(-3, fit)), numeric(1))
  expect_equal(shifted, expected)
  expect_true(all(shifted[c("norm", "logis")] < 0))
  expect_true(all(shifted[setdiff(names(shifted), c("norm", "logis"))] > 0))
  expect_true(all(shifted < margins$work_data$cutoffs))
  expect_identical(X[2:3, ], margins$input$X_cens[2:3, ])
  expect_identical(dimnames(X), dimnames(margins$input$X_cens))
  expect_equal(unname(attr(X, "numerics")), c(0, 0, 0))
  expect_true(all(X[1, ] < 0))
  for (family in c("norm", "logis")) {
    expect_equal(as.numeric(.inverse_scores(0, margins$fits[[family]])), 0)
  }
})

test_that("per-column validation rejects invalid positive draws and any material overshoot", {
  margins <- support_case()
  dat <- margins$work_data
  positive <- vapply(margins$fits, function(fit) .distribution(fit$family)$positive, logical(1))
  X <- dat$X_cens
  X[1, ] <- .5
  X[1, c("norm", "logis")] <- c(-2, 0)
  expect_equal(.check_completed(dat, X, positive)[1, c("norm", "logis")], c(norm = -2, logis = 0))
  for (family in names(positive)[positive]) {
    for (value in c(-.1, 0)) {
      invalid <- X
      invalid[1, family] <- value
      expect_error(.check_completed(dat, invalid, positive), "invalid censored")
    }
  }
  for (family in names(positive)) {
    invalid <- X
    invalid[1, family] <- 2.1
    expect_error(.check_completed(dat, invalid, positive), "below its cutoff")
    invalid[1, family] <- Inf
    expect_error(.check_completed(dat, invalid, positive), "nonfinite")
  }
})

test_that("inverse output is checked on the fitted scale before undoing the shift", {
  margins <- support_case()
  Z <- matrix(-3, 3, length(margins$fits))
  inverse <- .inverse_scores
  testthat::local_mocked_bindings(.inverse_scores = function(z, fit) {
    if (fit$family == "gamma") return(rep(-.1, length(z)))
    inverse(z, fit)
  })
  expect_error(.inverse_latent(Z, margins), "invalid censored")
})

test_that("BIC-selected normal and logistic margins retain their negative shifted tails", {
  for (family in c("norm", "logis")) {
    quantile <- .distribution(family)$q
    X <- quantile(stats::ppoints(401))
    cutoff <- X[41]
    observed <- as.integer(X >= cutoff)
    dat <- make_lod_data(matrix(pmax(X, cutoff), ncol = 1),
      matrix(observed, ncol = 1), cutoff)
    margins <- copmi_fit_margins(dat, margin_mode = "sd_shift",
      margin_candidates = c("norm", "logis"))
    expect_identical(unname(margins$family_selected), family)
    latent <- copmi_transform(margins)
    latent$Z[observed == 0] <- -12
    completed <- .inverse_latent(latent$Z, margins)
    expect_true(all(completed[observed == 0] < margins$shift$anchor))
    expect_identical(as.numeric(completed[observed == 1]), X[observed == 1])
  }
})

test_that("original-scale output stays positive for either selected support group", {
  margins <- support_case()
  margins$input_scale <- "raw"
  margins$input$X_cens <- exp(margins$analysis_data$X_cens)
  margins$input$cutoffs <- exp(margins$analysis_data$cutoffs)
  X <- .inverse_latent(matrix(-3, 3, length(margins$fits)), margins)
  expect_true(all(is.finite(X) & X > 0))
  expect_true(all(X[1, ] < margins$input$cutoffs))
  expect_identical(X[2:3, ], margins$input$X_cens[2:3, ])
})
