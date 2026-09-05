test_that("custom distributions preserve support, endpoints and invalid-parameter semantics", {
  parameters <- list(lomax = list(shape = 2, scale = 3), llogis = list(shape = 2, scale = 3),
    burr = list(shape1 = 2, shape2 = 3, scale = 2), gengamma = list(shape = 2, scale = 3, k = 2))
  for (family in names(parameters)) {
    dist <- .distribution(family)
    pars <- parameters[[family]]
    expect_equal(do.call(dist$d, c(list(x = c(-1, Inf)), pars)), c(0, 0), info = family)
    expect_equal(do.call(dist$d, c(list(x = -1, log = TRUE), pars)), -Inf, info = family)
    expect_equal(do.call(dist$p, c(list(q = c(-1, 0, Inf)), pars)), c(0, 0, 1), info = family)
    expect_equal(do.call(dist$q, c(list(p = c(0, 1)), pars)), c(0, Inf), info = family)
    invalid <- pars
    invalid[[1]] <- NA_real_
    expect_true(is.nan(do.call(dist$d, c(list(x = 1), invalid))))
    area <- stats::integrate(function(x) do.call(dist$d, c(list(x = x), pars)), 0, Inf)$value
    expect_equal(area, 1, tolerance = 1e-6, info = family)
    q <- c(0.01, 1, 10)
    derivative <- (do.call(dist$p, c(list(q = q + 1e-6), pars)) -
      do.call(dist$p, c(list(q = q - 1e-6), pars))) / 2e-6
    expect_equal(do.call(dist$d, c(list(x = q), pars)), derivative, tolerance = 1e-6, info = family)
  }
})

test_that("log tails avoid artificial probability floors", {
  expect_equal(pburr(exp(-600), 2, 3, 1, log.p = TRUE), log(3) - 1200)
  expect_equal(plomax(exp(-600), 2, 1, log.p = TRUE), log(2) - 600)
  expect_equal(dllogis(exp(-400), 4, 1, log = TRUE), log(4) - 1200)
  lp <- c(-700, -30, -1)
  expect_equal(plomax(qlomax(lp, 2, 3, log.p = TRUE), 2, 3, log.p = TRUE), lp)
  expect_equal(pburr(qburr(lp, 2, 3, 2, log.p = TRUE), 2, 3, 2, log.p = TRUE), lp)
})
