test_that("montecarlo_integrate input validation works", {
  expect_error(montecarlo_integrate(123, 0, 1, 1000), "F Must be a function")
  expect_error(montecarlo_integrate(function(x) sum(x), "a", 1, 1000), "Lower Bound Must be Numeric")
  expect_error(montecarlo_integrate(function(x) sum(x), 0, "b", 1000), "Upper Bound Must be Numeric")
  expect_error(montecarlo_integrate(function(x) sum(x), 0, 1, "1000"), "Number of Samples Must be Numeric")
  expect_error(montecarlo_integrate(function(x) sum(x), 0, 1, 1000, "dim"), "Dimension  Must be Numeric")
  expect_error(montecarlo_integrate(function(x) sum(x), 0, 1, 3000000000), "Number of Samples Must Be Less Than 2147483648")
})

test_that("montecarlo_integrate works for 1D function", {
  f <- function(x) x[1]^2
  set.seed(123)
  est <- montecarlo_integrate(f, 0, 1, n_samples = 1e5, dim = 1)
  expect_true(abs(est - 1/3) < 0.01)
})

test_that("montecarlo_integrate works for 2D function", {
  f <- function(x) x[1] + x[2]
  set.seed(123)
  est <- montecarlo_integrate(f, 0, 1, n_samples = 1e5, dim = 2)
  expect_true(abs(est - 1) < 0.02)
})

test_that("montecarlo_integrate handles negative bounds", {
  f <- function(x) x[1]
  set.seed(123)
  est <- montecarlo_integrate(f, -1, 1, n_samples = 1e5, dim = 1)
  expect_true(abs(est) < 0.01)
})
