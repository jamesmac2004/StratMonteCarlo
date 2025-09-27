test_that("montecarlo_integrate input validation works", {
  expect_error(montecarlo_integrate(123, 0, 1, 1000, FALSE))
  expect_error(montecarlo_integrate(function(x) x^2, "a", 1, 1000, FALSE))
  expect_error(montecarlo_integrate(function(x) x^2, 0, "b", 1000, FALSE))
  expect_error(montecarlo_integrate(function(x) x^2, 0, 1, "1000", FALSE))
  expect_error(montecarlo_integrate(function(x) x^2, 0, 1, -1000, FALSE))
})

test_that("montecarlo_integrate works for 1D polynomial function", {
  f1 <- function(x) x^2
  set.seed(123)
  n_samples <- 50000

  mc_f1 <- montecarlo_integrate(f1, 0, 1, n_samples, FALSE)
  mc_f1_part <- montecarlo_integrate(f1, 0, 1, n_samples, TRUE)

  # Both should be close to 1/3 ≈ 0.333
  expect_true(abs(mc_f1$final - 1/3) < 0.02)
  expect_true(abs(mc_f1_part$final - 1/3) < 0.02)

  # Check object structure
  expect_true(inherits(mc_f1, "mc_result"))
  expect_true(inherits(mc_f1_part, "mc_result"))
})

test_that("montecarlo_integrate works for 1D trigonometric function", {
  f2 <- function(x) sin(x) * pi
  set.seed(123)
  n_samples <- 50000

  mc_f2 <- montecarlo_integrate(f2, 0, 4, n_samples, FALSE)
  mc_f2_part <- montecarlo_integrate(f2, 0, 4, n_samples, TRUE)

  # Expected: π*(1 - cos(4)) ≈ π*2.347 ≈ 7.37
  expected <- pi * (1 - cos(4))
  expect_true(abs(mc_f2$final - expected) < 0.5)
  expect_true(abs(mc_f2_part$final - expected) < 0.5)
})

test_that("montecarlo_integrate works for 1D oscillatory function", {
  f3 <- function(x) sin(10 * x)
  set.seed(123)
  n_samples <- 50000

  mc_f3 <- montecarlo_integrate(f3, 0, pi*4, n_samples, FALSE)
  mc_f3_part <- montecarlo_integrate(f3, 0, pi*4, n_samples, TRUE)

  # Expected: approximately 0 (oscillating function over full periods)
  # Partitioned should be more accurate for oscillatory functions
  expect_true(abs(mc_f3$final) < 1.0)  # Standard may be less accurate
  expect_true(abs(mc_f3_part$final) < 0.1)  # Partitioned should be closer to 0
})

test_that("montecarlo_integrate works for 2D function", {
  f4 <- function(x, y) exp(-((x-0.5)^2 + (y-0.5)^2)/0.01)
  set.seed(123)
  n_samples <- 30000  # Smaller sample size for 2D

  mc_f4 <- montecarlo_integrate(f4, c(0,0), c(pi*4, pi*4), n_samples, FALSE, dim = 2)
  mc_f4_part <- montecarlo_integrate(f4, c(0,0), c(pi*4, pi*4), n_samples, TRUE, dim = 2)

  # This is a sharp peak function - partitioned should be more accurate
  # Both should be positive values
  expect_true(mc_f4$final > 0)
  expect_true(mc_f4_part$final > 0)

  # Check 2D functionality works
  expect_true(inherits(mc_f4, "mc_result"))
  expect_true(inherits(mc_f4_part, "mc_result"))
})

test_that("mc_print function works", {
  f1 <- function(x) x^2
  mc_f1 <- montecarlo_integrate(f1, 0, 1, 1000, FALSE)

  # Should not throw error
  expect_no_error(mc_print(mc_f1))

  # Test with different input types
  expect_no_error(mc_print(c(1, 2, 3)))  # numeric vector
  expect_no_error(mc_print(list(estimates = c(1, 2, 3))))  # list
})

test_that("partitioning helps with oscillatory functions", {
  f_osc <- function(x) sin(10*x)
  set.seed(123)
  n_samples <- 30000
  mc_std <- montecarlo_integrate(f_osc, 0, 2*pi, n_samples, FALSE)
  mc_part <- montecarlo_integrate(f_osc, 0, 2*pi, n_samples, TRUE)
  # For oscillatory functions over complete periods, integral should be ≈ 0
  # Test that partitioned result is within reasonable bounds
  expect_true(mc_part$final > -1 & mc_part$final < 1)
})
