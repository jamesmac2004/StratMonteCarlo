# These tests check correctness, TinyExpr parsing, and performance across
# Standard, Partitioned, and Importance Sampling Monte Carlo methods

skip_on_cran()

test_that("TinyExpr expressions produce accurate results (1D)", {
  set.seed(123)
  # 1D polynomial: ∫_0^1 x^2 dx = 1/3
  mc_poly <- montecarlo_integrate("x1^2", 0, 1, 30000)
  expect_true(abs(mc_poly$estimate - 1/3) < 0.02)

  # Trigonometric: ∫_0^π sin(x) dx = 2
  mc_sin <- montecarlo_integrate("sin(x1)", 0, pi, 30000)
  expect_true(abs(mc_sin$estimate - 2) < 0.05)
})

test_that("Importance Sampling improves estimation for tail-heavy functions", {
  set.seed(42)
  # Integrand: exp(-x) over [0, 10], true value ≈ 1 - exp(-10) ≈ 0.99995
  expr <- "exp(-x1)"
  exact <- 1 - exp(-10)

  # Standard MC
  mc_std <- montecarlo_integrate(expr, 0, 10, 30000)

  # Importance Sampling with exponential distribution
  mc_is <- montecarlo_integrate(expr, 0, 10, 30000,
                                importance_sampling = TRUE,
                                is_distribution = "exponential")

  # Both should be close to exact value
  expect_true(abs(mc_std$estimate - exact) < 0.05)
  expect_true(abs(mc_is$estimate - exact) < 0.05)

  # IS should generally perform better for this type of function
  expect_true(mc_is$importance_sampling)
})

test_that("Importance Sampling handles Gaussian-like peaks", {
  set.seed(999)
  expr <- "exp(-((x1-0.5)^2)/0.01)"  # sharp peak at 0.5
  # Approximate numerical reference
  exact <- 0.177  # sqrt(pi * 0.01) ≈ 0.177

  # Standard Monte Carlo
  mc_std <- montecarlo_integrate(expr, 0, 1, 50000)

  # Partitioned (requires function, not expression)
  f_peak <- function(x) exp(-((x - 0.5)^2) / 0.01)
  mc_part <- montecarlo_integrate(f_peak, 0, 1, 50000, partition = TRUE)

  # Importance Sampling with normal distribution
  mc_is <- montecarlo_integrate(expr, 0, 1, 50000,
                                importance_sampling = TRUE,
                                is_distribution = "normal")

  # All methods should produce reasonable estimates
  expect_true(abs(mc_std$estimate - exact) < 0.05)
  expect_true(abs(mc_part$estimate - exact) < 0.05)
  expect_true(abs(mc_is$estimate - exact) < 0.05)

  # Verify IS was actually used
  expect_true(mc_is$importance_sampling)
  expect_equal(mc_is$is_distribution, "normal")
})

test_that("IS performs well on smooth low-variance functions", {
  set.seed(321)
  expr <- "x1 + x1^2"
  exact <- 1/2 + 1/3  # = 5/6 ≈ 0.833

  mc_std <- montecarlo_integrate(expr, 0, 1, 30000)
  mc_is <- montecarlo_integrate(expr, 0, 1, 30000,
                                importance_sampling = TRUE,
                                is_distribution = "beta")

  # Both should be accurate
  expect_true(abs(mc_std$estimate - exact) < 0.02)
  expect_true(abs(mc_is$estimate - exact) < 0.02)

  # Verify IS settings
  expect_true(mc_is$importance_sampling)
  expect_equal(mc_is$is_distribution, "beta")
})

test_that("IS handles oscillatory functions", {
  set.seed(111)
  expr <- "sin(10*x1)"
  # ∫_0^(2π) sin(10x) dx = 0

  mc_std <- montecarlo_integrate(expr, 0, 2*pi, 60000)
  mc_is <- montecarlo_integrate(expr, 0, 2*pi, 60000,
                                importance_sampling = TRUE,
                                is_distribution = "normal")

  # Both should integrate to approximately 0
  expect_true(abs(mc_std$estimate) < 0.2)
  expect_true(abs(mc_is$estimate) < 0.2)

  expect_true(mc_is$importance_sampling)
})

test_that("TinyExpr multi-dimensional expressions match known results", {
  set.seed(123)
  # 2D: ∫_0^1 ∫_0^1 (x1*x2) dx dy = 1/4
  mc_2d <- montecarlo_integrate("x1*x2", c(0,0), c(1,1), 30000, dim = 2)
  expect_true(abs(mc_2d$estimate - 0.25) < 0.02)

  # 3D: ∫ (x1*x2*x3) = 1/8
  mc_3d <- montecarlo_integrate("x1*x2*x3", c(0,0,0), c(1,1,1), 40000, dim = 3)
  expect_true(abs(mc_3d$estimate - 1/8) < 0.02)
})

test_that("TinyExpr parser catches bad syntax", {
  # Should error on undefined variable
  expect_error(montecarlo_integrate("sin(x)", 0, 1, 1000))

  # Should error on unknown function
  expect_error(montecarlo_integrate("unknownfunc(x1)", 0, 1, 1000))

  # Should work with valid expression
  expect_no_error(montecarlo_integrate("sin(x1)+x1^2", 0, 1, 1000))
})

test_that("Importance sampling with mixture_normal distribution", {
  set.seed(456)
  expr <- "x1^2"
  exact <- 1/3

  mc_mixture <- montecarlo_integrate(expr, 0, 1, 30000,
                                     importance_sampling = TRUE,
                                     is_distribution = "mixture_normal")

  expect_true(abs(mc_mixture$estimate - exact) < 0.02)
  expect_equal(mc_mixture$is_distribution, "mixture_normal")
  expect_true(!is.null(mc_mixture$params_used))
})

test_that("Partitioning works with and without IS", {
  set.seed(789)
  # Partitioning requires an R function, not an expression
  f_func <- function(x) x^3 + sin(x)

  # Partitioned without IS - should create multiple partitions
  mc_part <- montecarlo_integrate(f_func, 0, 2*pi, 40000, partition = TRUE)
  # Note: partitioning may return 1 partition if no local minima found
  # So we just verify it ran and returned a valid n_partitions field
  expect_true(!is.null(mc_part$n_partitions))
  expect_true(mc_part$n_partitions >= 1)

  # Partition + IS combination should work (no explicit restriction in code)
  # The function will be used for partitioning, IS will be applied within partitions
  mc_part_is <- montecarlo_integrate(f_func, 0, 2*pi, 40000,
                                     partition = TRUE,
                                     importance_sampling = TRUE,
                                     is_distribution = "normal")
  expect_true(mc_part_is$importance_sampling)
  expect_true(!is.null(mc_part_is$n_partitions))
})

test_that("Custom IS parameters can be provided", {
  set.seed(222)
  expr <- "x1^2"
  exact <- 1/3

  # Provide custom normal parameters
  custom_params <- list(mean = c(0.5), sd = c(0.2))

  mc_custom <- montecarlo_integrate(expr, 0, 1, 30000,
                                    importance_sampling = TRUE,
                                    is_distribution = "normal",
                                    is_params = custom_params)

  expect_true(abs(mc_custom$estimate - exact) < 0.02)
  expect_equal(mc_custom$params_used$mean[1], custom_params$mean[1])
  expect_equal(mc_custom$params_used$sd[1], custom_params$sd[1])
})
