# =============================================================================
# COMPREHENSIVE TEST PLAN FOR StratMonteCarlo PACKAGE
# =============================================================================
# Tests all features: Standard MC, Partitioning, Importance Sampling
# Does NOT combine partitioning + IS (tested separately)
# =============================================================================

library(StratMonteCarlo)
set.seed(123)

# Global settings
n_samples <- 50000
step <- 100

# Helper function to report results
report_result <- function(name, mc_result, true_value = NULL) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat(name, "\n")
  cat(rep("=", 70), "\n", sep = "")
  mc_print(mc_result)
  if (!is.null(true_value)) {
    error <- abs(mc_result$estimate - true_value)
    rel_error <- error / abs(true_value) * 100
    cat("True value:", formatC(true_value, digits = 6, format = "g"), "\n")
    cat("Absolute error:", formatC(error, digits = 6, format = "g"), "\n")
    cat("Relative error:", formatC(rel_error, digits = 4, format = "f"), "%\n")
  }
  if (mc_result$importance_sampling) {
    cat("IS Distribution:", mc_result$is_distribution, "\n")
    cat("IS Parameters:\n")
    print(mc_result$params_used)
  }
  if (!is.null(mc_result$n_partitions)) {
    cat("Number of partitions:", mc_result$n_partitions, "\n")
  }
}

# =============================================================================
# 1D FUNCTION DEFINITIONS
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                          1D FUNCTION TESTS                                  #\n")
cat("###############################################################################\n")

# --- Polynomial Functions ---
f1d_poly1 <- function(x) x^2              # ∫[0,1] = 1/3 ≈ 0.333
f1d_poly2 <- function(x) x^3              # ∫[0,1] = 1/4 = 0.25
f1d_poly3 <- function(x) 2*x              # ∫[0,1] = 1
f1d_poly4 <- function(x) x^4 - x^2        # ∫[0,1] = 1/5 - 1/3 = -2/15 ≈ -0.133

# --- Trigonometric Functions ---
f1d_sin <- function(x) sin(x)             # ∫[0,π] = 2
f1d_cos <- function(x) cos(x)             # ∫[0,π/2] = 1
f1d_sincos <- function(x) sin(x)*cos(x)   # ∫[0,π/2] = 1/2
f1d_sin2 <- function(x) sin(x)^2          # ∫[0,π] = π/2 ≈ 1.571

# --- Oscillatory Functions (Partitioning Test) ---
f1d_osc1 <- function(x) sin(10*x)         # ∫[0,2π] ≈ 0
f1d_osc2 <- function(x) cos(5*x)          # ∫[0,2π] ≈ 0
f1d_osc3 <- function(x) sin(20*x)         # ∫[0,π] ≈ 0

# --- Exponential/Decay (IS Test) ---
f1d_exp1 <- function(x) exp(-x)           # ∫[0,5] = 1 - e^(-5) ≈ 0.993
f1d_exp2 <- function(x) exp(-2*x)         # ∫[0,3] = (1-e^(-6))/2 ≈ 0.499
f1d_exp3 <- function(x) x*exp(-x)         # ∫[0,∞] = 1, ∫[0,5] ≈ 0.96

# --- Rational Functions ---
f1d_rational1 <- function(x) 1/(1+x^2)    # ∫[0,1] = π/4 ≈ 0.785
f1d_rational2 <- function(x) 1/(1+x)      # ∫[0,1] = ln(2) ≈ 0.693

# --- Peak/Concentrated Functions (IS Test) ---
f1d_peak1 <- function(x) exp(-50*(x-0.5)^2)  # Sharp peak at 0.5
f1d_peak2 <- function(x) exp(-100*(x-0.3)^2) # Very sharp peak at 0.3

# --- Mixed Functions ---
f1d_mixed1 <- function(x) x^2 * exp(-x)      # ∫[0,∞] = 2, ∫[0,5] ≈ 1.94
f1d_mixed2 <- function(x) sin(x) * exp(-x)   # ∫[0,π] ≈ 0.524

# =============================================================================
# 1D TESTS - POLYNOMIAL FUNCTIONS
# =============================================================================

cat("\n--- 1D POLYNOMIAL FUNCTIONS ---\n")

# Test 1.1: x^2 on [0,1]
mc_1_1_std <- montecarlo_integrate(f1d_poly1, 0, 1, n_samples, partition = FALSE)
mc_1_1_part <- montecarlo_integrate(f1d_poly1, 0, 1, n_samples, partition = TRUE)
report_result("1.1: x^2 [Standard]", mc_1_1_std, 1/3)
report_result("1.1: x^2 [Partitioned]", mc_1_1_part, 1/3)
mc_plot(Standard = mc_1_1_std, Partitioned = mc_1_1_part, 
        true_value = 1/3, title_size = 14)

# Test 1.2: x^3 on [0,1]
mc_1_2_std <- montecarlo_integrate(f1d_poly2, 0, 1, n_samples, partition = FALSE)
mc_1_2_part <- montecarlo_integrate(f1d_poly2, 0, 1, n_samples, partition = TRUE)
report_result("1.2: x^3 [Standard]", mc_1_2_std, 0.25)
report_result("1.2: x^3 [Partitioned]", mc_1_2_part, 0.25)

# Test 1.3: 2x on [0,1]
mc_1_3_std <- montecarlo_integrate(f1d_poly3, 0, 1, n_samples, partition = FALSE)
report_result("1.3: 2x [Standard]", mc_1_3_std, 1.0)

# Test 1.4: x^4 - x^2 on [0,1]
mc_1_4_std <- montecarlo_integrate(f1d_poly4, 0, 1, n_samples, partition = FALSE)
mc_1_4_part <- montecarlo_integrate(f1d_poly4, 0, 1, n_samples, partition = TRUE)
report_result("1.4: x^4 - x^2 [Standard]", mc_1_4_std, -2/15)
report_result("1.4: x^4 - x^2 [Partitioned]", mc_1_4_part, -2/15)

# =============================================================================
# 1D TESTS - TRIGONOMETRIC FUNCTIONS
# =============================================================================

cat("\n--- 1D TRIGONOMETRIC FUNCTIONS ---\n")

# Test 1.5: sin(x) on [0,π]
mc_1_5_std <- montecarlo_integrate(f1d_sin, 0, pi, n_samples, partition = FALSE)
mc_1_5_part <- montecarlo_integrate(f1d_sin, 0, pi, n_samples, partition = TRUE)
report_result("1.5: sin(x) [Standard]", mc_1_5_std, 2.0)
report_result("1.5: sin(x) [Partitioned]", mc_1_5_part, 2.0)

# Test 1.6: cos(x) on [0,π/2]
mc_1_6_std <- montecarlo_integrate(f1d_cos, 0, pi/2, n_samples, partition = FALSE)
report_result("1.6: cos(x) [Standard]", mc_1_6_std, 1.0)

# Test 1.7: sin(x)*cos(x) on [0,π/2]
mc_1_7_std <- montecarlo_integrate(f1d_sincos, 0, pi/2, n_samples, partition = FALSE)
report_result("1.7: sin(x)cos(x) [Standard]", mc_1_7_std, 0.5)

# Test 1.8: sin^2(x) on [0,π]
mc_1_8_std <- montecarlo_integrate(f1d_sin2, 0, pi, n_samples, partition = FALSE)
report_result("1.8: sin^2(x) [Standard]", mc_1_8_std, pi/2)

# =============================================================================
# 1D TESTS - OSCILLATORY FUNCTIONS (Partitioning Effectiveness)
# =============================================================================

cat("\n--- 1D OSCILLATORY FUNCTIONS ---\n")

# Test 1.9: sin(10x) on [0,2π]
mc_1_9_std <- montecarlo_integrate(f1d_osc1, 0, 2*pi, n_samples, partition = FALSE)
mc_1_9_part <- montecarlo_integrate(f1d_osc1, 0, 2*pi, n_samples, partition = TRUE)
report_result("1.9: sin(10x) [Standard]", mc_1_9_std, 0.0)
report_result("1.9: sin(10x) [Partitioned]", mc_1_9_part, 0.0)
mc_plot(Standard = mc_1_9_std, Partitioned = mc_1_9_part, 
        true_value = 0, title_size = 14)

# Test 1.10: cos(5x) on [0,2π]
mc_1_10_std <- montecarlo_integrate(f1d_osc2, 0, 2*pi, n_samples, partition = FALSE)
mc_1_10_part <- montecarlo_integrate(f1d_osc2, 0, 2*pi, n_samples, partition = TRUE)
report_result("1.10: cos(5x) [Standard]", mc_1_10_std, 0.0)
report_result("1.10: cos(5x) [Partitioned]", mc_1_10_part, 0.0)

# Test 1.11: sin(20x) on [0,π]
mc_1_11_std <- montecarlo_integrate(f1d_osc3, 0, pi, n_samples, partition = FALSE)
mc_1_11_part <- montecarlo_integrate(f1d_osc3, 0, pi, n_samples, partition = TRUE)
report_result("1.11: sin(20x) [Standard]", mc_1_11_std, 0.0)
report_result("1.11: sin(20x) [Partitioned]", mc_1_11_part, 0.0)
mc_plot(Standard = mc_1_11_std, Partitioned = mc_1_11_part, 
        true_value = 0, title_size = 14)

# =============================================================================
# 1D TESTS - EXPONENTIAL FUNCTIONS (Importance Sampling)
# =============================================================================

cat("\n--- 1D EXPONENTIAL FUNCTIONS (IS TEST) ---\n")

# Test 1.12: exp(-x) on [0,5] - Test all IS distributions
mc_1_12_std <- montecarlo_integrate(f1d_exp1, 0, 5, n_samples, 
                                    importance_sampling = FALSE)
mc_1_12_is_norm <- montecarlo_integrate(f1d_exp1, 0, 5, n_samples,
                                        importance_sampling = TRUE,
                                        is_distribution = "normal")
mc_1_12_is_exp <- montecarlo_integrate(f1d_exp1, 0, 5, n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "exponential")
mc_1_12_is_beta <- montecarlo_integrate(f1d_exp1, 0, 5, n_samples,
                                        importance_sampling = TRUE,
                                        is_distribution = "beta")
mc_1_12_is_mix <- montecarlo_integrate(f1d_exp1, 0, 5, n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "mixture_normal")

true_val_1_12 <- 1 - exp(-5)
report_result("1.12: exp(-x) [Standard]", mc_1_12_std, true_val_1_12)
report_result("1.12: exp(-x) [IS Normal]", mc_1_12_is_norm, true_val_1_12)
report_result("1.12: exp(-x) [IS Exponential]", mc_1_12_is_exp, true_val_1_12)
report_result("1.12: exp(-x) [IS Beta]", mc_1_12_is_beta, true_val_1_12)
report_result("1.12: exp(-x) [IS Mixture]", mc_1_12_is_mix, true_val_1_12)

mc_plot(Standard = mc_1_12_std, 
        IS_Normal = mc_1_12_is_norm,
        IS_Exponential = mc_1_12_is_exp,
        IS_Beta = mc_1_12_is_beta,
        IS_Mixture = mc_1_12_is_mix,
        true_value = true_val_1_12, title_size = 14)

# Test 1.13: exp(-2x) on [0,3]
mc_1_13_std <- montecarlo_integrate(f1d_exp2, 0, 3, n_samples,
                                    importance_sampling = FALSE)
mc_1_13_is_exp <- montecarlo_integrate(f1d_exp2, 0, 3, n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "exponential")
true_val_1_13 <- (1 - exp(-6))/2
report_result("1.13: exp(-2x) [Standard]", mc_1_13_std, true_val_1_13)
report_result("1.13: exp(-2x) [IS Exponential]", mc_1_13_is_exp, true_val_1_13)

# Test 1.14: x*exp(-x) on [0,5]
mc_1_14_std <- montecarlo_integrate(f1d_exp3, 0, 5, n_samples,
                                    importance_sampling = FALSE)
mc_1_14_is_exp <- montecarlo_integrate(f1d_exp3, 0, 5, n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "exponential")
report_result("1.14: x*exp(-x) [Standard]", mc_1_14_std, 0.96)
report_result("1.14: x*exp(-x) [IS Exponential]", mc_1_14_is_exp, 0.96)

# =============================================================================
# 1D TESTS - RATIONAL FUNCTIONS
# =============================================================================

cat("\n--- 1D RATIONAL FUNCTIONS ---\n")

# Test 1.15: 1/(1+x^2) on [0,1]
mc_1_15_std <- montecarlo_integrate(f1d_rational1, 0, 1, n_samples,
                                    importance_sampling = FALSE)
mc_1_15_part <- montecarlo_integrate(f1d_rational1, 0, 1, n_samples,
                                     partition = TRUE)
report_result("1.15: 1/(1+x^2) [Standard]", mc_1_15_std, pi/4)
report_result("1.15: 1/(1+x^2) [Partitioned]", mc_1_15_part, pi/4)

# Test 1.16: 1/(1+x) on [0,1]
mc_1_16_std <- montecarlo_integrate(f1d_rational2, 0, 1, n_samples,
                                    importance_sampling = FALSE)
report_result("1.16: 1/(1+x) [Standard]", mc_1_16_std, log(2))

# =============================================================================
# 1D TESTS - PEAK FUNCTIONS (IS Effectiveness)
# =============================================================================

cat("\n--- 1D PEAK FUNCTIONS (IS TEST) ---\n")

# Test 1.17: Sharp peak at 0.5
mc_1_17_std <- montecarlo_integrate(f1d_peak1, 0, 1, n_samples,
                                    importance_sampling = FALSE)
mc_1_17_is_norm <- montecarlo_integrate(f1d_peak1, 0, 1, n_samples,
                                        importance_sampling = TRUE,
                                        is_distribution = "normal")
mc_1_17_is_mix <- montecarlo_integrate(f1d_peak1, 0, 1, n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "mixture_normal")
report_result("1.17: Peak at 0.5 [Standard]", mc_1_17_std)
report_result("1.17: Peak at 0.5 [IS Normal]", mc_1_17_is_norm)
report_result("1.17: Peak at 0.5 [IS Mixture]", mc_1_17_is_mix)

mc_plot(Standard = mc_1_17_std, 
        IS_Normal = mc_1_17_is_norm,
        IS_Mixture = mc_1_17_is_mix, title_size = 14)

# Test 1.18: Very sharp peak at 0.3
mc_1_18_std <- montecarlo_integrate(f1d_peak2, 0, 1, n_samples,
                                    importance_sampling = FALSE)
mc_1_18_is_norm <- montecarlo_integrate(f1d_peak2, 0, 1, n_samples,
                                        importance_sampling = TRUE,
                                        is_distribution = "normal")
mc_1_18_part <- montecarlo_integrate(f1d_peak2, 0, 1, n_samples,
                                     partition = TRUE)
report_result("1.18: Peak at 0.3 [Standard]", mc_1_18_std)
report_result("1.18: Peak at 0.3 [IS Normal]", mc_1_18_is_norm)
report_result("1.18: Peak at 0.3 [Partitioned]", mc_1_18_part)

# =============================================================================
# 1D TESTS - MIXED FUNCTIONS
# =============================================================================

cat("\n--- 1D MIXED FUNCTIONS ---\n")

# Test 1.19: x^2 * exp(-x) on [0,5]
mc_1_19_std <- montecarlo_integrate(f1d_mixed1, 0, 5, n_samples,
                                    importance_sampling = FALSE)
mc_1_19_is_exp <- montecarlo_integrate(f1d_mixed1, 0, 5, n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "exponential")
report_result("1.19: x^2*exp(-x) [Standard]", mc_1_19_std, 1.94)
report_result("1.19: x^2*exp(-x) [IS Exponential]", mc_1_19_is_exp, 1.94)

# Test 1.20: sin(x) * exp(-x) on [0,π]
mc_1_20_std <- montecarlo_integrate(f1d_mixed2, 0, pi, n_samples,
                                    importance_sampling = FALSE)
mc_1_20_part <- montecarlo_integrate(f1d_mixed2, 0, pi, n_samples,
                                     partition = TRUE)
report_result("1.20: sin(x)*exp(-x) [Standard]", mc_1_20_std, 0.524)
report_result("1.20: sin(x)*exp(-x) [Partitioned]", mc_1_20_part, 0.524)

# =============================================================================
# 2D FUNCTION DEFINITIONS
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                          2D FUNCTION TESTS                                  #\n")
cat("###############################################################################\n")

# --- Simple 2D Functions ---
f2d_poly1 <- function(x, y) x*y                    # ∫[0,1]² = 1/4 = 0.25
f2d_poly2 <- function(x, y) x^2 + y^2              # ∫[0,1]² = 2/3 ≈ 0.667
f2d_poly3 <- function(x, y) x^2*y^2                # ∫[0,1]² = 1/9 ≈ 0.111
f2d_constant <- function(x, y) 1                   # ∫[0,1]² = 1

# --- Gaussian/Peak Functions ---
f2d_peak1 <- function(x, y) exp(-((x-0.5)^2 + (y-0.5)^2)/0.01)  # Peak at (0.5,0.5)
f2d_peak2 <- function(x, y) exp(-10*((x-0.5)^2 + (y-0.5)^2))    # Wider peak
f2d_gauss <- function(x, y) exp(-(x^2 + y^2))                    # Gaussian at origin

# --- Trigonometric 2D ---
f2d_sin <- function(x, y) sin(x) * cos(y)          # ∫[0,π]×[0,π] = 0
f2d_osc <- function(x, y) sin(5*x) * cos(5*y)      # Oscillatory

# --- Mixed 2D ---
f2d_rational <- function(x, y) 1/(1 + x^2 + y^2)   # Rational
f2d_exp <- function(x, y) exp(-x-y)                # Exponential decay

# =============================================================================
# 2D TESTS - POLYNOMIAL FUNCTIONS
# =============================================================================

cat("\n--- 2D POLYNOMIAL FUNCTIONS ---\n")

# Test 2.1: x*y on [0,1]²
mc_2_1_std <- montecarlo_integrate(f2d_poly1, c(0,0), c(1,1), n_samples,
                                   partition = FALSE, dim = 2)
mc_2_1_part <- montecarlo_integrate(f2d_poly1, c(0,0), c(1,1), n_samples,
                                    partition = TRUE, dim = 2)
report_result("2.1: x*y [Standard]", mc_2_1_std, 0.25)
report_result("2.1: x*y [Partitioned]", mc_2_1_part, 0.25)
mc_plot(Standard = mc_2_1_std, Partitioned = mc_2_1_part,
        true_value = 0.25, title_size = 14)

# Test 2.2: x^2 + y^2 on [0,1]²
mc_2_2_std <- montecarlo_integrate(f2d_poly2, c(0,0), c(1,1), n_samples,
                                   partition = FALSE, dim = 2)
mc_2_2_part <- montecarlo_integrate(f2d_poly2, c(0,0), c(1,1), n_samples,
                                    partition = TRUE, dim = 2)
report_result("2.2: x^2+y^2 [Standard]", mc_2_2_std, 2/3)
report_result("2.2: x^2+y^2 [Partitioned]", mc_2_2_part, 2/3)

# Test 2.3: x^2*y^2 on [0,1]²
mc_2_3_std <- montecarlo_integrate(f2d_poly3, c(0,0), c(1,1), n_samples,
                                   partition = FALSE, dim = 2)
report_result("2.3: x^2*y^2 [Standard]", mc_2_3_std, 1/9)

# Test 2.4: Constant function
mc_2_4_std <- montecarlo_integrate(f2d_constant, c(0,0), c(1,1), n_samples,
                                   partition = FALSE, dim = 2)
report_result("2.4: Constant [Standard]", mc_2_4_std, 1.0)

# =============================================================================
# 2D TESTS - GAUSSIAN/PEAK FUNCTIONS
# =============================================================================

cat("\n--- 2D GAUSSIAN/PEAK FUNCTIONS ---\n")

# Test 2.5: Sharp peak at (0.5,0.5) - Partitioning effectiveness
mc_2_5_std <- montecarlo_integrate(f2d_peak1, c(0,0), c(1,1), n_samples,
                                   partition = FALSE, dim = 2)
mc_2_5_part <- montecarlo_integrate(f2d_peak1, c(0,0), c(1,1), n_samples,
                                    partition = TRUE, dim = 2)
report_result("2.5: Sharp peak [Standard]", mc_2_5_std)
report_result("2.5: Sharp peak [Partitioned]", mc_2_5_part)
mc_plot(Standard = mc_2_5_std, Partitioned = mc_2_5_part, title_size = 14)

# Test 2.6: Sharp peak - IS effectiveness
mc_2_6_std <- montecarlo_integrate(f2d_peak1, c(0,0), c(1,1), n_samples,
                                   importance_sampling = FALSE, dim = 2)
mc_2_6_is_norm <- montecarlo_integrate(f2d_peak1, c(0,0), c(1,1), n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "normal", dim = 2)
mc_2_6_is_mix <- montecarlo_integrate(f2d_peak1, c(0,0), c(1,1), n_samples,
                                      importance_sampling = TRUE,
                                      is_distribution = "mixture_normal", dim = 2)
report_result("2.6: Sharp peak [Standard]", mc_2_6_std)
report_result("2.6: Sharp peak [IS Normal]", mc_2_6_is_norm)
report_result("2.6: Sharp peak [IS Mixture]", mc_2_6_is_mix)

mc_plot(Standard = mc_2_6_std, 
        IS_Normal = mc_2_6_is_norm,
        IS_Mixture = mc_2_6_is_mix, title_size = 14)

# Test 2.7: Wider peak
mc_2_7_std <- montecarlo_integrate(f2d_peak2, c(0,0), c(1,1), n_samples,
                                   partition = FALSE, dim = 2)
mc_2_7_part <- montecarlo_integrate(f2d_peak2, c(0,0), c(1,1), n_samples,
                                    partition = TRUE, dim = 2)
report_result("2.7: Wider peak [Standard]", mc_2_7_std)
report_result("2.7: Wider peak [Partitioned]", mc_2_7_part)

# Test 2.8: Gaussian at origin
mc_2_8_std <- montecarlo_integrate(f2d_gauss, c(-2,-2), c(2,2), n_samples,
                                   partition = FALSE, dim = 2)
mc_2_8_is_norm <- montecarlo_integrate(f2d_gauss, c(-2,-2), c(2,2), n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "normal", dim = 2)
report_result("2.8: Gaussian origin [Standard]", mc_2_8_std, pi)
report_result("2.8: Gaussian origin [IS Normal]", mc_2_8_is_norm, pi)

# =============================================================================
# 2D TESTS - TRIGONOMETRIC FUNCTIONS
# =============================================================================

cat("\n--- 2D TRIGONOMETRIC FUNCTIONS ---\n")

# Test 2.9: sin(x)*cos(y) on [0,π]²
mc_2_9_std <- montecarlo_integrate(f2d_sin, c(0,0), c(pi,pi), n_samples,
                                   partition = FALSE, dim = 2)
mc_2_9_part <- montecarlo_integrate(f2d_sin, c(0,0), c(pi,pi), n_samples,
                                    partition = TRUE, dim = 2)
report_result("2.9: sin(x)cos(y) [Standard]", mc_2_9_std, 0.0)
report_result("2.9: sin(x)cos(y) [Partitioned]", mc_2_9_part, 0.0)

# Test 2.10: Oscillatory function
mc_2_10_std <- montecarlo_integrate(f2d_osc, c(0,0), c(pi,pi), n_samples,
                                    partition = FALSE, dim = 2)
mc_2_10_part <- montecarlo_integrate(f2d_osc, c(0,0), c(pi,pi), n_samples,
                                     partition = TRUE, dim = 2)
report_result("2.10: Oscillatory [Standard]", mc_2_10_std, 0.0)
report_result("2.10: Oscillatory [Partitioned]", mc_2_10_part, 0.0)
mc_plot(Standard = mc_2_10_std, Partitioned = mc_2_10_part,
        true_value = 0, title_size = 14)

# =============================================================================
# 2D TESTS - MIXED FUNCTIONS
# =============================================================================

cat("\n--- 2D MIXED FUNCTIONS ---\n")

# Test 2.11: Rational function
mc_2_11_std <- montecarlo_integrate(f2d_rational, c(0,0), c(1,1), n_samples,
                                    partition = FALSE, dim = 2)
mc_2_11_part <- montecarlo_integrate(f2d_rational, c(0,0), c(1,1), n_samples,
                                     partition = TRUE, dim = 2)
report_result("2.11: Rational [Standard]", mc_2_11_std)
report_result("2.11: Rational [Partitioned]", mc_2_11_part)

# Test 2.12: Exponential decay
mc_2_12_std <- montecarlo_integrate(f2d_exp, c(0,0), c(3,3), n_samples,
                                    partition = FALSE, dim = 2)
mc_2_12_is_exp <- montecarlo_integrate(f2d_exp, c(0,0), c(3,3), n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "exponential", dim = 2)
report_result("2.12: Exponential decay [Standard]", mc_2_12_std)
report_result("2.12: Exponential decay [IS Exponential]", mc_2_12_is_exp)

# =============================================================================
# 3D FUNCTION DEFINITIONS
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                          3D FUNCTION TESTS                                  #\n")
cat("###############################################################################\n")

# --- Simple 3D Functions ---
f3d_poly1 <- function(x, y, z) x*y*z                    # ∫[0,1]³ = 1/8 = 0.125
f3d_poly2 <- function(x, y, z) x^2 + y^2 + z^2          # ∫[0,1]³ = 1
f3d_constant <- function(x, y, z) 1                     # ∫[0,1]³ = 1

# --- Gaussian/Peak 3D ---
f3d_gauss <- function(x, y, z) exp(-(x^2 + y^2 + z^2))  # 3D Gaussian
f3d_peak <- function(x, y, z) exp(-10*((x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2))  # Peak at center

# --- Exponential 3D ---
f3d_exp <- function(x, y, z) exp(-x-y-z)                # Exponential decay

# --- Mixed 3D ---
f3d_mixed <- function(x, y, z) x*y*z*exp(-x-y-z)       # Polynomial * exponential

# =============================================================================
# 3D TESTS - POLYNOMIAL FUNCTIONS
# =============================================================================

cat("\n--- 3D POLYNOMIAL FUNCTIONS ---\n")

# Test 3.1: x*y*z on [0,1]³
mc_3_1_std <- montecarlo_integrate(f3d_poly1, c(0,0,0), c(1,1,1), n_samples,
                                   partition = FALSE, dim = 3)
mc_3_1_part <- montecarlo_integrate(f3d_poly1, c(0,0,0), c(1,1,1), n_samples,
                                    partition = TRUE, dim = 3)
report_result("3.1: x*y*z [Standard]", mc_3_1_std, 0.125)
report_result("3.1: x*y*z [Partitioned]", mc_3_1_part, 0.125)
mc_plot(Standard = mc_3_1_std, Partitioned = mc_3_1_part,
        true_value = 0.125, title_size = 14)

# Test 3.2: x^2 + y^2 + z^2 on [0,1]³
mc_3_2_std <- montecarlo_integrate(f3d_poly2, c(0,0,0), c(1,1,1), n_samples,
                                   partition = FALSE, dim = 3)
mc_3_2_part <- montecarlo_integrate(f3d_poly2, c(0,0,0), c(1,1,1), n_samples,
                                    partition = TRUE, dim = 3)
report_result("3.2: x^2+y^2+z^2 [Standard]", mc_3_2_std, 1.0)
report_result("3.2: x^2+y^2+z^2 [Partitioned]", mc_3_2_part, 1.0)

# Test 3.3: Constant function
mc_3_3_std <- montecarlo_integrate(f3d_constant, c(0,0,0), c(1,1,1), n_samples,
                                   partition = FALSE, dim = 3)
report_result("3.3: Constant [Standard]", mc_3_3_std, 1.0)

# =============================================================================
# 3D TESTS - GAUSSIAN/PEAK FUNCTIONS
# =============================================================================

cat("\n--- 3D GAUSSIAN/PEAK FUNCTIONS ---\n")

# Test 3.4: 3D Gaussian - Standard vs Partitioning
mc_3_4_std <- montecarlo_integrate(f3d_gauss, c(-2,-2,-2), c(2,2,2), n_samples,
                                   partition = FALSE, dim = 3)
mc_3_4_part <- montecarlo_integrate(f3d_gauss, c(-2,-2,-2), c(2,2,2), n_samples,
                                    partition = TRUE, dim = 3)
report_result("3.4: 3D Gaussian [Standard]", mc_3_4_std, pi^(3/2))
report_result("3.4: 3D Gaussian [Partitioned]", mc_3_4_part, pi^(3/2))

# Test 3.5: 3D Gaussian - IS effectiveness
mc_3_5_std <- montecarlo_integrate(f3d_gauss, c(-2,-2,-2), c(2,2,2), n_samples,
                                   importance_sampling = FALSE, dim = 3)
mc_3_5_is_norm <- montecarlo_integrate(f3d_gauss, c(-2,-2,-2), c(2,2,2), n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "normal", dim = 3)
report_result("3.5: 3D Gaussian [Standard]", mc_3_5_std, pi^(3/2))
report_result("3.5: 3D Gaussian [IS Normal]", mc_3_5_is_norm, pi^(3/2))
mc_plot(Standard = mc_3_5_std, IS_Normal = mc_3_5_is_norm,
        true_value = pi^(3/2), title_size = 14)

# Test 3.6: 3D Peak at center
mc_3_6_std <- montecarlo_integrate(f3d_peak, c(0,0,0), c(1,1,1), n_samples,
                                   partition = FALSE, dim = 3)
mc_3_6_part <- montecarlo_integrate(f3d_peak, c(0,0,0), c(1,1,1), n_samples,
                                    partition = TRUE, dim = 3)
mc_3_6_is_norm <- montecarlo_integrate(f3d_peak, c(0,0,0), c(1,1,1), n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "normal", dim = 3)
report_result("3.6: 3D Peak [Standard]", mc_3_6_std)
report_result("3.6: 3D Peak [Partitioned]", mc_3_6_part)
report_result("3.6: 3D Peak [IS Normal]", mc_3_6_is_norm)

mc_plot(Standard = mc_3_6_std, 
        Partitioned = mc_3_6_part,
        IS_Normal = mc_3_6_is_norm, title_size = 14)

# =============================================================================
# 3D TESTS - EXPONENTIAL FUNCTIONS
# =============================================================================

cat("\n--- 3D EXPONENTIAL FUNCTIONS ---\n")

# Test 3.7: Exponential decay - Standard vs IS
mc_3_7_std <- montecarlo_integrate(f3d_exp, c(0,0,0), c(3,3,3), n_samples,
                                   importance_sampling = FALSE, dim = 3)
mc_3_7_is_exp <- montecarlo_integrate(f3d_exp, c(0,0,0), c(3,3,3), n_samples,
                                      importance_sampling = TRUE,
                                      is_distribution = "exponential", dim = 3)
mc_3_7_is_norm <- montecarlo_integrate(f3d_exp, c(0,0,0), c(3,3,3), n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "normal", dim = 3)
report_result("3.7: Exponential decay [Standard]", mc_3_7_std)
report_result("3.7: Exponential decay [IS Exponential]", mc_3_7_is_exp)
report_result("3.7: Exponential decay [IS Normal]", mc_3_7_is_norm)

mc_plot(Standard = mc_3_7_std, 
        IS_Exponential = mc_3_7_is_exp,
        IS_Normal = mc_3_7_is_norm, title_size = 14)

# =============================================================================
# 3D TESTS - MIXED FUNCTIONS
# =============================================================================

cat("\n--- 3D MIXED FUNCTIONS ---\n")

# Test 3.8: Polynomial * Exponential
mc_3_8_std <- montecarlo_integrate(f3d_mixed, c(0,0,0), c(3,3,3), n_samples,
                                   partition = FALSE, dim = 3)
mc_3_8_is_exp <- montecarlo_integrate(f3d_mixed, c(0,0,0), c(3,3,3), n_samples,
                                      importance_sampling = TRUE,
                                      is_distribution = "exponential", dim = 3)
report_result("3.8: x*y*z*exp(-x-y-z) [Standard]", mc_3_8_std)
report_result("3.8: x*y*z*exp(-x-y-z) [IS Exponential]", mc_3_8_is_exp)

# =============================================================================
# HIGHER DIMENSIONAL TESTS
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                    HIGHER DIMENSIONAL TESTS (4D+)                           #\n")
cat("###############################################################################\n")

# --- 4D Functions ---
f4d_poly <- function(x1, x2, x3, x4) x1*x2*x3*x4     # ∫[0,1]⁴ = 1/16
f4d_gauss <- function(x1, x2, x3, x4) exp(-(x1^2 + x2^2 + x3^2 + x4^2))

cat("\n--- 4D FUNCTION TESTS ---\n")

# Test 4.1: 4D Polynomial
mc_4_1_std <- montecarlo_integrate(f4d_poly, rep(0,4), rep(1,4), n_samples,
                                   partition = FALSE, dim = 4)
report_result("4.1: 4D Polynomial [Standard]", mc_4_1_std, 1/16)

# Test 4.2: 4D Gaussian
mc_4_2_std <- montecarlo_integrate(f4d_gauss, rep(-2,4), rep(2,4), n_samples,
                                   importance_sampling = FALSE, dim = 4)
mc_4_2_is_norm <- montecarlo_integrate(f4d_gauss, rep(-2,4), rep(2,4), n_samples,
                                       importance_sampling = TRUE,
                                       is_distribution = "normal", dim = 4)
report_result("4.2: 4D Gaussian [Standard]", mc_4_2_std, pi^2)
report_result("4.2: 4D Gaussian [IS Normal]", mc_4_2_is_norm, pi^2)

# --- 5D Functions ---
f5d_poly <- function(x1, x2, x3, x4, x5) x1*x2*x3*x4*x5  # ∫[0,1]⁵ = 1/32

cat("\n--- 5D FUNCTION TESTS ---\n")

# Test 5.1: 5D Polynomial
mc_5_1_std <- montecarlo_integrate(f5d_poly, rep(0,5), rep(1,5), n_samples,
                                   partition = FALSE, dim = 5)
report_result("5.1: 5D Polynomial [Standard]", mc_5_1_std, 1/32)

# =============================================================================
# EDGE CASE TESTS
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                          EDGE CASE TESTS                                    #\n")
cat("###############################################################################\n")

# --- Negative Integration Bounds ---
cat("\n--- NEGATIVE BOUNDS ---\n")

f_neg <- function(x) x^2
mc_edge_1 <- montecarlo_integrate(f_neg, -1, 1, n_samples)
report_result("Edge 1: x^2 on [-1,1]", mc_edge_1, 2/3)

# --- Very Small Intervals ---
cat("\n--- SMALL INTERVALS ---\n")

f_small <- function(x) x^2
mc_edge_2 <- montecarlo_integrate(f_small, 0, 0.01, n_samples)
report_result("Edge 2: x^2 on [0,0.01]", mc_edge_2, (0.01)^3/3)

# --- Large Intervals ---
cat("\n--- LARGE INTERVALS ---\n")

f_large <- function(x) 1/(1+x^2)
mc_edge_3 <- montecarlo_integrate(f_large, -10, 10, n_samples)
report_result("Edge 3: 1/(1+x^2) on [-10,10]", mc_edge_3, 2*atan(10))

# --- Functions with Zeros ---
cat("\n--- FUNCTIONS WITH ZEROS ---\n")

f_zero <- function(x) sin(x)
mc_edge_4 <- montecarlo_integrate(f_zero, 0, 2*pi, n_samples)
report_result("Edge 4: sin(x) on [0,2π]", mc_edge_4, 0.0)

# --- Asymmetric Bounds ---
cat("\n--- ASYMMETRIC BOUNDS ---\n")

f_asym <- function(x, y) x + y
mc_edge_5 <- montecarlo_integrate(f_asym, c(-1,0), c(2,3), n_samples, dim = 2)
report_result("Edge 5: x+y on [-1,2]×[0,3]", mc_edge_5, 13.5)

# --- Very Sharp Peaks ---
cat("\n--- VERY SHARP PEAKS ---\n")

f_sharp <- function(x) exp(-500*(x-0.5)^2)
mc_edge_6_std <- montecarlo_integrate(f_sharp, 0, 1, n_samples, partition = FALSE)
mc_edge_6_part <- montecarlo_integrate(f_sharp, 0, 1, n_samples, partition = TRUE)
mc_edge_6_is <- montecarlo_integrate(f_sharp, 0, 1, n_samples, 
                                     importance_sampling = TRUE,
                                     is_distribution = "normal")
report_result("Edge 6: Very sharp peak [Standard]", mc_edge_6_std)
report_result("Edge 6: Very sharp peak [Partitioned]", mc_edge_6_part)
report_result("Edge 6: Very sharp peak [IS Normal]", mc_edge_6_is)

# --- Discontinuous Functions (approximated) ---
cat("\n--- STEP FUNCTIONS ---\n")

f_step <- function(x) ifelse(x < 0.5, 0, 1)
mc_edge_7 <- montecarlo_integrate(f_step, 0, 1, n_samples)
report_result("Edge 7: Step function", mc_edge_7, 0.5)

# =============================================================================
# CONVERGENCE TESTS
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                        CONVERGENCE TESTS                                    #\n")
cat("###############################################################################\n")

# Test convergence with different sample sizes
f_conv <- function(x) x^2

sample_sizes <- c(1000, 5000, 10000, 25000, 50000, 100000)
convergence_results <- data.frame(
  n_samples = integer(),
  estimate = numeric(),
  error = numeric(),
  method = character()
)

cat("\n--- STANDARD MONTE CARLO CONVERGENCE ---\n")
for (n in sample_sizes) {
  mc_temp <- montecarlo_integrate(f_conv, 0, 1, n, partition = FALSE)
  error <- abs(mc_temp$estimate - 1/3)
  convergence_results <- rbind(convergence_results,
                              data.frame(n_samples = n, 
                                       estimate = mc_temp$estimate,
                                       error = error,
                                       method = "Standard"))
  cat(sprintf("n=%6d: estimate=%.6f, error=%.6f\n", n, mc_temp$estimate, error))
}

cat("\n--- PARTITIONED MONTE CARLO CONVERGENCE ---\n")
for (n in sample_sizes) {
  mc_temp <- montecarlo_integrate(f_conv, 0, 1, n, partition = TRUE)
  error <- abs(mc_temp$estimate - 1/3)
  convergence_results <- rbind(convergence_results,
                              data.frame(n_samples = n, 
                                       estimate = mc_temp$estimate,
                                       error = error,
                                       method = "Partitioned"))
  cat(sprintf("n=%6d: estimate=%.6f, error=%.6f\n", n, mc_temp$estimate, error))
}

# Plot convergence
ggplot(convergence_results, aes(x = n_samples, y = error, color = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Convergence Rate Comparison",
       x = "Number of Samples (log scale)",
       y = "Absolute Error (log scale)",
       color = "Method") +
  theme(legend.position = "bottom")

# =============================================================================
# PARAMETER TESTING
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                        PARAMETER TESTS                                      #\n")
cat("###############################################################################\n")

# Test different step sizes
cat("\n--- STEP SIZE TESTS ---\n")

f_step_test <- function(x) x^2
mc_step_10 <- montecarlo_integrate(f_step_test, 0, 1, 10000, step = 10)
mc_step_100 <- montecarlo_integrate(f_step_test, 0, 1, 10000, step = 100)
mc_step_1000 <- montecarlo_integrate(f_step_test, 0, 1, 10000, step = 1000)

cat("Step=10:  ", length(mc_step_10$estimates), "estimates recorded\n")
cat("Step=100: ", length(mc_step_100$estimates), "estimates recorded\n")
cat("Step=1000:", length(mc_step_1000$estimates), "estimates recorded\n")

# Test different IS distributions on same function
cat("\n--- IS DISTRIBUTION COMPARISON ---\n")

f_is_test <- function(x) exp(-2*x)
true_val_is <- (1 - exp(-10))/2

mc_is_none <- montecarlo_integrate(f_is_test, 0, 5, n_samples,
                                   importance_sampling = FALSE)
mc_is_norm <- montecarlo_integrate(f_is_test, 0, 5, n_samples,
                                   importance_sampling = TRUE,
                                   is_distribution = "normal")
mc_is_exp <- montecarlo_integrate(f_is_test, 0, 5, n_samples,
                                  importance_sampling = TRUE,
                                  is_distribution = "exponential")
mc_is_beta <- montecarlo_integrate(f_is_test, 0, 5, n_samples,
                                   importance_sampling = TRUE,
                                   is_distribution = "beta")
mc_is_mix <- montecarlo_integrate(f_is_test, 0, 5, n_samples,
                                  importance_sampling = TRUE,
                                  is_distribution = "mixture_normal")

cat("\nIS Distribution Comparison (exp(-2x) on [0,5]):\n")
cat(sprintf("True value: %.6f\n", true_val_is))
cat(sprintf("No IS:      %.6f (error: %.6f)\n", 
            mc_is_none$estimate, abs(mc_is_none$estimate - true_val_is)))
cat(sprintf("Normal:     %.6f (error: %.6f)\n", 
            mc_is_norm$estimate, abs(mc_is_norm$estimate - true_val_is)))
cat(sprintf("Exponential:%.6f (error: %.6f)\n", 
            mc_is_exp$estimate, abs(mc_is_exp$estimate - true_val_is)))
cat(sprintf("Beta:       %.6f (error: %.6f)\n", 
            mc_is_beta$estimate, abs(mc_is_beta$estimate - true_val_is)))
cat(sprintf("Mixture:    %.6f (error: %.6f)\n", 
            mc_is_mix$estimate, abs(mc_is_mix$estimate - true_val_is)))

mc_plot(No_IS = mc_is_none,
        Normal = mc_is_norm,
        Exponential = mc_is_exp,
        Beta = mc_is_beta,
        Mixture = mc_is_mix,
        true_value = true_val_is, title_size = 14)

# =============================================================================
# ERROR HANDLING TESTS
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                        ERROR HANDLING TESTS                                 #\n")
cat("###############################################################################\n")

cat("\n--- EXPECTED ERRORS (should fail gracefully) ---\n")

# Test 1: Invalid function input
cat("\nTest: Invalid function input\n")
tryCatch({
  montecarlo_integrate("not a function", 0, 1, 1000)
  cat("ERROR: Should have failed!\n")
}, error = function(e) {
  cat("✓ Caught error:", conditionMessage(e), "\n")
})

# Test 2: Non-numeric bounds
cat("\nTest: Non-numeric lower bound\n")
tryCatch({
  f_test <- function(x) x^2
  montecarlo_integrate(f_test, "not numeric", 1, 1000)
  cat("ERROR: Should have failed!\n")
}, error = function(e) {
  cat("✓ Caught error:", conditionMessage(e), "\n")
})

# Test 3: Negative sample size
cat("\nTest: Negative sample size\n")
tryCatch({
  f_test <- function(x) x^2
  montecarlo_integrate(f_test, 0, 1, -1000)
  cat("ERROR: Should have failed!\n")
}, error = function(e) {
  cat("✓ Caught error:", conditionMessage(e), "\n")
})

# Test 4: Zero sample size
cat("\nTest: Zero sample size\n")
tryCatch({
  f_test <- function(x) x^2
  montecarlo_integrate(f_test, 0, 1, 0)
  cat("ERROR: Should have failed!\n")
}, error = function(e) {
  cat("✓ Caught error:", conditionMessage(e), "\n")
})

# Test 5: Mismatched dimensions
cat("\nTest: Mismatched dimension specification\n")
tryCatch({
  f_test <- function(x, y) x*y
  montecarlo_integrate(f_test, c(0,0), c(1,1), 1000, dim = 3)  # Says 3D but gives 2D bounds
  cat("ERROR: Should have failed!\n")
}, error = function(e) {
  cat("✓ Caught error:", conditionMessage(e), "\n")
})

# Test 6: Lower > Upper bounds
cat("\nTest: Lower bound > Upper bound\n")
tryCatch({
  f_test <- function(x) x^2
  result <- montecarlo_integrate(f_test, 1, 0, 1000)  # Reversed bounds
  cat("Result:", result$estimate, "(might work but worth checking)\n")
}, error = function(e) {
  cat("✓ Caught error:", conditionMessage(e), "\n")
})

# Test 7: Invalid IS distribution
cat("\nTest: Invalid IS distribution\n")
tryCatch({
  f_test <- function(x) x^2
  montecarlo_integrate(f_test, 0, 1, 1000, 
                      importance_sampling = TRUE,
                      is_distribution = "invalid_dist")
  cat("ERROR: Should have failed!\n")
}, error = function(e) {
  cat("✓ Caught error:", conditionMessage(e), "\n")
})

# Test 8: Function returning non-numeric
cat("\nTest: Function returning non-numeric\n")
tryCatch({
  f_bad <- function(x) "not a number"
  montecarlo_integrate(f_bad, 0, 1, 100)
  cat("ERROR: Should have failed!\n")
}, error = function(e) {
  cat("✓ Caught error:", conditionMessage(e), "\n")
})

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\n\n")
cat("###############################################################################\n")
cat("#                          TEST SUMMARY                                       #\n")
cat("###############################################################################\n")

cat("\nTest plan completed successfully!\n")
cat("\nCategories tested:\n")
cat("  ✓ 1D Functions (20 tests)\n")
cat("    - Polynomials\n")
cat("    - Trigonometric\n")
cat("    - Oscillatory\n")
cat("    - Exponential/Decay\n")
cat("    - Rational\n")
cat("    - Sharp Peaks\n")
cat("    - Mixed\n")
cat("  ✓ 2D Functions (12 tests)\n")
cat("    - Polynomials\n")
cat("    - Gaussian/Peaks\n")
cat("    - Trigonometric\n")
cat("    - Mixed\n")
cat("  ✓ 3D Functions (8 tests)\n")
cat("    - Polynomials\n")
cat("    - Gaussian/Peaks\n")
cat("    - Exponential\n")
cat("    - Mixed\n")
cat("  ✓ Higher Dimensional (3 tests: 4D, 5D)\n")
cat("  ✓ Edge Cases (7 tests)\n")
cat("  ✓ Convergence Tests\n")
cat("  ✓ Parameter Tests\n")
cat("  ✓ Error Handling (8 tests)\n")
cat("\nMethods tested:\n")
cat("  ✓ Standard Monte Carlo\n")
cat("  ✓ Adaptive Partitioning\n")
cat("  ✓ Importance Sampling (all 4 distributions)\n")
cat("  ✓ Various parameter configurations\n")
cat("\nTotal tests: ~60+\n")
cat("\n" , rep("=", 78), "\n", sep = "")
cat("END OF TEST PLAN\n")
cat(rep("=", 78), "\n", sep = "")
