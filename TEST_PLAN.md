##SETUP
library(StratMonteCarlo)  # or devtools::load_all()
library(ggplot2)
set.seed(123)
n_samples <- 50000  # Use consistent sample size

##1D FUNCTIONS
# Simple polynomial functions
f1 <- function(x) x^2              # ∫[0,1] = 1/3 ≈ 0.333
f2 <- function(x) x^3              # ∫[0,1] = 1/4 = 0.25  
f3 <- function(x) 2*x              # ∫[0,1] = 1
f4 <- function(x) 1                # ∫[0,1] = 1

# Trigonometric functions
f5 <- function(x) sin(x)           # ∫[0,π] = 2
f6 <- function(x) cos(x)           # ∫[0,π/2] = 1
f7 <- function(x) sin(x)*cos(x)    # ∫[0,π/2] = 1/2

# Oscillatory (test partitioning)
f8 <- function(x) sin(10*x)        # ∫[0,2π] ≈ 0
f9 <- function(x) cos(5*x)         # ∫[0,2π] ≈ 0

# Exponential/rational
f10 <- function(x) exp(-x)         # ∫[0,∞] = 1, ∫[0,5] ≈ 0.993
f11 <- function(x) 1/(1+x^2)       # ∫[-1,1] = π/2 ≈ 1.571


##2D FUNCTIONS
# Simple 2D functions
f2d_1 <- function(x, y) x*y              # ∫[0,1]×[0,1] = 1/4 = 0.25
f2d_2 <- function(x, y) x^2 + y^2        # ∫[0,1]×[0,1] = 2/3 ≈ 0.667
f2d_3 <- function(x, y) 1                # ∫[0,1]×[0,1] = 1

# Gaussian-type (test partitioning effectiveness)
f2d_4 <- function(x, y) exp(-((x-0.5)^2 + (y-0.5)^2)/0.01)  # Sharp peak at (0.5,0.5)

## 3D FUNCTIONS
f3d_1 <- function(x, y, z) x*y*z         # ∫[0,1]³ = 1/8 = 0.125
f3d_2 <- function(x, y, z) 1             # ∫[0,1]³ = 1

## HOW TO RUN
#1D
mc_f1 <- montecarlo_integrate(f1, 0, 1, n_samples, FALSE)
mc_f1_part <- montecarlo_integrate(f1, 0, 1, n_samples, TRUE)
mc_print(mc_f1)
mc_print(mc_f1_part)
mc_plot(mc_f1 = mc_f1, mc_f1_part = mc_f1_part, true_value = 1/3)

#MULTI DIMENSIONAL
mc_f2d_1 <- montecarlo_integrate(f2d_1, c(0,0), c(pi*4, pi*4), n_samples, FALSE, dim = 2)
mc_f2d_1_part <- montecarlo_integrate(f2d_1, c(0,0), c(pi*4, pi*4), n_samples, TRUE, dim = 2)
mc_print(mc_f2d_1)
mc_print(mc_f2d_1_part)

## ERROR TESTING
montecarlo_integrate("not a function", 0, 1, 1000)
montecarlo_integrate(f1, "not numeric", 1, 1000)  
montecarlo_integrate(f1, 0, 1, -1000)  
