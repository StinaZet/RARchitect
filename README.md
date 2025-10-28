
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RARchitect

<!-- badges: start -->

<!-- badges: end -->

A framework for simulating Response-Adaptive Randomization (RAR)
clinical trial designs for different types of outcomes. Its purpose is
to assist researchers and statisticians in designing RAR trials by
allowing them to test various configurations and parameters, such as
prior settings, allocation tuning and clipping strategies, and
recruitment rates, through simulations.

## Installation

You can install the development version of RARchitect from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("StinaZet/RARchitect")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(RARchitect)
# Simulate a Binary Outcome Trial
# Example 1: Simulate a Binary Outcome Trial
set.seed(101)
results_binary <- simulate_brar_trial(
 outcome_type = "binary",
 distribution = "bernoulli",
 arms = 2, N = 100, blocksize = 10,
 modelpar = c(0.6, 0.4),
 priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
 tuning = 1, clipping = 0, burnin = 0, randmethod = "coin",
 postprobmethod = "exact", # Or "simulation"
 recruitment_rate = 5,
 observation_delay = 30)
head(results_binary)
#>   Batch Recruitment time Outcome time Arm Outcome AP arm 1 AP arm 2
#> 1     1                0           30   2       1      0.5      0.5
#> 2     1                0           30   2       1      0.5      0.5
#> 3     1                0           30   1       0      0.5      0.5
#> 4     1                0           30   1       0      0.5      0.5
#> 5     1                0           30   2       0      0.5      0.5
#> 6     1                0           30   2       0      0.5      0.5

# Example 2: Simulate a Normal Outcome Trial
set.seed(102)
model_params_norm <- matrix(c(10, 8, 2, 2), nrow = 2, byrow = TRUE)
prior_params_norm <- matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE)

results_normal <- simulate_brar_trial(
 outcome_type = "cont",
 distribution = "normal",
 arms = 2, N = 200, blocksize = 10, known_var = TRUE,
 priors = prior_params_norm,
 modelpar = model_params_norm,
 tuning = 1, clipping = 0, burnin = 10, randmethod = "urn",
 postprobmethod = "exact", # Or "simulation"
 recruitment_rate = 10,
 observation_delay = 15)
head(results_normal)
#>   Batch Recruitment time Outcome time Arm   Outcome AP arm 1 AP arm 2
#> 1     1                1           16   1  9.680793      0.5      0.5
#> 2     1                1           16   1 12.656393      0.5      0.5
#> 3     1                1           16   2  7.304414      0.5      0.5
#> 4     1                1           16   2 10.045077      0.5      0.5
#> 5     1                1           16   1 10.576071      0.5      0.5
#> 6     1                1           16   2  8.608396      0.5      0.5

# Example 3: Simulate an Exponential Outcome Trial
# True rates for arms: lambda1 = 0.2 (mean = 5), lambda2 = 0.1 (mean = 10)
set.seed(105)
prior_params_exp <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE) # Gamma(1,1) priors
true_rates_exp <- c(0.2, 0.1)

results_exponential <- simulate_brar_trial(
 outcome_type = "cont",
 distribution = "exponential",
 arms = 2, N = 120, blocksize = 1,
 priors = prior_params_exp,
 modelpar = true_rates_exp,
 tuning = 1, clipping = 0, burnin = 0,
 randmethod = "coin",
 postprobmethod = "simulation",
 recruitment_rate = 6,
 observation_delay = 5)
head(results_exponential)
#>   Batch Recruitment time Outcome time Arm     Outcome AP arm 1 AP arm 2
#> 1     1                0            5   2 25.94989670   0.5008   0.4992
#> 2     2                0            5   1  0.46815252   0.9294   0.0706
#> 3     3                0            5   1  0.49257058   0.9915   0.0085
#> 4     4                0            5   1  6.81570902   0.9983   0.0017
#> 5     5                0            5   1  0.02656646   0.9846   0.0154
#> 6     6                0            5   1  1.82558216   0.9959   0.0041
```
