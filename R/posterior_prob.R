# Monte Carlo Estimate of Thompson Sampling/BRAR Arm 1 Allocation Probability
# This helper function estimates the probability of selecting Arm 1 under TS
# for a two-arm binary outcome case using Monte Carlo simulations.
posterior_bin_sim = function(M = 10000, alphas, betas)
{
  # Simulate M samples from the posterior Beta distribution for Arm 1
  arm1 = stats::rbeta(M, alphas[1], betas[1])
  # Simulate M samples from the posterior Beta distribution for Arm 2
  arm2 = stats::rbeta(M, alphas[2], betas[2])

  # Calculate the proportion of times Arm 1's samples are greater than Arm 2's
  # This proportion is the Monte Carlo estimate of the probability of selecting Arm 1.
  ap_arm1 = mean(arm1 > arm2)
  ap_arm2 = 1 - ap_arm1
  return(c(ap_arm1, ap_arm2))
}

# Calculate the Exact Posterior Probability of Superiority
# Implements the closed-form Bayesian solution to calculate the probability
# that Treatment 2 is superior to Treatment 1 for binary outcomes.
posterior_bin_exact = function(alphas, betas)
{
  alpha_1 = alphas[1]
  beta_1 = betas[1]
  alpha_2 = alphas[2]
  beta_2 = betas[2]

  # This is the closed-form solution for P(T1 > T2)
  ap_arm2 = 0
  for (iii in 0:(alpha_2 - 1))
  {
    ap_arm2 = ap_arm2 +
      exp(lbeta(alpha_1 + iii, beta_1 + beta_2) - log(beta_2 + iii) - lbeta(1 + iii, beta_2) - lbeta(alpha_1, beta_1))
  }

  ap_arm1 = 1 - ap_arm2
  ap_arm2 = ap_arm2

  return(c(ap_arm1, ap_arm2))
}


# Helper function to estimate Thompson Sampling allocation probabilities for Normal outcomes
# with known population variance (Normal-Normal conjugate model).
# This function calculates allocation probabilities by performing Monte Carlo
# simulations from the posterior distributions of the means for each arm and
# determining which arm's sampled mean is the highest.
posterior_norm_sim = function(M = 10000, means, sds)
{
  arm1 = stats::rnorm(M, means[1], sds[1])
  arm2 = stats::rnorm(M, means[2], sds[2])

  # Calculate the proportion of times Arm 1's samples are greater than Arm 2's
  # This proportion is the Monte Carlo estimate of the probability of selecting Arm 1.
  ap_arm1 = mean(arm1 > arm2)
  ap_arm2 = 1 - ap_arm1
  return(c(ap_arm1, ap_arm2))
}


# Helper function to estimate Thompson Sampling allocation probabilities for Normal outcomes
# with known population variance (Normal-Normal conjugate model).
# This function calculates allocation probabilities exactly.
posterior_norm_exact = function(means, sds)
{
  # Calculate the posterior probability that arm 1 is better than arm 2.
  ap_arm1 = stats::pnorm(((means[1] - means[2]) / sqrt(sds[1]^2 + sds[2]^2)))
  ap_arm2 = 1 - ap_arm1
  return(c(ap_arm1, ap_arm2))
}

# Helper function to estimate Thompson Sampling allocation probabilities for Normal outcomes
# with unknown population variance (Normal-Inverse-Gamma prior).
# This function calculates allocation probabilities by performing Monte Carlo
# simulations from the posterior distributions of the means for each arm and
# determining which arm's sampled mean is the highest.
posterior_norm_unknownvar_sim <- function(M = 10000, mu_n, kappa_n, alpha_n, beta_n) {
  K = length(mu_n)   # number of arms
  samples = matrix(NA, nrow = M, ncol = K)

  # Sample posterior means for each arm
  for (k in 1:K) {
    scale_k = sqrt(beta_n[k] / (alpha_n[k] * kappa_n[k]))
    df_k = 2 * alpha_n[k]
    samples[, k] = mu_n[k] + scale_k * stats::rt(M, df = df_k)
  }

  # For each Monte Carlo draw, find which arm has the highest sampled mean
  winners = max.col(samples, ties.method = "first")

  # Compute allocation probabilities = frequency each arm is best
  ap = tabulate(winners, nbins = K) / M
  return(ap)
}



# Helper function to estimate Thompson Sampling allocation probabilities for
# exponential outcomes. This function calculates allocation probabilities by
# performing Monte Carlo simulations from the posterior distributions of the
# means for each arm and determining which arm's sampled mean is the highest.
posterior_exp_sim <- function(M = 10000, shapes, rates) {
  arm1 <- rgamma(M, shape = shapes[1], rate = rates[1])
  arm2 <- rgamma(M, shape = shapes[2], rate = rates[2])
  ap_arm1 <- mean(arm1 > arm2)
  ap_arm2 <- 1 - ap_arm1
  return(c(ap_arm1, ap_arm2))
}

# Helper function to calculate Thompson Sampling allocation probabilities for
# exponential outcomes with a Gamma prior. This function calculates allocation
# probabilities exactly, but it only works for integer alpha and beta parameters.
posterior_exp_exact <- function(shapes, rates) {
  a1 <- shapes[1]
  a2 <- shapes[2]
  b1 <- rates[1]
  b2 <- rates[2]
  prob <- 0
  for (k in 0:(a1-1)) {
    prob <- prob + choose(a2 + k - 1, k) *
      (b1 / (b1 + b2))^k *
      (b2 / (b1 + b2))^a2
  }
  ap_arm1 <- prob
  ap_arm2 <- 1 - prob
  return(c(ap_arm1, ap_arm2))
}

