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
  ap_arm1 = sum(arm1 > arm2) / M
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
  total_prob_T2_best = 0
  for (i in 0:(alpha_2 - 1)) {
    total_prob_T2_best = total_prob_T2_best +
      exp(lbeta(alpha_1 + i, beta_1 + beta_2) - log(beta_2 + i) - lbeta(1 + i, beta_2) - lbeta(alpha_1, beta_1))
  }

  prob_T1_is_best = 1 - total_prob_T2_best
  prob_T2_is_best = total_prob_T2_best

  return(c(prob_T1_is_best,prob_T2_is_best))
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
  ap_arm1 = sum(arm1 > arm2) / M
  ap_arm2 = 1 - ap_arm1
  return(c(ap_arm1, ap_arm2))
}
