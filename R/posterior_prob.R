# Monte Carlo Estimate of Thompson Sampling/BRAR Arm 1 Allocation Probability
# This helper function estimates the probability of selecting Arm 1 under TS
# for a two-arm binary outcome case using Monte Carlo simulations.
posterior_bin_sim <- function(M = 10000, alphas, betas, direction = c("higher", "lower")) {
  # Number of arms.
  K = length(alphas)

  # Simulate M posterior samples for each arm.
  # Create an M x K matrix where each column is rbeta(M, alpha_k, beta_k)
  samples = replicate(K, stats::rbeta(M, alphas, betas))

  if (direction == "lower") {
    samples = -samples
  }

  # For each simulation (row), find which arm had the largest sampled value
  best_arm = max.col(samples, ties.method = "random")

  # Compute selection probability for each arm
  ap = tabulate(best_arm, nbins = K) / M
  return(ap)
}

# Calculate the Exact Posterior Probability of Superiority
# Implements the closed-form Bayesian solution to calculate the probability
# that Treatment 2 is superior to Treatment 1 for binary outcomes.
posterior_bin_exact = function(alphas, betas, direction = c("higher", "lower"))
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

  if (direction == "lower") {
    tmp = ap_arm1
    ap_arm1 = ap_arm2
    ap_arm2 = tmp
  }

  return(c(ap_arm1, ap_arm2))
}


# Helper function to estimate Thompson Sampling allocation probabilities for Normal outcomes
# with known population variance (Normal-Normal conjugate model).
# This function calculates allocation probabilities by performing Monte Carlo
# simulations from the posterior distributions of the means for each arm and
# determining which arm's sampled mean is the highest.
posterior_norm_sim <- function(M = 10000, means, sds, direction = c("higher", "lower")) {
  # Number of arms.
  K = length(means)

  # Simulate M posterior samples for each arm
  samples = replicate(K, stats::rnorm(M, means, sds))

  if (direction == "lower") {
    samples = -samples
  }

  # For each simulation (row), find which arm had the highest sampled value
  best_arm = max.col(samples, ties.method = "random")

  # Compute selection probability for each arm
  ap = tabulate(best_arm, nbins = K) / M
  return(ap)
}

# Helper function to estimate Thompson Sampling allocation probabilities for Normal outcomes
# with known population variance (Normal-Normal conjugate model).
# This function calculates allocation probabilities exactly.
posterior_norm_exact = function(means, sds, direction = c("higher", "lower"))
{
  # Calculate the posterior probability that arm 1 is better than arm 2.
  ap_arm1 = stats::pnorm(((means[1] - means[2]) / sqrt(sds[1]^2 + sds[2]^2)))
  ap_arm2 = 1 - ap_arm1

  if (direction == "lower") {
    tmp = ap_arm1
    ap_arm1 = ap_arm2
    ap_arm2 = tmp
  }

  return(c(ap_arm1, ap_arm2))
}

# Helper function to estimate Thompson Sampling allocation probabilities for Normal outcomes
# with unknown population variance (Normal-Inverse-Gamma prior).
# This function calculates allocation probabilities by performing Monte Carlo
# simulations from the posterior distributions of the means for each arm and
# determining which arm's sampled mean is the highest.
posterior_norm_unknownvar_sim <- function(M = 10000, mu_n, kappa_n, alpha_n, beta_n, direction = c("higher", "lower")) {
  K = length(mu_n)   # number of arms
  samples = matrix(NA, nrow = M, ncol = K)

  # Sample posterior means for each arm
  for (k in 1:K) {
    scale_k = sqrt(beta_n[k] / (alpha_n[k] * kappa_n[k]))
    df_k = 2 * alpha_n[k]
    samples[, k] = mu_n[k] + scale_k * stats::rt(M, df = df_k)
  }

  if (direction == "lower") {
    samples = -samples
  }

  # For each Monte Carlo draw, find which arm has the highest sampled mean
  winners = max.col(samples, ties.method = "random")

  # Compute allocation probabilities = frequency each arm is best
  ap = tabulate(winners, nbins = K) / M
  return(ap)
}



# Helper function to estimate Thompson Sampling allocation probabilities for
# exponential outcomes. This function calculates allocation probabilities by
# performing Monte Carlo simulations from the posterior distributions of the
# means for each arm and determining which arm's sampled mean is the highest.
posterior_exp_sim <- function(M = 10000, shapes, rates, direction = c("higher", "lower")) {
  # Number of treatment arms.
  K = length(shapes)

  # Simulate M posterior samples for each arm
  samples = replicate(K, stats::rgamma(M, shape = shapes, rate = rates))

  if (direction == "lower") {
    samples = -samples
  }

  # For each simulation (row), find which arm had the highest sampled value
  best_arm = max.col(samples, ties.method = "random")

  # Compute posterior selection probabilities
  ap = tabulate(best_arm, nbins = K) / M
  return(ap)
}

# Helper function to calculate Thompson Sampling allocation probabilities for
# exponential outcomes with a Gamma prior. This function calculates allocation
# probabilities exactly, but it only works for integer alpha and beta parameters.
posterior_exp_exact <- function(shapes, rates, direction = c("higher", "lower")) {
  a1 = shapes[1]
  a2 = shapes[2]
  b1 = rates[1]
  b2 = rates[2]
  prob = 0
  for (k in 0:(a1-1)) {
    prob = prob + choose(a2 + k - 1, k) *
      (b1 / (b1 + b2))^k *
      (b2 / (b1 + b2))^a2
  }
  ap_arm1 = prob
  ap_arm2 = 1 - prob

  if (direction == "lower") {
    tmp = ap_arm1
    ap_arm1 = ap_arm2
    ap_arm2 = tmp
  }

  return(c(ap_arm1, ap_arm2))
}


# Helper function to calculate Top 2 Thompson Sampling allocation probabilities.
alloc_probs_T2TS <- function(alloc_probs_raw, beta) {
  K = length(alloc_probs_raw)
  pi_new = numeric(K)

  for (k in seq_len(K)) {
    # Term for all k' ≠ k
    k_others = seq_len(K)[-k]
    sum_term = sum((alloc_probs_raw[k] / (1 - alloc_probs_raw[k_others])) * alloc_probs_raw[k_others])

    # apply formula
    pi_new[k] = beta * alloc_probs_raw[k] + (1 - beta) * sum_term
  }

  # Normalize to ensure they sum to 1 (numerical stability)
  pi_new = pi_new / sum(pi_new)

  return(pi_new)
}
