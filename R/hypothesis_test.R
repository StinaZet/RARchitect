# Randomization-based hypothesis test.
#' @title Perform a Randomization-Based Permutation Test

perform_permutation_test <- function(y_c, n_c, y_e, n_e,
                                      alternative = "two.sided", B = 10000) {

  # 1. Calculate the observed test statistic (Difference in Proportions)
  # Handle division by zero if n_c or n_e is 0 (though analyze_brar_trial checks this)
  if (n_e == 0 || n_c == 0) return(NA_real_)
  p_obs = y_e/n_e - y_c/n_c

  # 2. Pool the data (Outcomes under H0)
  n_total = n_c + n_e
  y_total = y_c + y_e

  outcomes = c(rep(1, y_total), rep(0, n_total - y_total))

  perm_stats = numeric(B)

  # 3. Permutation Loop
  for (i in 1:B) {
    perm_outcomes = sample(outcomes)
    y_e_perm = sum(perm_outcomes[1:n_e])
    y_c_perm = sum(perm_outcomes[(n_e + 1):n_total])

    perm_stats[i] = y_e_perm/n_e - y_c_perm/n_c
  }

  # 4. Calculate p-value
  if (alternative == "greater") {
    # One-sided test for p_e > p_c (positive difference)
    p_value = (sum(perm_stats >= p_obs) + 1) / (B + 1) # Add +1 for obs
  } else if (alternative == "less") {
    # One-sided test for p_e < p_c (negative difference)
    p_value = (sum(perm_stats <= p_obs) + 1) / (B + 1) # Add +1 for obs
  } else {
    # Two-sided p-value: P(|T_perm| >= |T_obs|)
    p_value = (sum(abs(perm_stats) >= abs(p_obs)) + 1) / (B + 1) # Add +1 for obs
  }

  return(p_value)
}

# Perform a Monte Carlo Test for Two Proportions
perform_mc_test <- function(y_c, n_c, y_e, n_e,
                             alternative = "two.sided", B = 10000) {

  # 1. Calculate the observed test statistic (Difference in Proportions)
  # This check is technically already done in the parent function
  if (n_e == 0 || n_c == 0) return(NA_real_)

  p_e_obs = y_e / n_e
  p_c_obs = y_c / n_c
  T_obs  = p_e_obs - p_c_obs

  # 2. Estimate common proportion under H0
  pooled_p = (y_c + y_e) / (n_c + n_e)

  # Handle edge case where pooled_p is 0 or 1 (no variance)
  if (pooled_p == 0 || pooled_p == 1) {
    # If T_obs is also 0, any sim will also be 0. p = 1.
    # If T_obs is not 0, it's impossible under H0. p = 0 (or 1/(B+1)).
    return(if (T_obs == 0) 1.0 else 1 / (B + 1))
  }

  # 3. Monte Carlo simulation under H0
  # Simulate the *number of successes* in each arm from the pooled prob
  # This is much faster than simulating individual patient outcomes
  y_e_sim = stats::rbinom(B, n_e, pooled_p)
  y_c_sim = stats::rbinom(B, n_c, pooled_p)

  # Statistic for simulated data
  T_sim = (y_e_sim / n_e) - (y_c_sim / n_c)

  # 4. Compute p-value
  if (alternative == "two.sided") {
    p_value = (sum(abs(T_sim) >= abs(T_obs)) + 1) / (B + 1)
  } else if (alternative == "greater") {
    p_value = (sum(T_sim >= T_obs) + 1) / (B + 1)
  } else { # "less"
    p_value = (sum(T_sim <= T_obs) + 1) / (B + 1)
  }

  return(p_value)
}
