simulate_brar_ci <- function(trial_data, priors, blocksize, postprobmethod,
                             multiarm_method, tuning, clipping, randmethod,
                             urn_alpha, burnin, arms_to_test, alpha, onesided,
                             direction, B) {

  # --- Extract trial structure ---
  arms = max(trial_data$Arm)
  N = length(trial_data$Outcome)

  # --- Storage: B × K matrix ---
  K = length(arms_to_test)
  T_sim = matrix(NA, nrow = B, ncol = K)

  # --- Bootstrap loop ---
  for (bbb in seq_len(B))
  {
    # Simulate bootstrap dataset
    mc_data = simulate_brar_trial(
      outcome_type = "binary",
      distribution = "bernoulli",
      arms = arms,
      N = N,
      blocksize = blocksize,
      priors = priors,
      modelpar = rep(mean(trial_data$Outcome), arms),  # pooled mean under H0
      direction = direction,
      randmethod = randmethod,
      tuning = tuning,
      clipping = clipping,
      postprobmethod = postprobmethod,
      multiarm_method = multiarm_method,
      recruitment_rate = 1e6,
      observation_delay = 0)

    # Extract simulated control outcomes
    y_c_sim = mc_data$Outcome[mc_data$Arm == 1]

    # Compute bootstrap statistics for each arm
    for (jjj in seq_along(arms_to_test))
    {
      k = arms_to_test[jjj]
      y_e_sim = mc_data$Outcome[mc_data$Arm == k]
      T_sim[bbb, jjj] = mean(y_e_sim) - mean(y_c_sim)
    }
  }

  # --- Compute CI bounds for each arm ---
  ci = matrix(NA, nrow = 2, ncol = K)  # row1=low, row2=high

  for (jjj in seq_along(arms_to_test))
  {
    # The arm to calculate the confidence interval for.
    col = T_sim[, jjj]

    # The different sets of confidence intervals.
    if (onesided)
    {
      if (direction == "higher") {
        ci[1, jjj] = stats::quantile(col, alpha, na.rm = TRUE)
        ci[2, jjj] = Inf
      } else {
        ci[1, jjj] = -Inf
        ci[2, jjj] = stats::quantile(col, 1 - alpha, na.rm = TRUE)
      }
    } else {  # two-sided
      ci[1, jjj] = stats::quantile(col, alpha / 2, na.rm = TRUE)
      ci[2, jjj] = stats::quantile(col, 1 - alpha / 2, na.rm = TRUE)
    }
  }

  colnames(ci) = paste0("Arm", arms_to_test)
  rownames(ci) = c("lower", "upper")

  return(ci)
}
