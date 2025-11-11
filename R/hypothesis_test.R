# Randomization-based hypothesis test.
randomization_test_brar <- function(trial_data, priors, blocksize, arms_to_test = 2:ncol(priors),
                                    direction = "higher", B = 10000) {
  control_outcome = trial_data$Outcome[trial_data$Arm == 1]
  n_control = length(control_outcome)
  p_values = numeric(length(arms_to_test))

  for (i in seq_along(arms_to_test)) {
    k = arms_to_test[i]
    exp_outcome = trial_data$Outcome[trial_data$Arm == k]
    n_exp = length(exp_outcome)

    p_obs = mean(exp_outcome) - mean(control_outcome)
    perm_stats = numeric(B)
    for (b in 1:B) {
      trial_perm = trial_data
      trial_perm$Arm = brar_randomization(trial_data, priors, blocksize, return = "allocations")
      perm_exp = trial_perm$Outcome[trial_perm$Arm == k]
      perm_control = trial_perm$Outcome[trial_perm$Arm == 1]
      perm_stats[b] = mean(perm_exp) - mean(perm_control)
    }

    if (direction == "higher") {
      p_values[i] = mean(perm_stats >= p_obs)
    } else {
      p_values[i] = mean(perm_stats <= p_obs)
    }
  }

  names(p_values) = paste0("Arm", arms_to_test)
  return(p_values)
}


# Monte Carlo Simulation Test for Multiple Arms
monte_carlo_test_brar <- function(y_c, n_c, y_e_list, n_e_list, alternative, B = 10000) {
  n_arms = length(y_e_list)
  p_values = numeric(n_arms)

  for (i in seq_len(n_arms)) {
    y_e = y_e_list[[i]]
    n_e = n_e_list[[i]]

    # Observed statistic (risk difference)
    T_obs = y_e / n_e - y_c / n_c

    # Pooled proportion under H0 (no difference)
    pooled_p = (y_c + y_e) / (n_c + n_e)

    # Simulate under H0 using pooled probability
    y_c_sim = stats::rbinom(B, n_c, pooled_p)
    y_e_sim = stats::rbinom(B, n_e, pooled_p)
    T_sim = y_e_sim / n_e - y_c_sim / n_c

    # Compute p-value depending on alternative
    if (alternative == "greater") {
      p_values[i] = mean(T_sim >= T_obs)
    } else if (alternative == "less") {
      p_values[i] = mean(T_sim <= T_obs)
    } else if (alternative == "two.sided") {
      # Two-sided: based on absolute deviation
      p_values[i] = mean(abs(T_sim) >= abs(T_obs))
    }
  }

  names(p_values) = paste0("Arm", seq_len(n_arms))
  return(p_values)
}




# The AP test for binary data with BRAR.
ap_test_brar <- function(
    trial_data, priors, blocksize,
    direction = "higher", multiple_tests = FALSE,
    onesided = TRUE, alpha = 0.05, B = 10000, modelpar,
    N = NULL,  randmethod = "coin", tuning = 1, clipping = 0,
    postprobmethod = "simulation", multiarm_method, ...
) {

  arms = length(unique(trial_data$Arm))
  if (is.null(N)) N = nrow(trial_data)

  # --- Determine which arms to test ---
  if (multiple_tests) {
    arms_to_test = 2:arms
  } else {
    # Select best experimental arm based on allocation frequency
    alloc_freq = table(trial_data$Arm) / nrow(trial_data)
    control_prob = alloc_freq[1]
    exp_probs = alloc_freq[-1]
    if (direction == "higher") {
      best_arm = which.max(exp_probs) + 1
    } else {
      best_arm = which.min(exp_probs) + 1
    }
    arms_to_test = best_arm
  }

  # --- Observed test statistic: how often each arm was favored ---
  # The test statistic = number of batches where arm's allocation prob > 1/arms
  favored_counts = numeric(length(arms_to_test))
  for (i in seq_along(arms_to_test)) {
    k = arms_to_test[i]
    favored_counts[i] = sum(tapply(trial_data$Arm == k, trial_data$Batch, mean) > 1 / arms)
  }

  # --- Simulate null distribution ---
  null_matrix = matrix(NA, nrow = B, ncol = length(arms_to_test))
  for (b in seq_len(B)) {
    sim_trial = simulate_brar_trial(outcome_type = c("binary"),
      distribution = c("bernoulli"),
      arms = arms, N = N, blocksize = blocksize,
      priors = priors, modelpar = rep(mean(trial_data$Outcome), arms),
      direction = direction, randmethod = randmethod,
      tuning = tuning, clipping = clipping,
      postprobmethod = postprobmethod,
      multiarm_method = multiarm_method,
      recruitment_rate = 100000,
      observation_delay = 0,...
    )

    for (i in seq_along(arms_to_test)) {
      k = arms_to_test[i]
      null_matrix[b, i] = sum(tapply(sim_trial$Arm == k, sim_trial$Batch, mean) > 1 / arms)
    }
  }

  # --- Compute conservative integer critical values and p-values ---
  p_values = numeric(length(arms_to_test))
  critical_values = integer(length(arms_to_test))

  for (i in seq_along(arms_to_test)) {
    obs_stat = favored_counts[i]
    null_dist = null_matrix[, i]
    unique_vals = sort(unique(null_dist))
    tail_probs = sapply(unique_vals, function(v) mean(null_dist >= v))
    low_tail_probs = sapply(unique_vals, function(v) mean(null_dist <= v))

    if (onesided) {
      # Conservative integer cutoff
      valid_vals = unique_vals[tail_probs <= alpha]
      if (length(valid_vals) > 0) {
        crit_val = max(valid_vals)
      } else {
        nonzero_vals = unique_vals[tail_probs > 0]
        crit_val = max(nonzero_vals)
      }

      p_val = mean(null_dist >= obs_stat)
    } else {
      # Two-sided: both tails (non-symmetric)
      valid_high = unique_vals[tail_probs <= alpha / 2]
      valid_low = unique_vals[low_tail_probs <= alpha / 2]

      if (length(valid_high) > 0) {
        crit_high = max(valid_high)
      } else {
        nonzero_high = unique_vals[tail_probs > 0]
        crit_high = max(nonzero_high)
      }
      if (length(valid_low) > 0) {
        crit_low = min(valid_low)
      } else {
        nonzero_low = unique_vals[low_tail_probs > 0]
        crit_low = min(nonzero_low)
      }

      crit_val = c(low = crit_low, high = crit_high)
      p_val = mean(null_dist <= obs_stat) + mean(null_dist >= obs_stat)
      p_val = min(p_val, 1)
    }

    p_values[i] = p_val
    critical_values[i] = if (onesided) crit_val else NA_integer_
  }

  names(p_values) = paste0("Arm", arms_to_test)
  return(list(
    p_values = p_values,
    critical_values = critical_values,
    test_statistics = favored_counts
  ))
}
