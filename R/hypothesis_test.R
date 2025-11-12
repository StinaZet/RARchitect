# Randomization-based hypothesis test.
randomization_test_brar <- function(trial_data, priors, blocksize, postprobmethod,
                                    multiarm_method, tuning, clipping, randmethod,
                                    urn_alpha, burnin, alternative, arms_to_test, B = 10000) {

  # Extract things from the dataset.
  N = length(trial_data$Outcome)
  blocks = max(trial_data$Batch)
  control_outcome = trial_data$Outcome[trial_data$Arm == 1]
  n_control = length(control_outcome)

  # For storing.
  perm_stats = matrix(NA, nrow = B, ncol = max(trial_data$Arm) - 1)
  p_values = numeric(length(arms_to_test))


  # Loop over all the re-rendomizations.
  for (bbb in 1:B)
  {
    # Create the permuted dataset.
    trial_perm = trial_data

    # Start by re-randomizing the burnin.
    if (burnin > 0 && burnin <= N) {
      burnin_idx = 1:burnin

      # Generate roughly balanced allocation
      full_cycles = floor(burnin / arms) # Number of full cycles
      remainder = burnin %% arms # Leftover participants

      # Repeat each arm for full cycles
      arm_assignments = rep(1:arms, times = full_cycles)

      # Add remaining participants randomly among the arms
      if (remainder > 0) {
        arm_assignments = c(arm_assignments, sample(1:arms, remainder))
      }

      # Shuffle to avoid any ordering bias.
      # Change the treatment assigments in the permuted dataset.
      trial_perm$Arm[burnin_idx][burnin_idx] = sample(arm_assignments, burnin)

      # Update the prior parameters.
      current_alpha_params = priors[1, ]
      current_beta_params = priors[2, ]

      for(i in 1:arms) {
        priors[1, i] = current_alpha_params[i] + sum(trial_perm$Outcome[trial_perm$Arm[burnin_idx] == i])
        priors[2, i] = current_beta_params[i] + sum(trial_perm$Arm[burnin_idx] == i) - sum(trial_perm$Outcome[trial_perm$Arm[burnin_idx] == i])
      }
    }

    # Re-randomize the rest of the trial, block by block.
    for (lll in 1:blocks)
    {
      index1 = 1:burnin + lll
      index2 = burnin + lll: burnin + blocksize + lll
      trial_perm$Arm[index2] = brar_randomization(trial_data = trial_perm[index1, ], priors = priors,
                                                  blocksize = blocksize, postprobmethod = postprobmethod,
                                                  multiarm_method = multiarm_method, tuning = tuning,
                                                  clipping = clipping, randmethod = randmethod,
                                                  urn_alpha = urn_alpha, return = "allocations")
    }

    for (jjj in 2:max[trial_data$Arms])
    {
      # For the Wald statistic for the re-randomized datasets for all arms.
      p1 = mean(trial_perm$Outcome[trial_perm$Arm == jjj])
      p0 = mean(trial_perm$Outcome[trial_perm$Arm == 1])
      n1 = length(trial_perm$Arm[trial_perm$Arm==jjj])
      n0 = length(trial_perm$Arm[trial_perm$Arm==1])
      se = sqrt(p1 * (1 - p1) / n_1 + p0 * (1 - p0) / n_0)

      perm_stats[bbb, jjj - 1] = (p1 - p0) / se
    }
  }

  for (i in seq_along(arms_to_test))
  {
    k = arms_to_test[i]
    exp_outcome = trial_data$Outcome[trial_data$Arm == k]
    n_exp = length(exp_outcome)

    p1obs = mean(exp_outcome)
    p0obs = mean(control_outcome)
    seObs = sqrt(p1obs * (1 - p1obs) / n_exp + p0obs * (1 - p0obs) / n_control)

    p_obs = (p1obs - p0obs) / seObs

    if (alternative == "greater") {
      p_values[i] = mean(perm_stats[,i] >= p_obs)
    } else if(alternative == "less") {
      p_values[i] = mean(perm_stats[,i] <= p_obs)
    } else if (alternative == "two.sided") {
      p_values[i] = mean(abs(perm_stats[,i]) >= abs(p_obs))
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
