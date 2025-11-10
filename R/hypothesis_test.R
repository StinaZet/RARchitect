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
monte_carlo_test_brar <- function(y_c, n_c, y_e_list, n_e_list,
                                  direction = "higher", B = 10000) {
  n_arms = length(y_e_list)
  p_values = numeric(n_arms)

  pooled_p = (y_c + unlist(y_e_list)) / (n_c + unlist(n_e_list))

  for (i in seq_len(n_arms)) {
    n_e = n_e_list[[i]]
    y_e = y_e_list[[i]]
    T_obs = y_e / n_e - y_c / n_c

    y_c_sim = rbinom(B, n_c, pooled_p[i])
    y_e_sim = rbinom(B, n_e, pooled_p[i])
    T_sim = y_e_sim / n_e - y_c_sim / n_c

    if (direction == "higher") {
      p_values[i] = mean(T_sim >= T_obs)
    } else {
      p_values[i] = mean(T_sim <= T_obs)
    }
  }

  names(p_values) = paste0("Arm", seq_len(n_arms))
  return(p_values)
}



# The AP test for binary data with BRAR.
ap_test_brar <- function(trial_data, priors, blocksize, direction = "higher",
                         multiple_tests = FALSE, B = 10000) {
  arms = length(unique(trial_data$Arm))
  if (multiple_tests) {
    arms_to_test = 2:arms
  } else {
    exp_means = tapply(trial_data$Outcome, trial_data$Arm, mean)[-1]
    if (direction == "higher") {
      best_exp = which.max(exp_means) + 1
    } else {
      best_exp = which.min(exp_means) + 1
    }
    arms_to_test = best_exp
  }

  pooled_mean = mean(trial_data$Outcome)
  obs_effects = sapply(arms_to_test, function(k) mean(trial_data$Outcome[trial_data$Arm==k]) -
                          mean(trial_data$Outcome[trial_data$Arm==1]))

  null_matrix = matrix(NA, nrow=B, ncol=length(arms_to_test))
  for (b in 1:B) {
    sim_trial = .simulate_brar_trial_binary(arms=arms, N=nrow(trial_data),
                                             blocksize=blocksize, priors=priors,
                                             modelpar=rep(pooled_mean, arms))
    null_matrix[b, ] = sapply(arms_to_test, function(k) mean(sim_trial$Outcome[sim_trial$Arm==k]) -
                                 mean(sim_trial$Outcome[sim_trial$Arm==1]))
  }

  p_values = numeric(length(arms_to_test))
  for (i in seq_along(arms_to_test)) {
    if (direction == "higher") {
      p_values[i] = mean(null_matrix[,i] >= obs_effects[i])
    } else {
      p_values[i] = mean(null_matrix[,i] <= obs_effects[i])
    }
  }

  names(p_values) = paste0("Arm", arms_to_test)
  return(p_values)
}
