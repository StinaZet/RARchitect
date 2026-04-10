# Randomization-based hypothesis test.
randomization_test_brar <- function(trial_data, priors, blocksize, postprobmethod,
                                    multiarm_method, tuning, clipping, randmethod,
                                    urn_alpha, burnin, direction, alternative,
                                    arms_to_test, B, ...) {

  # Extract things from the dataset.
  N = length(trial_data$Outcome)
  blocks = max(trial_data$Batch)
  control_outcome = trial_data$Outcome[trial_data$Arm == 1]
  n_control = length(control_outcome)
  arms = max(trial_data$Arm)

  # For storing.
  perm_stats = matrix(NA, nrow = B, ncol = arms - 1)
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
      # Change the treatment assignments in the permuted dataset.
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
    for (lll in 0:blocks)
    {
      index1 = 1:(burnin + lll * blocksize)
      index2 = (burnin + lll * blocksize):(burnin + blocksize + lll * blocksize)
      trial_perm$Arm[index2] = brar_randomization(trial_data = trial_perm[index1, ], priors = priors,
                                                  blocksize = blocksize, postprobmethod = postprobmethod,
                                                  multiarm_method = multiarm_method, tuning = tuning,
                                                  clipping = clipping, randmethod = randmethod, direction = direction,
                                                  urn_alpha = urn_alpha, return = "allocations")
    }

    # Calculate the test statistic for all arms (even if the arm will not be tested).
    for (jjj in 2:arms)
    {
      # For the Wald statistic for the re-randomized datasets for all arms.
      p1 = mean(trial_perm$Outcome[trial_perm$Arm == jjj])
      p0 = mean(trial_perm$Outcome[trial_perm$Arm == 1])
      n1 = length(trial_perm$Arm[trial_perm$Arm==jjj])
      n0 = length(trial_perm$Arm[trial_perm$Arm==1])
      se = sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)

      perm_stats[bbb, jjj - 1] = (p1 - p0) / se
    }
  }

  # Compute the values of the test statistics and p-values for the observed data.
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
# First function is for finding the null distribution.
monte_carlo_null_brar <- function(trial_data, priors, blocksize,
                                 postprobmethod, multiarm_method,
                                 tuning, clipping, randmethod,
                                 urn_alpha, burnin, direction,
                                 B, alpha) {

  # --- Extract trial structure ---
  N = length(trial_data$Outcome)
  arms = max(trial_data$Arm)

  # --- Storage ---
  mc_stats = matrix(NA, nrow = B, ncol = arms - 1)

  # --- Simulate under null hypothesis ---
  for (b in seq_len(B)) {
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
      recruitment_rate = 100000,
      observation_delay = 0
    )

    # --- Compute Wald test statistic for each experimental arm ---
    for (k in 2:arms) {
      p1 = mean(mc_data$Outcome[mc_data$Arm == k])
      p0 = mean(mc_data$Outcome[mc_data$Arm == 1])
      n1 = sum(mc_data$Arm == k)
      n0 = sum(mc_data$Arm == 1)
      se = sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
      mc_stats[b, k - 1] = (p1 - p0) / se
    }
  }

  # --- Compute null critical values for chosen alpha ---
  crit_values = list(
    greater = apply(mc_stats, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE),
    less = apply(mc_stats, 2, stats::quantile, probs = alpha, na.rm = TRUE),
    two.sided = apply(mc_stats, 2, function(x)
      stats::quantile(abs(x), probs = 1 - alpha / 2, na.rm = TRUE))
  )

  return(list(
    mc_stats = mc_stats,
    crit_values = crit_values,
    alpha = alpha
  ))
}

# Second function is for performing the test.
monte_carlo_test_brar <- function(trial_data, arms_to_test, critval = NULL,
                                 alternative, null_distribution = NULL) {


  # If the critical value is given, just compute the test decision.
  if(is.numeric(critval))
  {
    test_results = numeric(length(arms_to_test))

    control_outcome = trial_data$Outcome[trial_data$Arm == 1]
    n_control = length(control_outcome)

    for (i in seq_along(arms_to_test)) {
      # Compute the observed test statistic.
      k = arms_to_test[i]
      exp_outcome = trial_data$Outcome[trial_data$Arm == k]
      n_exp = length(exp_outcome)

      p1obs = mean(exp_outcome)
      p0obs = mean(control_outcome)
      seObs = sqrt(p1obs * (1 - p1obs) / n_exp + p0obs * (1 - p0obs) / n_control)
      T_obs = (p1obs - p0obs) / seObs

      # --- Compute test results from empirical critical values ---
      if (alternative == "greater") {
        test_results[i] = ifelse(T_obs >= critval[i], 1, 0)
      } else if (alternative == "less") {
        test_results[i] = ifelse(T_obs <= critval[i], 1, 0)
      } else if (alternative == "two.sided") {
        test_results[i] = ifelse(abs(T_obs) >= critval[i], 1, 0)
      }
    }

    names(test_results) = paste0("Arm", arms_to_test)

    return(list(test_results = test_results, critval = critval))
  } else{ # Compute p-values with the full null distribution.
    mc_stats = null_distribution$mc_stats
    crit_values = null_distribution$crit_values
    alpha = null_distribution$alpha

    # --- Extract control arm data ---
    control_outcome = trial_data$Outcome[trial_data$Arm == 1]
    n_control = length(control_outcome)

    # --- Compute observed test statistics ---
    p_values = numeric(length(arms_to_test))
    critval = numeric(length(arms_to_test))

    for (i in seq_along(arms_to_test)) {
      k = arms_to_test[i]
      exp_outcome = trial_data$Outcome[trial_data$Arm == k]
      n_exp = length(exp_outcome)

      p1obs = mean(exp_outcome)
      p0obs = mean(control_outcome)
      seObs = sqrt(p1obs * (1 - p1obs) / n_exp + p0obs * (1 - p0obs) / n_control)
      T_obs = (p1obs - p0obs) / seObs

      # --- Compute p-values from empirical null distribution ---
      if (alternative == "greater") {
        p_values[i] = mean(mc_stats[, i] >= T_obs)
        critval[i] = crit_values$greater[i]
      } else if (alternative == "less") {
        p_values[i] = mean(mc_stats[, i] <= T_obs)
        critval[i] = crit_values$less[i]
      } else if (alternative == "two.sided") {
        p_values[i] = mean(abs(mc_stats[, i]) >= abs(T_obs))
        critval[i] = crit_values$two.sided[i]
      }
    }

    names(p_values) = paste0("Arm", arms_to_test)

    return(list(p_values = p_values, critval = critval, alpha = alpha))
  }
}


# The AP test for binary data with BRAR.
# First function is for finding the null distribution.
ap_null_brar <- function(trial_data, priors, blocksize,
                         postprobmethod, multiarm_method,
                         tuning, clipping, randmethod,
                         urn_alpha, burnin, direction,
                         B, alpha) {

  # --- Extract trial structure ---
  N = length(trial_data$Outcome)
  arms = max(trial_data$Arm)

  # --- Storage ---
  ap_stats = matrix(NA, nrow = B, ncol = arms - 1)

  # --- Simulate under null hypothesis ---
  for (b in seq_len(B)) {
    ap_data = simulate_brar_trial(
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
      recruitment_rate = 100000,
      observation_delay = 0
    )

    # Remove the burn-in and then extract the first individual in every block.
    # Remove the first three columns (batch, treatment time, and outcome time).
    ap_data = ap_data[-(1:burnin), -(1:3)][seq(1, N - burnin, by = blocksize), ]

    # --- Compute AP test statistic for each arm ---
    for (k in 1:(arms - 1)) {

      # Value for equal randomization.
      ER = 1 / arms

      # AP test statistics if direction is "higher".
      if(direction == "higher"){
        ap_stats[b, k] = sum(ap_data[, 3 + k] > ER) # The 3 is because the first two columns are treatment and outcome, and the third is the allocation probabilities for the control arm.
      } else if(direction == "lower"){
        ap_stats[b, k] = sum(ap_data[, 3 + k] < ER) # The 3 is because the first two columns are treatment and outcome, and the third is the allocation probabilities for the control arm.
      }
    }
  }

  # --- Compute null critical values for chosen alpha ---
  crit_values = list(
    if(direction == "higher"){
      greater = apply(ap_stats, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
    } else if(direction == "lower"){
      less = apply(ap_stats, 2, stats::quantile, probs = alpha, na.rm = TRUE)
    }
  )

  return(list(
    ap_stats = ap_stats,
    crit_values = crit_values,
    alpha = alpha
  ))
}


# Second function is for performing the test.
ap_test_brar <- function(trial_data, burnin, blocksize, arms_to_test, critval = NULL,
                         alternative, null_distribution = NULL) {


  # ER probability to compare allocation probabilities against.
  arms = max(trial_data$Arm)
  ER = 1 / arms

  # Prepare the data for AP test.
  N = length(trial_data[,1])

  # Remove the duplicates in term of allocation probabilities
  trial_data = trial_data[-(1:burnin), ][seq(1, N - burnin, by = blocksize), ]

  # If the critical value is given, just compute the test decision.
  if(is.numeric(critval))
  {
    test_results = numeric(length(arms_to_test))

    for (i in seq_along(arms_to_test)) {
      # Compute the observed test statistic.
      if (alternative == "greater") {
        ap_obs = sum(trial_data[, 3 + i] > ER) # The 3 is because the first two columns are treatment and outcome, and the third is the allocation probabilities for the control arm.
        test_results[i] = ifelse(ap_obs > critval[i], 1, 0)
      } else if (alternative == "less") {
        ap_obs = sum(trial_data[, 3 + i] < ER) # The 3 is because the first two columns are treatment and outcome, and the third is the allocation probabilities for the control arm.
        test_results[i] = ifelse(ap_obs < critval[i], 1, 0)
      }
    }

    names(test_results) = paste0("Arm", arms_to_test)

    return(list(test_results = test_results, critval = critval))
  } else{ # Compute p-values with the full null distribution.
    ap_stats = null_distribution$ap_stats
    crit_values = null_distribution$crit_values
    alpha = null_distribution$alpha


    # --- Compute observed test statistics ---
    p_values = numeric(length(arms_to_test))
    critval = numeric(length(arms_to_test))

    for (i in seq_along(arms_to_test)) {
      k = arms_to_test[i]
      if (alternative == "greater") {
        ap_obs = sum(trial_data[, 3 + i] > ER) # The 3 is because the first two columns are treatment and outcome, and the third is the allocation probabilities for the control arm.
      } else if (alternative == "less") {
        ap_obs = sum(trial_data[, 3 + i] < ER) # The 3 is because the first two columns are treatment and outcome, and the third is the allocation probabilities for the control arm.
      }

      # --- Compute p-values from empirical null distribution ---
      if (alternative == "greater") {
        p_values[i] = mean(ap_stats[, i] > ap_obs)
        if(p_values[i] == 0) p_values[i] = mean(ap_stats[, i] > ap_obs - 1)
        critval[i] = crit_values$greater[i]
      } else if (alternative == "less") {
        p_values[i] = mean(ap_stats[, i] < ap_obs)
        if(p_values[i] == 0) p_values[i] = mean(ap_stats[, i] < ap_obs + 1)
        critval[i] = crit_values$less[i]
      }
    }

    names(p_values) = paste0("Arm", arms_to_test)

    return(list(p_values = p_values, critval = critval, alpha = alpha))
  }
}


