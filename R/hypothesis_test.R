# Randomization-based hypothesis test.
randomization_test_brar <- function(
    trial_data, priors, blocksize,
    direction = c("higher", "lower"),
    postprobmethod = c("simulation", "exact"),
    multiarm_method = c("top2", "fixed"),
    tuning = 1, clipping = 0, urn_alpha = 0,
    randmethod = c("coin", "block", "urn"),
    arms_to_test = 2:ncol(priors),
    B = 10000,
    alternative = c("greater", "less", "two.sided")
) {
  direction <- match.arg(direction)
  postprobmethod <- match.arg(postprobmethod)
  multiarm_method <- match.arg(multiarm_method)
  randmethod <- match.arg(randmethod)
  alternative <- match.arg(alternative)

  control_outcome <- trial_data$Outcome[trial_data$Arm == 1]
  n_control <- length(control_outcome)

  p_values <- numeric(length(arms_to_test))

  for (i in seq_along(arms_to_test)) {
    k <- arms_to_test[i]
    exp_outcome <- trial_data$Outcome[trial_data$Arm == k]
    n_exp <- length(exp_outcome)

    # Observed test statistic (difference in proportions)
    p_obs <- mean(exp_outcome) - mean(control_outcome)

    # --- Generate null distribution ---
    perm_stats <- numeric(B)
    for (b in 1:B) {
      # Generate BRAR sequence for this block configuration
      trial_data_perm <- trial_data
      trial_data_perm$Arm <- generate_brar_randomization(
        trial_data, priors, blocksize,
        direction, postprobmethod, multiarm_method,
        tuning, clipping, urn_alpha, randmethod,
        return = "allocations"
      )
      # Recompute outcome sums for control vs k
      perm_exp <- trial_data_perm$Outcome[trial_data_perm$Arm == k]
      perm_control <- trial_data_perm$Outcome[trial_data_perm$Arm == 1]

      # Compute permuted test statistic
      perm_stats[b] <- mean(perm_exp) - mean(perm_control)
    }

    # --- Compute p-value ---
    if (alternative == "greater") {
      p_values[i] <- (sum(perm_stats >= p_obs) + 1) / (B + 1)
    } else if (alternative == "less") {
      p_values[i] <- (sum(perm_stats <= p_obs) + 1) / (B + 1)
    } else {
      p_values[i] <- (sum(abs(perm_stats) >= abs(p_obs)) + 1) / (B + 1)
    }
  }

  return(p_values)
}



# Monte Carlo Simulation Test for Multiple Arms
monte_carlo_test_brar <- function(y_c, n_c, y_e_list, n_e_list,
                                  alternative = "two.sided", B = 10000) {
  # y_e_list: list of observed successes for each experimental arm
  # n_e_list: list of number of patients per experimental arm
  n_arms <- length(y_e_list)
  p_values <- numeric(n_arms)

  pooled_y <- y_c + unlist(y_e_list)
  pooled_n <- n_c + unlist(n_e_list)
  pooled_p <- pooled_y / pooled_n

  for (i in seq_len(n_arms)) {
    n_e <- n_e_list[[i]]
    y_e <- y_e_list[[i]]

    # Observed statistic
    T_obs <- y_e / n_e - y_c / n_c

    # Monte Carlo simulation
    y_c_sim <- stats::rbinom(B, n_c, pooled_p[i])
    y_e_sim <- stats::rbinom(B, n_e, pooled_p[i])
    T_sim <- y_e_sim / n_e - y_c_sim / n_c

    # p-value
    if (alternative == "greater") {
      p_values[i] <- (sum(T_sim >= T_obs) + 1) / (B + 1)
    } else if (alternative == "less") {
      p_values[i] <- (sum(T_sim <= T_obs) + 1) / (B + 1)
    } else { # two-sided
      p_values[i] <- (sum(abs(T_sim) >= abs(T_obs)) + 1) / (B + 1)
    }
  }
  return(p_values)
}



# The AP test for binary data with BRAR.
ap_test_brar <- function(trial_data, direction = "higher",
                         multiple_tests = FALSE, onesided = TRUE, B = 10000,
                         burnin = 0, blocksize = 1,
                         priors, postprobmethod = "simulation",
                         multiarm_method = "top2", randmethod = "coin",
                         urn_alpha = 0) {

  arms <- length(unique(trial_data$Arm))

  # Determine arms to test
  if (multiple_tests) {
    arms_to_test <- 2:arms
  } else {
    exp_arm_estimates <- tapply(trial_data$Outcome, trial_data$Arm, mean)[-1]
    if (direction == "higher") {
      best_exp <- which.max(exp_arm_estimates) + 1
    } else {
      best_exp <- which.min(exp_arm_estimates) + 1
    }
    arms_to_test <- best_exp
  }

  # Pooled mean across arms (excluding burn-in)
  if (burnin > 0) {
    pooled_mean <- mean(trial_data$Outcome[-(1:burnin)])
  } else {
    pooled_mean <- mean(trial_data$Outcome)
  }

  # Store observed effect estimates
  obs_effects <- sapply(arms_to_test, function(k) {
    mean(trial_data$Outcome[trial_data$Arm == k]) -
      mean(trial_data$Outcome[trial_data$Arm == 1])
  })

  # Simulate null distribution
  null_matrix <- matrix(NA, nrow = B, ncol = length(arms_to_test))

  for (b in 1:B) {
    sim_trial <- .simulate_brar_trial_binary(
      direction = direction, arms = arms, N = nrow(trial_data),
      blocksize = blocksize, priors = priors, modelpar = rep(pooled_mean, arms),
      tuning = 1, clipping = 0, burnin = burnin,
      postprobmethod = postprobmethod,
      randmethod = randmethod, urn_alpha = urn_alpha,
      multiarm_method = multiarm_method
    )
    null_matrix[b, ] <- sapply(arms_to_test, function(k) {
      mean(sim_trial$Outcome[sim_trial$Arm == k]) -
        mean(sim_trial$Outcome[sim_trial$Arm == 1])
    })
  }

  # Compute p-values
  p_values <- numeric(length(arms_to_test))
  for (i in seq_along(arms_to_test)) {
    null_dist <- null_matrix[, i]
    obs <- obs_effects[i]
    if (onesided) {
      if (direction == "higher") {
        p_values[i] <- mean(null_dist >= obs)
      } else {
        p_values[i] <- mean(null_dist <= obs)
      }
    } else {
      # Two-sided: use absolute deviations from mean
      p_values[i] <- mean(abs(null_dist) >= abs(obs))
    }
  }

  names(p_values) <- paste0("Arm", arms_to_test)
  return(p_values)
}
