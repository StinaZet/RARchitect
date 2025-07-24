#' Simulate a Normal Response Adaptive Randomization (NRAR) Trial (Internal)
#'
#' This is an internal helper function for `simulate_brar_trial`. It simulates a
#' multi-arm clinical trial using Thompson Sampling with a Normal-Normal conjugate
#' model for normal outcomes (assuming known population variance for each arm).
#'
#' @param arms Numeric. Number of arms in the trial.
#' @param N Numeric. Total sample size for the trial.
#' @param blocksize Numeric. (Fixed) size of each block of participants.
#' @param priors Matrix. Normal prior parameters (mu0, tau0) for the mean of each arm.
#' @param modelpar Matrix. True mean and standard deviation for each normal arm.
#' @param tuning Numeric. Tuning parameter for allocation probabilities.
#' @param clipping Character or Numeric. Clipping parameter for allocation probabilities,
#'   can be numeric (fixed) or "adaptive".
#' @param burnin Numeric. Number of initial participants for burn-in.
#' @param ensure_all_arms_sampled Logical. If TRUE, ensures at least one arm is sampled per block.
#' @keywords internal
.simulate_brar_trial_normal <- function(arms = 2, N, blocksize,
                                        priors,
                                        modelpar, tuning = 1, clipping = 0, burnin = 0,
                                        ensure_all_arms_sampled = FALSE)
{
  # --- Input Validation and Setup ---
  # Only specific validation relevant to this internal function.
  # Broader validation is handled by the main `simulate_brar_trial` function.
  if (!is.matrix(modelpar) || nrow(modelpar) != 2 || ncol(modelpar) != arms) {
    stop("'modelpar' must be a 2-row matrix with 'arms' columns (first row: true means, second row: true sds).")
  }
  if (!is.matrix(priors) || nrow(priors) != 2 || ncol(priors) != arms) {
    stop("'priors' must be a 2-row matrix with 'arms' columns (first row: mu0, second row: tau0).")
  }
  if (is.numeric(clipping) && (clipping < 0 || clipping >= 1)) {
    stop("The numeric 'clipping' parameter must be between 0 and 1 (exclusive of 1).")
  }
  if (is.character(clipping) && clipping != "adaptive") {
    stop("The character 'clipping' parameter must be 'adaptive'.")
  }
  if (is.numeric(clipping) && clipping > 0 && (clipping * arms) > 1) {
    stop("Invalid numeric 'clipping' value: clipping * arms must be <= 1 to allow for consistent bounds across all arms and sum of probabilities to be 1.")
  }

  # Extract true parameters from modelpar for convenience
  true_means <- modelpar[1, ]
  true_sds <- modelpar[2, ]
  true_vars <- true_sds^2 # Convert true sds to true variances

  # Determine block sizes for each iteration, considering burn-in
  if (burnin > 0) {
    Nblocks <- 1 + floor((N - burnin) / blocksize)
    block_sizes <- c(burnin, rep(blocksize, Nblocks - 1))
  } else {
    Nblocks <- floor(N / blocksize)
    block_sizes <- rep(blocksize, Nblocks)
  }

  # Initialize vectors/matrix to store results for all N participants
  outcomes <- numeric(N)
  selected_arm <- numeric(N)
  batch_number <- numeric(N)
  allocation_probs_matrix <- matrix(NA, nrow = N, ncol = arms)
  colnames(allocation_probs_matrix) <- paste0("AlloProb_Arm", 1:arms)

  # Initialize Normal-Normal posterior parameters for each arm
  current_posterior_mu_params <- priors[1, ] # mu0
  current_posterior_tau2_params <- priors[2, ]^2 # tau0^2 (note: priors[2, ] should be tau0, not tau0^2)

  # Count of observations for each arm (n_k)
  n_k_counts <- rep(0, arms)
  # Sum of outcomes for each arm (sum_y_k)
  sum_y_k <- rep(0, arms)

  current_range_end <- 0

  # --- Simulation for burn-in period ---
  if (burnin > 0 && burnin <= N) {
    burnin_indices <- 1:burnin
    batch_number[burnin_indices] <- 1

    selected_arm[burnin_indices] <- blockrand::blockrand(burnin, arms, block.sizes = 1, levels = seq(1, arms, by = 1))$treatment

    outcomes[burnin_indices] <- rnorm(burnin,
                                      mean = true_means[selected_arm[burnin_indices]],
                                      sd = true_sds[selected_arm[burnin_indices]])

    allocation_probs_matrix[burnin_indices, ] <- matrix(
      rep(1 / arms, each = burnin),
      ncol = arms, byrow = TRUE
    )

    for (k in 1:arms) {
      arm_k_burnin_outcomes <- outcomes[selected_arm[burnin_indices] == k]
      n_k_counts[k] <- n_k_counts[k] + length(arm_k_burnin_outcomes)
      sum_y_k[k] <- sum_y_k[k] + sum(arm_k_burnin_outcomes)
    }

    # Update posterior parameters after burn-in
    for (k in 1:arms) {
      if (n_k_counts[k] > 0) {
        mu0_k <- priors[1, k]
        tau2_0k_orig <- priors[2, k] # This is tau0, need to square for variance
        sigma2_k <- true_vars[k] # Known population variance for arm k

        current_posterior_tau2_params[k] <- 1 / (1 / (tau2_0k_orig^2) + n_k_counts[k] / sigma2_k) # Use tau0^2 here
        current_posterior_mu_params[k] <- current_posterior_tau2_params[k] *
          (mu0_k / (tau2_0k_orig^2) + sum_y_k[k] / sigma2_k) # Use tau0^2 here
      }
    }

    current_range_end <- burnin
  }

  # --- Main Simulation Loop (Block-wise) ---
  start_block_idx <- ifelse(burnin > 0 && burnin <= N, 2, 1)
  if (N == 0) start_block_idx = 1


  for (i in start_block_idx:Nblocks) {
    current_block_size <- block_sizes[i]

    # Handle cases where remaining N is smaller than blocksize
    if (current_range_end + current_block_size > N) {
      current_block_size = N - current_range_end
      if (current_block_size <= 0) break
    }

    current_block_indices <- (current_range_end + 1):(current_range_end + current_block_size)
    batch_number[current_block_indices] <- i

    current_posterior_sds_of_mean <- sqrt(current_posterior_tau2_params)

    # Assumes posterior_norm_sim is available elsewhere in package
    alloc_probs_raw <- posterior_norm_sim(
      means = current_posterior_mu_params,
      sds = current_posterior_sds_of_mean
    )

    # --- Apply tuning parameter (c) ---
    if (tuning == 0) {
      alloc_probs_tuned <- rep(1 / arms, arms)
    } else {
      numerator_vec <- alloc_probs_raw ^ tuning
      denominator_sum <- sum(numerator_vec)
      if (denominator_sum == 0) {
        alloc_probs_tuned <- rep(1 / arms, arms)
      } else {
        alloc_probs_tuned <- numerator_vec / denominator_sum
      }
    }

    # --- Apply clipping ---
    current_clipping_value <- 0
    if (is.numeric(clipping) && clipping > 0) {
      current_clipping_value <- clipping
    } else if (is.character(clipping) && clipping == "adaptive") {
      adaptive_batch_num <- i
      current_clipping_value <- (1 / arms) * (adaptive_batch_num)^(-0.7)
      current_clipping_value <- min(current_clipping_value, 1/arms)
      current_clipping_value <- max(current_clipping_value, 1e-6)
    }

    if (current_clipping_value > 0) {
      lower_bound_per_arm <- current_clipping_value
      upper_bound_per_arm <- 1 - (arms - 1) * current_clipping_value

      alloc_probs_temp <- pmax(alloc_probs_tuned, lower_bound_per_arm)
      alloc_probs_temp <- pmin(alloc_probs_temp, upper_bound_per_arm)

      sum_temp_probs <- sum(alloc_probs_temp)
      if (sum_temp_probs == 0) {
        alloc_probs_final <- rep(1 / arms, arms)
      } else {
        alloc_probs_final <- alloc_probs_temp / sum_temp_probs
      }
    } else {
      alloc_probs_final <- alloc_probs_tuned
    }

    allocation_probs_matrix[current_block_indices, ] <- matrix(
      rep(alloc_probs_final, each = current_block_size),
      ncol = arms, byrow = TRUE
    )

    selected_arm[current_block_indices] <- sample(
      1:arms, current_block_size, prob = alloc_probs_final, replace = TRUE
    )

    # --- Ensure All Arms are Sampled (Exploration Guarantee) ---
    if (ensure_all_arms_sampled && length(unique(selected_arm[current_block_indices])) != arms && current_block_size >= arms) {
      missing_arms <- setdiff(1:arms, unique(selected_arm[current_block_indices]))
      if (length(missing_arms) > 0 && length(missing_arms) <= current_block_size) {
        selected_arm[sample(current_block_indices, length(missing_arms), replace = FALSE)] <- missing_arms
      }
    }

    # --- Simulate Outcomes ---
    outcomes[current_block_indices] <- rnorm(
      current_block_size,
      mean = true_means[selected_arm[current_block_indices]],
      sd = true_sds[selected_arm[current_block_indices]]
    )

    # --- Update Posterior Parameters for the Next Block ---
    if (i < Nblocks) {
      for (k in 1:arms) {
        arm_k_current_block_outcomes <- outcomes[current_block_indices][selected_arm[current_block_indices] == k]
        n_k_new_obs <- length(arm_k_current_block_outcomes)

        if (n_k_new_obs > 0) {
          n_k_counts[k] <- n_k_counts[k] + n_k_new_obs
          sum_y_k[k] <- sum_y_k[k] + sum(arm_k_current_block_outcomes)

          mu0_k <- priors[1, k]
          tau2_0k_orig <- priors[2, k] # This is tau0, not tau0^2
          sigma2_k <- true_vars[k]

          current_posterior_tau2_params[k] <- 1 / (1 / (tau2_0k_orig^2) + n_k_counts[k] / sigma2_k)
          current_posterior_mu_params[k] <- current_posterior_tau2_params[k] *
            (mu0_k / (tau2_0k_orig^2) + sum_y_k[k] / sigma2_k)
        }
      }
    }

    current_range_end <- current_range_end + current_block_size
  }

  return(
    data.frame(
      Batch = batch_number,
      Arm = selected_arm,
      Outcome = outcomes,
      allocation_probs_matrix
    )
  )
}
