#' Simulate a Binary Response Adaptive Randomization (BRAR) Trial (Internal)
#'
#' This is an internal helper function for `simulate_brar_trial`. It simulates a
#' multi-arm clinical trial using Thompson Sampling with a Beta-Bernoulli model
#' for binary outcomes.
#'
#' @param arms Numeric. Number of arms in the trial.
#' @param N Numeric. Total sample size for the trial.
#' @param blocksize Numeric. (Fixed) size of each block of participants.
#' @param priors Matrix. A 2-row matrix where the first row contains alpha parameters
#'   and the second row contains beta parameters for the Beta distributions of each arm.
#' @param modelpar Numeric vector. True success probabilities for each Bernoulli arm.
#' @param tuning Numeric. Tuning parameter for allocation probabilities.
#' @param clipping Numeric. Clipping parameter for allocation probabilities.
#' @param burnin Numeric. Number of initial participants for burn-in.
#' @param ensure_all_arms_sampled Logical. If TRUE, ensures at least one arm is sampled per block.
#' @param postprobmethod Character. Method for calculating posterior probabilities ("simulation" or "exact").
#' @keywords internal
.simulate_brar_trial_binary = function(arms = 2, N, blocksize,
                                       priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
                                       modelpar, tuning = 1, clipping = 0, burnin = 0,
                                       postprobmethod, ensure_all_arms_sampled = FALSE)
{
  # --- Input Validation and Setup ---
  # Only specific validation relevant to this internal function.
  # Broader validation is handled by the main `simulate_brar_trial` function.
  if (length(modelpar) != arms) {
    stop("Length of 'modelpar' must match 'arms'.")
  }
  if (!is.matrix(priors) || nrow(priors) != 2 || ncol(priors) != arms) {
    stop("'priors' must be a 2-row matrix with 'arms' columns (first row: alpha, second row: beta).")
  }
  if (clipping < 0 || clipping >= 1) { # Numeric clipping validation
    stop("The 'clipping' parameter must be between 0 and 1 (exclusive of 1).")
  }
  if (clipping > 0 && (clipping * arms) > 1) {
    stop("Invalid 'clipping' value: clipping * arms must be <= 1 to allow for consistent bounds across all arms and sum of probabilities to be 1.")
  }

  # Determine block sizes for each iteration, considering burn-in
  if (burnin > 0) {
    Nblocks = 1 + floor((N - burnin) / blocksize)
    block_sizes = c(burnin, rep(blocksize, Nblocks - 1))
  } else {
    Nblocks = floor(N / blocksize)
    block_sizes = rep(blocksize, Nblocks)
  }

  # Initialize vectors/matrix to store results for all N participants
  rewards = numeric(N)
  selected_arm = numeric(N)
  batch_number = numeric(N)
  # Matrix to store allocation probabilities for all arms
  allocation_probs_matrix = matrix(NA, nrow = N, ncol = arms)
  colnames(allocation_probs_matrix) = paste0("AlloProb_Arm", 1:arms)

  # --- Simulation for burn-in ---
  if (burnin > 0 && burnin <= N) { # Ensure burnin makes sense for N
    burnin_indices = 1:burnin
    selected_arm[burnin_indices] = blockrand::blockrand(burnin, arms, block.sizes = 1, levels = seq(1, arms, by = 1))$treatment
    rewards[burnin_indices] = rbinom(burnin, 1, modelpar[selected_arm[burnin_indices]])
    batch_number[burnin_indices] = 1

    # Store the allocation probabilities for the current block (equal for burn-in)
    allocation_probs_matrix[burnin_indices, ] = matrix(
      rep(1/arms, each = burnin),
      ncol = arms, byrow = TRUE
    )

    # Initialize Beta distribution parameters (alpha and beta) for each arm
    current_alpha_params = priors[1, ]
    current_beta_params = priors[2, ]

    for(i in 1:arms) {
      current_alpha_params[i] = current_alpha_params[i] + sum(rewards[selected_arm[burnin_indices] == i])
      current_beta_params[i] = current_beta_params[i] + sum(selected_arm[burnin_indices] == i) - sum(rewards[selected_arm[burnin_indices] == i])
    }

    current_range_end = burnin
  } else {
    current_alpha_params = priors[1, ]
    current_beta_params = priors[2, ]
    current_range_end = 0
  }


  # --- Main Simulation Loop (Block-wise) ---
  # Start from block 1 if no burn-in, otherwise from block 2
  start_block_idx <- ifelse(burnin > 0 && burnin <= N, 2, 1) # Ensure we don't start at block 2 if burnin = N
  if (N == 0) start_block_idx = 1 # No blocks if N is 0

  for (i in start_block_idx:Nblocks) {
    current_block_size = block_sizes[i]

    # Handle cases where remaining N is smaller than blocksize
    if (current_range_end + current_block_size > N) {
      current_block_size = N - current_range_end
      if (current_block_size <= 0) break # No more patients to enroll
    }

    current_block_indices = (current_range_end + 1):(current_range_end + current_block_size)
    batch_number[current_block_indices] = i

    # Calculate the raw allocation probabilities for each arm
    # Assumes posterior_bin_sim and posterior_bin_exact are available elsewhere in package
    if(postprobmethod == "simulation") {
      alloc_probs_raw = posterior_bin_sim(alphas = current_alpha_params, betas = current_beta_params)
    } else if(postprobmethod == "exact") {
      alloc_probs_raw = posterior_bin_exact(alphas = current_alpha_params, betas = current_beta_params)
    } else {
      # This case should be caught by main function validation
      stop("Internal Error: Invalid postprobmethod.")
    }

    # --- Apply tuning parameter (c) from Wathen & Thall (2017) ---
    if (tuning == 0) {
      alloc_probs_tuned = rep(1 / arms, arms)
    } else {
      numerator_vec = alloc_probs_raw ^ tuning
      denominator_sum = sum(numerator_vec)
      if (denominator_sum == 0) {
        alloc_probs_tuned = rep(1 / arms, arms) # Fallback to equal if probabilities vanish
      } else {
        alloc_probs_tuned = numerator_vec / denominator_sum
      }
    }

    # --- Apply clipping ---
    if (clipping > 0) {
      lower_bound_per_arm = clipping
      upper_bound_per_arm = 1 - (arms - 1) * clipping

      alloc_probs_temp = pmax(alloc_probs_tuned, lower_bound_per_arm)
      alloc_probs_temp = pmin(alloc_probs_temp, upper_bound_per_arm)

      sum_temp_probs = sum(alloc_probs_temp)
      if (sum_temp_probs == 0) {
        alloc_probs_final = rep(1 / arms, arms)
      } else {
        alloc_probs_final = alloc_probs_temp / sum_temp_probs
      }
    } else {
      alloc_probs_final = alloc_probs_tuned
    }

    # Store the allocation probabilities for the current block
    allocation_probs_matrix[current_block_indices, ] = matrix(
      rep(alloc_probs_final, each = current_block_size),
      ncol = arms, byrow = TRUE
    )

    # Sample arms for the current block
    selected_arm[current_block_indices] = sample(
      1:arms, current_block_size, prob = alloc_probs_final, replace = TRUE
    )

    # --- Ensure All Arms are Sampled (Exploration Guarantee) ---
    if (ensure_all_arms_sampled && length(unique(selected_arm[current_block_indices])) != arms && current_block_size >= arms) {
      missing_arms = setdiff(1:arms, unique(selected_arm[current_block_indices]))
      if (length(missing_arms) > 0 && length(missing_arms) <= current_block_size) {
        selected_arm[sample(current_block_indices, length(missing_arms), replace = FALSE)] = missing_arms
      }
    }

    # --- Simulate Rewards ---
    rewards[current_block_indices] = rbinom(
      current_block_size, size = 1, prob = modelpar[selected_arm[current_block_indices]]
    )

    # --- Update Beta Priors for the Next Block ---
    if (i < Nblocks) {
      for (k in 1:arms) {
        arm_k_indices_in_batch = current_block_indices[selected_arm[current_block_indices] == k]
        successes_arm_k = sum(rewards[arm_k_indices_in_batch])
        failures_arm_k = length(arm_k_indices_in_batch) - successes_arm_k

        current_alpha_params[k] = current_alpha_params[k] + successes_arm_k
        current_beta_params[k] = current_beta_params[k] + failures_arm_k
      }
    }

    current_range_end = current_range_end + current_block_size
  }

  return(
    data.frame(
      Batch = batch_number,
      Arm = selected_arm,
      Outcome = rewards,
      allocation_probs_matrix
    )
  )
}
