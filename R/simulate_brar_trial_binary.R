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
#' @param randmethod Character. The randomisation method when blocksize > 1 and
#' arms = 2. Defaults to "coin".
#' @param postprobmethod Character. Method for calculating posterior probabilities ("simulation" or "exact").
#' @keywords internal
.simulate_brar_trial_binary = function(arms = 2, N, blocksize, priors, modelpar,
                                       tuning = 1, clipping = 0, burnin = 0,
                                       postprobmethod, randmethod = "coin")
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

    # Shuffle to avoid any ordering bias
    selected_arm[burnin_idx] = sample(arm_assignments, burnin)

    # Simulate rewards for burn-in participants
    rewards[burnin_idx] = stats::rbinom(burnin, 1, modelpar[selected_arm[burnin_idx]])
    batch_number[burnin_idx] = 1

    # Store the allocation probabilities for burn-in (equal)
    allocation_probs_matrix[burnin_idx, ] = matrix(1/arms, nrow = burnin, ncol = arms)

    # Initialize Beta distribution parameters (alpha and beta) for each arm
    current_alpha_params = priors[1, ]
    current_beta_params = priors[2, ]

    for(i in 1:arms) {
      current_alpha_params[i] = current_alpha_params[i] + sum(rewards[selected_arm[burnin_idx] == i])
      current_beta_params[i] = current_beta_params[i] + sum(selected_arm[burnin_idx] == i) - sum(rewards[selected_arm[burnin_idx] == i])
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
      lower_bound_per_arm = current_clipping_value
      upper_bound_per_arm = 1 - (arms - 1) * current_clipping_value

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

    alloc_probs_final = round(alloc_probs_final, digits = 10)


    if(randmethod == "block")
    {
      # From Proper, Connett, and Murray (2021). https://journals.sagepub.com/doi/full/10.1177/17407745211010139
      # Store the allocation probabilities for the current block. These are
      # the same as what they would be with the coin design according to
      # the original work.
      allocation_probs_matrix[current_block_indices, ] = matrix(
        rep(alloc_probs_final, each = current_block_size),
        ncol = arms, byrow = FALSE
      )

      # The target allocation ratio.
      target = alloc_probs_final[1] * blocksize
      # The floor, defined as target - 1 if target is an integer.
      below = ifelse(target %% 1 == 0, target - 1, floor(target))
      # The ceiling of target.
      above = ceiling(target)

      # Randomise if the floor or ceiling is used.
      u = stats::rbinom(1, 1, (target - below))
      e = u * above + (1 - u) * floor

      # The number of patients on each arm.
      arm_assignments = c(rep(1, times = e), rep(2, times = blocksize - e))

      # Shuffle to avoid any ordering bias
      selected_arm[current_block_indices] = sample(arm_assignments, blocksize)

      # --- Simulate Rewards ---
      rewards[current_block_indices] = stats::rbinom(
        current_block_size, size = 1, prob = modelpar[selected_arm[current_block_indices]]
      )

    } else if(randmethod == "urn"){
      # From Zhao (2015). https://www.sciencedirect.com/science/article/pii/S1551714415300264?via%3Dihub

      # The first values of the probabilities are the usual probabilities.
      urnprob = alloc_probs_final[1]

      # The alpha value for the urn-design.
      alpha = 3
      for (iii in 1:blocksize)
      {
        # Simulate the treatment and outcome
        treatment[iii] = 1 + stats::rbinom(1, 1, urnprob[iii])
        outcome[iii] = stats::rbinom(1, 1, modelpar[treatment[iii]])

        # Update the allocation probabilities.
        term1 = max(alpha * alloc_probs_final[1] - sum(outcome) + (iii - 1) * alloc_probs_final[1], 0)
        term2 = max(alpha * (1 - alloc_probs_final[1]) - (length(outcome) - sum(outcome)) + (iii - 1) * (1 - alloc_probs_final[1]), 0)
        urnprob[iii + 1] = term1 / (term1 + term2)
      }

      # Save in the relevant matrices for output.
      # Store the allocation probabilities for the current block
      allocation_probs_matrix[current_block_indices, ] = matrix(
        c(urnprob, (1 - urnprob)),
        ncol = arms, byrow = FALSE
      )

      # The selected arms.
      selected_arm[current_block_indices] = treatment

      rewards[current_block_indices] = outcome

    } else{

      # Store the allocation probabilities for the current block
      allocation_probs_matrix[current_block_indices, ] = matrix(
        rep(alloc_probs_final, each = current_block_size),
        ncol = arms, byrow = FALSE
      )

      # The selected arms.
      selected_arm[current_block_indices] = sample(
        1:arms, current_block_size, prob = alloc_probs_final, replace = TRUE
      )

      # --- Simulate Rewards ---
      rewards[current_block_indices] = stats::rbinom(
        current_block_size, size = 1, prob = modelpar[selected_arm[current_block_indices]]
      )
    }


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
