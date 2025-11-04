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
#' @param multiarm_method Character. Method for handling >2 arms. Either `"top2"`
#' for Top 2 Thompson Sampling, or `"fixed"` for a fixed ratio to the control arm.
#' @keywords internal
.simulate_brar_trial_binary <- function(direction, arms, N, blocksize, priors, modelpar,
                                       tuning = 1, clipping = 0, burnin = 0,
                                       postprobmethod, randmethod = "coin",
                                       urn_alpha, multiarm_method)
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
  start_block_idx = ifelse(burnin > 0 && burnin <= N, 2, 1) # Ensure we don't start at block 2 if burnin = N
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
    if(postprobmethod == "simulation") {
      alloc_probs_raw = posterior_bin_sim(alphas = current_alpha_params, betas = current_beta_params, direction = direction)
    } else if(postprobmethod == "exact") {
      alloc_probs_raw = posterior_bin_exact(alphas = current_alpha_params, betas = current_beta_params, direction = direction)
    } else {
      # This case should be caught by main function validation
      stop("Internal Error: Invalid postprobmethod.")
    }

    # Other methods for allocation probabilities if there are more than 2 arms.
    if (arms > 2)
    {
      if (multiarm_method == "fixed")
      {
        if (alloc_probs_raw[1] >= 1 / arms)
        {
          # Don't change anything if the allocation probability to the control
          # arm is larger than or equal to 1 / arms.
          alloc_probs_raw = alloc_probs_raw
        } else{
          # The sum of the other probabilities
          remaining_mass = 1 - 1 / arms

          # Normalize the other arms (2:K) to sum to 1
          other_probs = alloc_probs_raw[-1]
          other_probs = other_probs / sum(other_probs)

          # Rescale them to fit into the remaining mass
          other_probs = other_probs * remaining_mass

          # Replace in the original vector
          alloc_probs_raw = c(1 / arms, other_probs)
        }
      } else if (multiarm_method == "top2"){
        # The beta parameter for Top 2 Thompson Sampling.
        top2beta = 0.5
        alloc_probs_raw = alloc_probs_T2TS(alloc_probs_raw, top2beta)
      } else {
        # This case should be caught by main function validation
        stop("Internal Error: Invalid multiarm_method.")
      }
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
    current_clipping_value = 0
    if (is.numeric(clipping) && clipping > 0) {
      current_clipping_value = clipping
    } else if (is.character(clipping) && clipping == "adaptive") {
      adaptive_batch_num = i
      current_clipping_value = (1 / arms) * (adaptive_batch_num)^(-0.7)
      current_clipping_value = min(current_clipping_value, 1/arms)
      current_clipping_value = max(current_clipping_value, 1e-6)
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


    # Determine the current block size
    current_block_size = length(current_block_indices)

    if(randmethod == "block")
    {
      # https://journals.sagepub.com/doi/full/10.1177/17407745211010139
      # Multi-Arm Block Randomization (BAR compliant: Fixed counts for the block)

      # Store the allocation probabilities for the current block.
      allocation_probs_matrix[current_block_indices, ] = matrix(
        rep(alloc_probs_final, each = current_block_size),
        ncol = arms, byrow = FALSE
      )

      # --- 1. Calculate integer assignment counts for all K arms (Quota Sampling) ---
      raw_counts = alloc_probs_final * current_block_size
      base_counts = floor(raw_counts)
      remainder = current_block_size - sum(base_counts)

      fractional_parts = raw_counts - base_counts

      if (remainder > 0) {
        if (all(fractional_parts == 0)) fractional_parts = rep(1/arms, arms)

        # Distribute the 'remainder' slots probabilistically
        add_indices = sample(
          1:arms,
          size = remainder,
          prob = fractional_parts,
          replace = FALSE
        )
        add_tab = tabulate(add_indices, nbins = arms)
        counts_final = base_counts + add_tab
      } else {
        counts_final = base_counts
      }

      # --- 2. Construct the assignment vector and shuffle ---
      arm_assignments = rep(1:arms, times = counts_final)

      # Shuffle assignments within block to prevent predictability
      selected_arm[current_block_indices] = sample(arm_assignments, size = current_block_size, replace = FALSE)

    } else if(randmethod == "urn"){
      # https://www.sciencedirect.com/science/article/pii/S1551714415300264?via%3Dihub
      # Multi-Arm Design-Adaptive Urn Randomization (DAR - Outcome-Blind, Zhao 2015 generalized)

      # Urn Mass parameter (alpha).
      alpha_urn = urn_alpha

      # Target probabilities for the block (pi_k)
      target_probs = alloc_probs_final

      # Initialize tracking within the current block: N_k(i-1)
      # This tracks the number of patients *assigned* to each arm in the block so far.
      count_in_block = rep(0L, arms)

      selected_arm_block = numeric(current_block_size)
      rewards_block = numeric(current_block_size)

      # Matrix to store the calculated allocation probabilities for each draw in the block
      probs_per_draw = matrix(NA, nrow = current_block_size, ncol = arms)

      # --- Patient-by-Patient Urn Process within the block ---
      for (iii in 1:current_block_size)
      {
        # 1. Calculate Urn Weights W_k(i) based on Zhao (2015) Eq 7a:
        # W_k(i) = max(alpha * pi_k + (i-1) * pi_k - N_k(i-1), 1)

        # Term 1: alpha * pi_k + (i-1) * pi_k
        term1 = alpha_urn * target_probs + (iii - 1) * target_probs

        # Weight W_k: Weights are adjusted down based on number of patients assigned (count_in_block).
        weights = pmax(term1 - count_in_block, 1)

        # 2. Calculate current allocation probabilities
        probs_now = weights / sum(weights)
        probs_per_draw[iii, ] = probs_now

        # 3. Draw arm
        draw = sample(1:arms, size = 1, prob = probs_now)
        selected_arm_block[iii] = draw

        # 4. Simulate the outcome (NOTE: The outcome is NOT used for the *next* draw)
        outcome = stats::rbinom(1, 1, modelpar[draw])
        rewards_block[iii] = outcome

        # 5. Update assignment counts for the next patient (i+1)
        count_in_block[draw] = count_in_block[draw] + 1L # Increment N_k for the drawn arm
      }

      # --- 6. Save results ---
      # Store the probabilities used for randomization (probs_per_draw)
      allocation_probs_matrix[current_block_indices, ] = probs_per_draw

      selected_arm[current_block_indices] = selected_arm_block

      rewards[current_block_indices] = rewards_block

    } else { # This is the "coin" or default Thompson Sampling allocation
      # Coin Randomization (BAR compliant: Simple random sampling/Multinomial)

      # Store the allocation probabilities for the current block
      allocation_probs_matrix[current_block_indices, ] = matrix(
        rep(alloc_probs_final, each = current_block_size),
        ncol = arms, byrow = FALSE
      )

      # The selected arms.
      selected_arm[current_block_indices] = sample(
        1:arms, current_block_size, prob = alloc_probs_final, replace = TRUE
      )
    }

    # --- Simulate Rewards (Only needed for "block" and "coin" as "urn" already simulated them) ---
    if (randmethod != "urn") {
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
