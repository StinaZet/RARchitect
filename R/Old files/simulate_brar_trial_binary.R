#' Simulate a Binary Response Adaptive Randomization (BRAR) Trial
#'
#' This function simulates a multi-arm clinical trial using Thompson Sampling
#' with a Beta-Bernoulli model for binary outcomes. It supports fixed block sizes,
#' prior specifications, a tuning parameter for allocation probabilities,
#' clipping of allocation probabilities, burn-in periods, and ensuring at
#' least one arm is sampled per block.
#'
#' @param arms Numeric. Number of arms in the trial.
#' @param N Numeric. Total sample size for the trial.
#' @param blocksize Numeric. (Fixed) size of each block of participants.
#' @param priors Matrix. A 2-row matrix where the first row contains the alpha
#'   parameters and the second row contains the beta parameters for the Beta
#'   distributions of each arm. The number of columns must match `arms`.
#'   e.g., `matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)` for two arms with uniform priors.
#' @param modelpar Numeric vector. A vector with the true success
#'probabilities for each Bernoulli arm. e.g., `c(0.6, 0.4)`.
#' @param tuning Numeric. A parameter (referred to as 'c' in the Wathen and Thall, 2017 paper)
#'that shrinks allocation probabilities towards equal randomization.
#'- If `tuning = 1`, no shrinking occurs; original adaptive randomization probabilities are used.
#'- As `tuning` approaches 0, probabilities shrink towards `1/arms` (equal randomization for `arms` arms).
#'- A common value is 0.5. Defaults to 1.
#' @param clipping Numeric. Forces a minimum and maximum probability of selecting each arm.
#'If `clipping = 0`, no clipping is applied. If `clipping > 0`, each arm's
#'allocation probability is forced to be at least `clipping` and at most
#'`1 - (arms - 1) * clipping`. Probabilities are then re-normalized to sum to 1.
#'It is **required** that `clipping * arms <= 1` for consistent bounds.
#' @param burnin Numeric. The number of initial participants allocated
#'before standard block processing begins. These participants are treated
#'as a single initial block. Defaults to 0.
#' @param ensure_all_arms_sampled Logical. If TRUE, ensures that at least
#'one participant is allocated to each arm within each block (if `blocksize >= arms`).
#'This is useful for ensuring initial exploration. Defaults to FALSE.
#'
#' @return A data frame with `N` rows and the following columns:
#'\itemize{
#'\item `Batch`: The block number for each participant.
#'\item `Arm`: The arm selected for the participant (1 to `arms`).
#'\item `Outcome`: The binary reward (0 or 1) associated with the selected arm.
#'\item `AlloProb_ArmX`: For each arm X, the allocation probability of Arm X for that block.
#'}
#' @export
#'
#' @examples
#'
#' # Example 1: Basic TS with no clipping or burn-in (2 arms)
#' results1 = simulate_brar_trial(arms = 2, N = 100, blocksize = 10, modelpar = c(0.6, 0.4))
#' head(results1)
#' tail(results1)
#'
#' # Example 2: TS with clipping to 0.1 (2 arms)
#' results2 = simulate_brar_trial(arms = 2, N = 100, blocksize = 10, modelpar = c(0.6, 0.4), clipping = 0.1)
#' head(results2)
#'
#' # Example 3: TS with a burn-in period of 20 participants (2 arms)
#' results3 = simulate_brar_trial(arms = 2, N = 100, blocksize = 10, modelpar = c(0.6, 0.4), burnin = 20)
#' head(results3)
#'
#' # Example 4: Ensure all arms are sampled in each block (2 arms)
#' results4 = simulate_brar_trial(arms = 2, N = 100, blocksize = 10, modelpar = c(0.6, 0.4), ensure_all_arms_sampled = TRUE)
#' head(results4)
#'
#' # Example 5: TS with tuning parameter set to 0.5 (more shrinkage towards equal randomization, 2 arms)
#' results5 = simulate_brar_trial(arms = 2, N = 100, blocksize = 10, modelpar = c(0.6, 0.4), tuning = 0.5)
#' head(results5)
#'
#' # Example 6: Using the matrix format for priors (2 arms)
#' priors_matrix = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
#' results6 = simulate_brar_trial(arms = 2, N = 100, blocksize = 10, modelpar = c(0.6, 0.4), priors = priors_matrix)
#' head(results6)
#'
#' # Example 7: Simulation with 3 arms, no tuning, no clipping
#' results7 = simulate_brar_trial(
#'   arms = 3, N = 150, blocksize = 15,
#'   priors = matrix(c(1, 1, 1, 1, 1, 1), nrow = 2, byrow = TRUE),
#'   modelpar = c(0.7, 0.5, 0.3)
#' )
#' head(results7)
#' tail(results7)
#' # Check observed success rates for each arm
#' print(paste("Observed success rate Arm 1:", mean(results7$Outcome[results7$Arm == 1])))
#' print(paste("Observed success rate Arm 2:", mean(results7$Outcome[results7$Arm == 2])))
#' print(paste("Observed success rate Arm 3:", mean(results7$Outcome[results7$Arm == 3])))
#'
#' # Example 8: Simulation with 3 arms, tuning, and clipping (both lower and upper bounds applied)
#' results8 = simulate_brar_trial(
#'   arms = 3, N = 150, blocksize = 15,
#'   priors = matrix(c(1, 1, 1, 1, 1, 1), nrow = 2, byrow = TRUE),
#'   modelpar = c(0.7, 0.5, 0.3),
#'   tuning = 0.5,
#'   clipping = 0.05, # Lower bound for each arm is 0.05
#'                   # Upper bound for each arm is 1 - (3-1)*0.05 = 1 - 0.1 = 0.9
#'   ensure_all_arms_sampled = TRUE
#' )
#' head(results8)
#' tail(results8)
#' # Check average allocation probabilities over the trial
#' apply(results8[, grepl("AlloProb_Arm", names(results8))], 2, mean)
#'
simulate_brar_trial_binary = function(arms = 2, N, blocksize,
                                priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
                                modelpar, tuning = 1, clipping = 0, burnin = 0,
                                postprobmethod, ensure_all_arms_sampled = FALSE)
{

  # --- Input Validation and Setup ---
  if (length(modelpar) != arms) {
    stop("Length of 'modelpar' must match 'arms'.")
  }
  if (!is.matrix(priors) || nrow(priors) != 2 || ncol(priors) != arms) {
    stop("'priors' must be a 2-row matrix with 'arms' columns (first row: alpha, second row: beta).")
  }
  if (tuning < 0) {
    stop("The 'tuning' parameter must be non-negative.")
  }
  if (clipping < 0 || clipping >= 1) {
    stop("The 'clipping' parameter must be between 0 and 1 (exclusive of 1).")
  }
  # Validation for consistent clipping bounds for multiple arms
  if (clipping > 0 && (clipping * arms) > 1) {
    stop("Invalid 'clipping' value: clipping * arms must be <= 1 to allow for consistent bounds across all arms and sum of probabilities to be 1.")
  }
  if (ensure_all_arms_sampled && blocksize < arms && blocksize > 0) {
    stop("If 'ensure_all_arms_sampled' is TRUE, 'blocksize' must be at least 'arms' (or 0 if burnin covers initial).")
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
  selected_arm[1:burnin] = blockrand::blockrand(burnin, arms, block.sizes = 1, levels = seq(1, arms, by = 1))$treatment
  rewards[1:burnin] = rbinom(burnin, 1, modelpar[selected_arm[1:burnin]])
  batch_number[1:burnin] = 1

  # Store the allocation probabilities for the current block
  allocation_probs_matrix[1:burnin, ] = matrix(
    rep(1/arms, each = burnin),
    ncol = arms, byrow = TRUE
  )

  # Initialize Beta distribution parameters (alpha and beta) for each arm
  current_alpha_params = c()
  current_beta_params = c()
  for(i in 1:arms)
  {
    current_alpha_params[i] = priors[1, i] + sum(rewards[selected_arm == i])
    current_beta_params[i] = priors[2, i] + sum(selected_arm == i) - sum(rewards[selected_arm == i])
  }

  # `current_range_end` tracks the index of the last participant allocated
  current_range_end = burnin

  # --- Main Simulation Loop (Block-wise) ---
  for (i in 2:Nblocks) {
    # Define the indices for the current block
    current_block_indices = (current_range_end + 1):(current_range_end + block_sizes[i])
    batch_number[current_block_indices] = i
    current_block_size = block_sizes[i]

    # Calculate the raw allocation probabilities for each arm using posterior_bin_sim.
    if(postprobmethod == "simulation")
    {
      alloc_probs_raw = posterior_bin_sim(alphas = current_alpha_params, betas = current_beta_params)
    }

    # Calculate the raw allocation probabilities for each arm using posterior_bin_exact.
    if(postprobmethod == "exact")
    {
      alloc_probs_raw = posterior_bin_exact(alphas = current_alpha_params, betas = current_beta_params)
    }


    # --- Apply tuning parameter (c) from Wathen & Thall (2017) ---
    # Formula: p_k_tuned = (p_k ^ tuning) / sum(p_j ^ tuning)
    if (tuning == 0) {
      # If tuning is 0, probabilities should become 1/arms (equal randomization)
      alloc_probs_tuned = rep(1 / arms, arms)
    } else {
      numerator_vec = alloc_probs_raw ^ tuning
      denominator_sum = sum(numerator_vec)

      # Avoid division by zero if all probabilities are effectively zero (highly unlikely with Beta priors > 0)
      if (denominator_sum == 0) {
        alloc_probs_tuned = rep(1 / arms, arms) # Fallback to equal if probabilities vanish
      } else {
        alloc_probs_tuned = numerator_vec / denominator_sum
      }
    }

    # --- Apply clipping: Forces a minimum and maximum probability of selecting each arm ---
    if (clipping > 0) {
      lower_bound_per_arm = clipping
      # Upper bound: If all other (arms-1) arms get 'clipping', this arm can get at most 1 - (arms-1)*clipping
      upper_bound_per_arm = 1 - (arms - 1) * clipping

      # Apply lower bound first
      alloc_probs_temp = pmax(alloc_probs_tuned, lower_bound_per_arm)

      # Then apply upper bound
      alloc_probs_temp = pmin(alloc_probs_temp, upper_bound_per_arm)

      # Finally, re-normalize to ensure sum is 1.
      # Note: This renormalization might cause probabilities to slightly fall outside
      # the strict bounds in rare cases if the sum significantly deviates after clipping.
      sum_temp_probs = sum(alloc_probs_temp)
      if (sum_temp_probs == 0) {
        alloc_probs_final = rep(1 / arms, arms) # Fallback if all probs become zero
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

    # Sample arms for the current block based on the calculated 'alloc_probs_final'
    selected_arm[current_block_indices] = sample(
      1:arms, current_block_size, prob = alloc_probs_final, replace = TRUE
    )

    # --- Ensure All Arms are Sampled (Exploration Guarantee) ---
    if (ensure_all_arms_sampled && length(unique(selected_arm[current_block_indices])) != arms && current_block_size >= arms) {
      missing_arms = setdiff(1:arms, unique(selected_arm[current_block_indices]))
      # Randomly reassign a few participants from the current block to the missing arms
      # This prioritizes ensuring all arms get sampled over strict TS allocation for a few participants.
      if (length(missing_arms) <= current_block_size) {
        selected_arm[sample(current_block_indices, length(missing_arms))] = missing_arms
      }
    }

    # --- Simulate Rewards ---
    # Generate binary rewards (0 or 1) based on the true success probabilities
    # of the selected arms.
    rewards[current_block_indices] = rbinom(
      current_block_size, size = 1, prob = modelpar[selected_arm[current_block_indices]]
    )

    # --- Update Beta Priors for the Next Block ---
    # Only update if there are more blocks to come.
    if (i < Nblocks) {
      for (k in 1:arms) {
        arm_k_indices_in_batch = current_block_indices[selected_arm[current_block_indices] == k]
        successes_arm_k = sum(rewards[arm_k_indices_in_batch])
        failures_arm_k = length(arm_k_indices_in_batch) - successes_arm_k

        current_alpha_params[k] = current_alpha_params[k] + successes_arm_k
        current_beta_params[k] = current_beta_params[k] + failures_arm_k
      }
    }

    # Update the end index for the next iteration
    current_range_end = current_range_end + current_block_size
  }

  # Return results as a data frame
  return(
    data.frame(
      Batch = batch_number,
      Arm = selected_arm,
      Outcome = rewards,
      allocation_probs_matrix # This automatically adds columns named AlloProb_ArmX
    )
  )
}
