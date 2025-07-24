#' Simulate a Normal Response Adaptive Randomization (NRAR) Trial
#'
#' This function simulates a multi-arm clinical trial using Thompson Sampling
#' with a Normal-Normal conjugate model for normal outcomes (assuming known
#' population variance for each arm). It supports fixed block sizes, prior
#' specifications, a tuning parameter for allocation probabilities,
#' clipping of allocation probabilities (including adaptive clipping),
#' burn-in periods, and ensuring at least one arm is sampled per block.
#'
#' @param arms Numeric. Number of arms in the trial.
#' @param N Numeric. Total sample size for the trial.
#' @param blocksize Numeric. (Fixed) size of each block of participants.
#' @param priors Matrix. A 2-row matrix containing the Normal prior parameters
#'   for the mean of each arm (assuming known population variance). Each column
#'   corresponds to an arm.
#'   Rows should be:
#'   1. `mu0`: Prior mean for the arm's mean.
#'   2. `tau0`: Prior standard deviation of the arm's mean (often denoted as $\sigma_0$ or $\lambda_0^{-1}$).
#'   Example for two arms with N(0,1) priors for means:
#'   `matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE, dimnames = list(c("mu0", "tau2_0"), paste0("Arm", 1:2)))`
#' @param modelpar Matrix. A 2-row matrix specifying the true parameters
#'   for each normal arm. Each column corresponds to an arm.
#'   Rows should be:
#'   1. `mean`: True mean for the arm.
#'   2. `sd`: True standard deviation for the arm. This `sd` is assumed to be
#'      the known population standard deviation for the purpose of posterior updates.
#'   Example for two arms: `matrix(c(10, 8, 2, 2), nrow = 2, byrow = TRUE,
#'                              dimnames = list(c("mean", "sd"), paste0("Arm", 1:2)))`
#' @param tuning Numeric. A parameter (referred to as 'c' or 'gamma' in the
#'   Wathen and Thall, 2017 paper) that shrinks allocation probabilities
#'   towards equal randomization.
#'   - If `tuning = 1`, no shrinking occurs; original adaptive randomization probabilities are used.
#'   - As `tuning` approaches 0, probabilities shrink towards `1/arms` (equal randomization for `arms` arms).
#'   - A common value is 0.5. Defaults to 1.
#' @param clipping Character or Numeric. Forces a minimum and maximum probability
#'   of selecting each arm.
#'   - If `clipping = 0` (Numeric), no clipping is applied.
#'   - If `clipping > 0` (Numeric), each arm's allocation probability is forced
#'     to be at least `clipping` and at most `1 - (arms - 1) * clipping`.
#'     Probabilities are then re-normalized to sum to 1.
#'     It is **required** that `clipping * arms <= 1` for consistent bounds.
#'   - If `clipping = "adaptive"` (Character), adaptive clipping is applied as
#'     `rho_min,t = 1/arms * (batch_number)^(-0.7)` (as in Hadad et al., 2021).
#' @param burnin Numeric. The number of initial participants allocated
#'   before standard block processing begins. These participants are treated
#'   as a single initial block, allocated using equal randomization. Defaults to 0.
#' @param ensure_all_arms_sampled Logical. If TRUE, ensures that at least
#'   one participant is allocated to each arm within each block (if `blocksize >= arms`).
#'   This is useful for ensuring initial exploration. Defaults to FALSE.
#'
#' @return A data frame with `N` rows and the following columns:
#'   \itemize{
#'     \item `Batch`: The block number for each participant.
#'     \item `Arm`: The arm selected for the participant (1 to `arms`).
#'     \item `Outcome`: The continuous outcome associated with the selected arm.
#'     \item `AlloProb_ArmX`: For each arm X, the allocation probability of Arm X for that block.
#'   }
#' @export
#'
#' @examples
#' # Example 1: Basic TS with normal outcomes (2 arms, known variance)
#' # True parameters: Arm 1 (mean=10, sd=2), Arm 2 (mean=8, sd=2)
#' model_params_norm <- matrix(c(10, 8, 2, 2), nrow = 2, byrow = TRUE,
#'                             dimnames = list(c("mean", "sd"), paste0("Arm", 1:2)))
#'
#' # Prior parameters (Normal-Normal):
#' # For each arm: mu0=0 (prior mean), tau2_0=1 (prior variance of the mean)
#' prior_params_norm <- matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE,
#'                             dimnames = list(c("mu0", "tau2_0"), paste0("Arm", 1:2)))
#'
#' results_norm1 <- simulate_brar_trial_normal(
#'   arms = 2, N = 200, blocksize = 20,
#'   priors = prior_params_norm,
#'   modelpar = model_params_norm,
#'   tuning = 1,
#'   clipping = 0,
#'   burnin = 0,
#'   ensure_all_arms_sampled = TRUE,
#' )
#'
#' head(results_norm1)
#' tail(results_norm1)
#'
#' # Check observed means for each arm
#' print(paste("Observed mean Arm 1:", mean(results_norm1$Outcome[results_norm1$Arm == 1])))
#' print(paste("Observed mean Arm 2:", mean(results_norm1$Outcome[results_norm1$Arm == 2])))
#'
#' # Check average allocation probabilities over the trial
#' apply(results_norm1[, grepl("AlloProb_Arm", names(results_norm1))], 2, mean)
#'
#' # Example 2: Simulation with 3 arms, tuning, and fixed clipping
#' model_params_3arms <- matrix(c(10, 8, 6, 2, 2, 2), nrow = 2, byrow = TRUE,
#'                              dimnames = list(c("mean", "sd"), paste0("Arm", 1:3)))
#' prior_params_3arms <- matrix(c(0, 0, 0, 1, 1, 1), nrow = 2, byrow = TRUE,
#'                              dimnames = list(c("mu0", "tau2_0"), paste0("Arm", 1:3)))
#'
#' results_norm2 <- simulate_brar_trial_normal(
#'   arms = 3, N = 300, blocksize = 30,
#'   priors = prior_params_3arms,
#'   modelpar = model_params_3arms,
#'   tuning = 0.5,
#'   clipping = 0.05, # Lower bound for each arm is 0.05; upper bound is 1 - (3-1)*0.05 = 0.9
#'   burnin = 0,
#'   ensure_all_arms_sampled = TRUE,
#' )
#'
#' head(results_norm2)
#' tail(results_norm2)
#' apply(results_norm2[, grepl("AlloProb_Arm", names(results_norm2))], 2, mean)
#'
#' # Example 3: Simulation with adaptive clipping
#' results_norm3 <- simulate_brar_trial_normal(
#'   arms = 2, N = 200, blocksize = 20,
#'   priors = prior_params_norm,
#'   modelpar = model_params_norm,
#'   tuning = 1,
#'   clipping = "adaptive", # Adaptive clipping applied
#'   burnin = 0,
#'   ensure_all_arms_sampled = FALSE,
#' )
#' head(results_norm3)
#'
simulate_brar_trial_normal <- function(arms = 2, N, blocksize,
                                       priors,
                                       modelpar, tuning = 1, clipping = 0, burnin = 0,
                                       ensure_all_arms_sampled = FALSE)
{

  # --- Input Validation and Setup ---
  if (!is.matrix(modelpar) || nrow(modelpar) != 2 || ncol(modelpar) != arms) {
    stop("'modelpar' must be a 2-row matrix with 'arms' columns (first row: true means, second row: true sds).")
  }
  if (!is.matrix(priors) || nrow(priors) != 2 || ncol(priors) != arms) {
    stop("'priors' must be a 2-row matrix with 'arms' columns (first row: mu0, second row: tau2_0).")
  }
  if (tuning < 0) {
    stop("The 'tuning' parameter must be non-negative.")
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
  if (ensure_all_arms_sampled && blocksize < arms && blocksize > 0) {
    stop("If 'ensure_all_arms_sampled' is TRUE, 'blocksize' must be at least 'arms' (or 0 if burnin covers initial).")
  }


  # Extract true parameters from modelpar for convenience
  true_means <- modelpar[1, ]
  true_sds <- modelpar[2, ]
  true_vars <- true_sds^2 # Convert true sds to true variances

  # Determine block sizes for each iteration, considering burn-in
  if (burnin > 0) {
    Nblocks <- 1 + floor((N - burnin) / blocksize)
    # The first 'block' is the burn-in period
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
  # These will be updated after each block.
  # current_mu_n (posterior mean of the mean)
  # current_tau2_n (posterior sd of the mean)
  current_posterior_mu_params <- priors[1, ] # mu0
  current_posterior_tau2_params <- priors[2, ]^2 # tau0^2

  # Count of observations for each arm (n_k)
  n_k_counts <- rep(0, arms)
  # Sum of outcomes for each arm (sum_y_k)
  sum_y_k <- rep(0, arms)

  # `current_range_end` tracks the index of the last participant allocated
  current_range_end <- 0

  # --- Simulation for burn-in period ---
  if (burnin > 0) {
    burnin_indices <- 1:burnin
    batch_number[burnin_indices] <- 1

    # For burn-in, allocate equally (using sampling to fill the exact count)
    selected_arm[burnin_indices] <- blockrand::blockrand(burnin, arms, block.sizes = 1, levels = seq(1, arms, by = 1))$treatment

    # Simulate outcomes for burn-in participants
    outcomes[burnin_indices] <- rnorm(burnin,
                                      mean = true_means[selected_arm[burnin_indices]],
                                      sd = true_sds[selected_arm[burnin_indices]])

    # Store equal allocation probabilities for burn-in period
    allocation_probs_matrix[burnin_indices, ] <- matrix(
      rep(1 / arms, each = burnin),
      ncol = arms, byrow = TRUE
    )

    # Update sufficient statistics based on burn-in data
    for (k in 1:arms) {
      arm_k_burnin_outcomes <- outcomes[selected_arm[burnin_indices] == k]
      n_k_counts[k] <- n_k_counts[k] + length(arm_k_burnin_outcomes)
      sum_y_k[k] <- sum_y_k[k] + sum(arm_k_burnin_outcomes)
    }

    # Update posterior parameters after burn-in
    for (k in 1:arms) {
      if (n_k_counts[k] > 0) {
        mu0_k <- priors[1, k]
        tau2_0k <- priors[2, k]
        sigma2_k <- true_vars[k] # Known population variance for arm k

        # Posterior variance of the mean (tau2_nk)
        current_posterior_tau2_params[k] <- 1 / (1 / tau2_0k + n_k_counts[k] / sigma2_k)
        # Posterior mean of the mean (mu_nk)
        current_posterior_mu_params[k] <- current_posterior_tau2_params[k] *
          (mu0_k / tau2_0k + sum_y_k[k] / sigma2_k)
      }
    }

    current_range_end <- burnin
  }

  # --- Main Simulation Loop (Block-wise) ---
  # Start from block 1 if no burn-in, otherwise from block 2
  start_block_idx <- ifelse(burnin > 0, 2, 1)

  for (i in start_block_idx:Nblocks) {
    # Define the indices for the current block
    current_block_indices <- (current_range_end + 1):(current_range_end + block_sizes[i])
    batch_number[current_block_indices] <- i
    current_block_size <- block_sizes[i]

    # Calculate current posterior standard deviations of the mean for `posterior_norm_sim`
    current_posterior_sds_of_mean <- sqrt(current_posterior_tau2_params)

    # Calculate the raw allocation probabilities for each arm
    alloc_probs_raw <- posterior_norm_sim(
      means = current_posterior_mu_params,
      sds = current_posterior_sds_of_mean
    )

    # --- Apply tuning parameter (c) from Wathen & Thall (2017) ---
    if (tuning == 0) {
      alloc_probs_tuned <- rep(1 / arms, arms)
    } else {
      numerator_vec <- alloc_probs_raw ^ tuning
      denominator_sum <- sum(numerator_vec)

      if (denominator_sum == 0) { # Fallback if all probabilities vanish (should be rare)
        alloc_probs_tuned <- rep(1 / arms, arms)
      } else {
        alloc_probs_tuned <- numerator_vec / denominator_sum
      }
    }

    # --- Apply clipping: Forces a minimum and maximum probability of selecting each arm ---
    current_clipping_value <- 0 # Default if no clipping or adaptive clipping not applied yet

    if (is.numeric(clipping) && clipping > 0) {
      current_clipping_value <- clipping
    } else if (is.character(clipping) && clipping == "adaptive") {
      # Use the batch number for adaptive clipping, adjusted if burnin shifted batch numbers
      adaptive_batch_num <- i # The batch index 'i' corresponds to the block in sequence
      current_clipping_value <- (1 / arms) * (adaptive_batch_num)^(-0.7)
      # Ensure adaptive clipping doesn't exceed 1/arms or fall too low
      current_clipping_value <- min(current_clipping_value, 1/arms) # Max value is 1/arms (initial equal allocation)
      current_clipping_value <- max(current_clipping_value, 1e-6) # Keep a small minimum to avoid zero
    }

    if (current_clipping_value > 0) {
      lower_bound_per_arm <- current_clipping_value
      upper_bound_per_arm <- 1 - (arms - 1) * current_clipping_value

      # Apply lower bound first
      alloc_probs_temp <- pmax(alloc_probs_tuned, lower_bound_per_arm)

      # Then apply upper bound
      alloc_probs_temp <- pmin(alloc_probs_temp, upper_bound_per_arm)

      # Finally, re-normalize to ensure sum is 1.
      sum_temp_probs <- sum(alloc_probs_temp)
      if (sum_temp_probs == 0) { # Fallback if all probs become zero after clipping
        alloc_probs_final <- rep(1 / arms, arms)
      } else {
        alloc_probs_final <- alloc_probs_temp / sum_temp_probs
      }
    } else {
      alloc_probs_final <- alloc_probs_tuned
    }

    # Store the allocation probabilities for the current block
    allocation_probs_matrix[current_block_indices, ] <- matrix(
      rep(alloc_probs_final, each = current_block_size),
      ncol = arms, byrow = TRUE
    )

    # Sample arms for the current block based on the calculated 'alloc_probs_final'
    selected_arm[current_block_indices] <- sample(
      1:arms, current_block_size, prob = alloc_probs_final, replace = TRUE
    )

    # --- Ensure All Arms are Sampled (Exploration Guarantee) ---
    if (ensure_all_arms_sampled && length(unique(selected_arm[current_block_indices])) != arms && current_block_size >= arms) {
      missing_arms <- setdiff(1:arms, unique(selected_arm[current_block_indices]))
      if (length(missing_arms) > 0 && length(missing_arms) <= current_block_size) {
        # Randomly reassign participants to ensure missing arms are sampled
        selected_arm[sample(current_block_indices, length(missing_arms), replace = FALSE)] <- missing_arms
      }
    }
    # Error handling for blocksize == 1 with ensure_all_arms_sampled
    if (ensure_all_arms_sampled && blocksize == 1 && arms > 1) {
      stop("You have a block size of 1 and more than one arm, but you asked to ensure each arm is sampled at least once in each block! This is impossible. Please change the block size or set 'ensure_all_arms_sampled = FALSE'.")
    }

    # --- Simulate Outcomes ---
    # Generate normal outcomes based on the true means and sds of the selected arms.
    outcomes[current_block_indices] <- rnorm(
      current_block_size,
      mean = true_means[selected_arm[current_block_indices]],
      sd = true_sds[selected_arm[current_block_indices]]
    )

    # --- Update Posterior Parameters for the Next Block ---
    # Only update if there are more blocks to come.
    if (i < Nblocks) {
      for (k in 1:arms) {
        arm_k_current_block_outcomes <- outcomes[current_block_indices][selected_arm[current_block_indices] == k]
        n_k_new_obs <- length(arm_k_current_block_outcomes)

        if (n_k_new_obs > 0) {
          # Update sufficient statistics
          n_k_counts[k] <- n_k_counts[k] + n_k_new_obs
          sum_y_k[k] <- sum_y_k[k] + sum(arm_k_current_block_outcomes)

          mu0_k <- priors[1, k]
          tau2_0k <- priors[2, k]
          sigma2_k <- true_vars[k] # Known population variance for arm k

          # Update posterior parameters using Normal-Normal conjugate formulas
          current_posterior_tau2_params[k] <- 1 / (1 / tau2_0k + n_k_counts[k] / sigma2_k)
          current_posterior_mu_params[k] <- current_posterior_tau2_params[k] *
            (mu0_k / tau2_0k + sum_y_k[k] / sigma2_k)
        }
      }
    }

    # Update the end index for the next iteration
    current_range_end <- current_range_end + current_block_size
  }

  # Return results as a data frame
  return(
    data.frame(
      Batch = batch_number,
      Arm = selected_arm,
      Outcome = outcomes,
      allocation_probs_matrix # This automatically adds columns named AlloProb_ArmX
    )
  )
}
