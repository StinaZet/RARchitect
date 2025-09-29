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
.simulate_brar_trial_normal <- function(arms = 2, N, blocksize, priors, modelpar,
                                        tuning = 1, clipping = 0, burnin = 0,
                                        postprobmethod, ensure_all_arms_sampled = FALSE,
                                        known_var = FALSE)
{
  # --- Input Validation and Setup ---
  # Only specific validation relevant to this internal function.
  # Broader validation is handled by the main `simulate_brar_trial` function.
  if (!is.matrix(modelpar) || nrow(modelpar) != 2 || ncol(modelpar) != arms) {
    stop("'modelpar' must be a 2-row matrix with 'arms' columns (first row: true means, second row: true sds).")
  }
  if (known_var) {
    if (!is.matrix(priors) || nrow(priors) != 2 || ncol(priors) != arms) {
      stop("'priors' must be a 2-row matrix (mu0, tau0) for known-variance case.")
    }
  } else {
    if (!is.matrix(priors) || nrow(priors) != 4 || ncol(priors) != arms) {
      stop("'priors' must be a 4-row matrix (mu0, kappa0, alpha0, beta0) for unknown-variance case.")
    }
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

  # --- Initialize Priors ---
  if (known_var) {
    # Normal-Normal
    current_mu <- priors[1, ]
    current_tau2 <- priors[2, ]^2
  } else {
    # Normal-Inverse-Gamma
    current_mu    <- priors[1, ]
    current_kappa <- priors[2, ]
    current_alpha <- priors[3, ]
    current_beta  <- priors[4, ]
  }

  # Count of observations for each arm (n_k)
  n_k <- rep(0, arms)
  # Sum of outcomes for each arm (sum_y_k)
  sum_y_k <- rep(0, arms)
  # Sum of squared outcomes for each arm (sum_y_k)
  sum_y2_k <- rep(0, arms)

  current_range_end <- 0

  # --- Simulation for burn-in period ---
  # Initialize vectors/matrix to store results for all N participants
  rewards = numeric(N)           # continuous outcomes now
  selected_arm = numeric(N)
  batch_number = numeric(N)
  allocation_probs_matrix = matrix(NA, nrow = N, ncol = arms)
  colnames(allocation_probs_matrix) = paste0("AlloProb_Arm", 1:arms)

  # --- Simulation for burn-in ---
  if (burnin > 0 && burnin <= N) {
    burnin_idx = 1:burnin

    # Generate roughly balanced allocation
    full_cycles = floor(burnin / arms)       # number of full cycles
    remainder = burnin %% arms               # leftover participants

    # Repeat each arm for full cycles
    arm_assignments = rep(1:arms, times = full_cycles)

    # Add remaining participants randomly among the arms
    if (remainder > 0) {
      arm_assignments = c(arm_assignments, sample(1:arms, remainder))
    }

    # Shuffle to avoid ordering bias
    selected_arm[burnin_idx] = sample(arm_assignments, burnin)

    # Simulate outcomes for burn-in participants
    outcomes[burnin_idx] <- stats::rnorm(burnin,
                                         mean = true_means[selected_arm[burnin_idx]],
                                         sd = true_sds[selected_arm[burnin_idx]])


    batch_number[burnin_idx] = 1

    # Store allocation probabilities (equal for burn-in)
    allocation_probs_matrix[burnin_idx, ] = matrix(
      rep(1/arms, each = burnin),
      ncol = arms, byrow = TRUE
    )

    for (k in 1:arms) {
      yk <- outcomes[burnin_idx][selected_arm[burnin_idx] == k]
      n_k[k] <- length(yk)
      sum_y_k[k] <- sum(yk)
      sum_y2_k[k] <- sum(yk^2)
    }

    # Update posterior parameters after burn-in
    for (k in 1:arms) {
      if (n_k[k] > 0) {
        if (known_var) {
          # Known variance case
          mu0 <- priors[1, k] # Prior mean
          tau0 <- priors[2, k] # Prior sd
          sigma2 <- true_vars[k] # True variance
          current_tau2[k] <- 1 / (1 / (tau0^2) + n_k[k] / sigma2) # Posterior sd
          current_mu[k]   <- current_tau2[k] * (mu0 / (tau0^2) + sum_y_k[k] / sigma2) # Posterior mean
        } else {
          # Unknown variance case
          # Prior values
          mu0 <- priors[1, k]
          kappa0 <- priors[2, k]
          alpha0 <- priors[3, k]
          beta0 <- priors[4, k]
          # Data from the burn-in
          ybar <- sum_y_k[k] / n_k[k]
          S <- sum_y2_k[k] - n_k[k] * ybar^2
          # Posterior values
          current_kappa[k] <- kappa0 + n_k[k]
          current_mu[k]    <- (kappa0 * mu0 + n_k[k] * ybar) / current_kappa[k]
          current_alpha[k] <- alpha0 + n_k[k] / 2
          current_beta[k]  <- beta0 + 0.5 * S +
            (kappa0 * n_k[k]) / (2 * current_kappa[k]) * (ybar - mu0)^2
        }
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

    idx <- (current_range_end + 1):(current_range_end + current_block_size)
    batch_number[idx] <- i


    # Calculate the raw allocation probabilities for each arm
    if (known_var) {
      post_sds <- sqrt(current_tau2)
      if (postprobmethod == "simulation") {
        alloc_probs_raw <- posterior_norm_sim(current_mu, post_sds)
      } else if(postprobmethod == "exact") {
        alloc_probs_raw <- posterior_norm_exact(current_mu, post_sds)
      } else {
        # This case should be caught by main function validation
        stop("Internal Error: Invalid postprobmethod.")
      }
    } else {
      alloc_probs_raw <- posterior_norm_unknownvar_sim(
        mu_n = current_mu,
        kappa_n = current_kappa,
        alpha_n = current_alpha,
        beta_n = current_beta
      )
    }

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

    alloc_probs_final = round(alloc_probs_final, digits = 10)

    allocation_probs_matrix[idx, ] <- matrix(
      rep(alloc_probs_final, each = current_block_size),
      ncol = arms, byrow = FALSE
    )

    selected_arm[idx] <- sample(
      1:arms, current_block_size, prob = alloc_probs_final, replace = TRUE
    )

    # --- Ensure All Arms are Sampled (Exploration Guarantee) ---
    if (ensure_all_arms_sampled && length(unique(selected_arm[idx])) != arms && current_block_size >= arms) {
      missing_arms <- setdiff(1:arms, unique(selected_arm[idx]))
      if (length(missing_arms) > 0 && length(missing_arms) <= current_block_size) {
        selected_arm[sample(idx, length(missing_arms), replace = FALSE)] <- missing_arms
      }
    }

    # --- Simulate Outcomes ---
    outcomes[idx] <- stats::rnorm(
      current_block_size,
      mean = true_means[selected_arm[idx]],
      sd = true_sds[selected_arm[idx]]
    )

    # --- Update Posterior Parameters for the Next Block ---
    if (i < Nblocks) {
        for (k in 1:arms) {
          yk <- outcomes[idx][selected_arm[idx] == k]
          if (length(yk) > 0) {
            n_k[k] <- n_k[k] + length(yk)
            sum_y_k[k] <- sum_y_k[k] + sum(yk)
            sum_y2_k[k] <- sum_y2_k[k] + sum(yk^2)
            if (known_var) {
              mu0 <- priors[1, k]
              tau0 <- priors[2, k]
              sigma2 <- true_vars[k]
              current_tau2[k] <- 1 / (1 / (tau0^2) + n_k[k] / sigma2)
              current_mu[k]   <- current_tau2[k] * (mu0 / (tau0^2) + sum_y_k[k] / sigma2)
            } else {
              mu0 <- priors[1, k]
              kappa0 <- priors[2, k]
              alpha0 <- priors[3, k]
              beta0 <- priors[4, k]
              ybar <- sum_y_k[k] / n_k[k]
              S <- sum_y2_k[k] - n_k[k] * ybar^2
              current_kappa[k] <- kappa0 + n_k[k]
              current_mu[k]    <- (kappa0 * mu0 + n_k[k] * ybar) / current_kappa[k]
              current_alpha[k] <- alpha0 + n_k[k] / 2
              current_beta[k]  <- beta0 + 0.5 * S +
                                  (kappa0 * n_k[k]) / (2 * current_kappa[k]) * (ybar - mu0)^2
            }
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
