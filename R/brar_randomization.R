#' Generate Randomization Sequences for a BRAR Trial
#'
#' Generates randomization sequences or allocation probabilities for a BRAR trial
#' given observed outcomes, priors, and allocation settings.
#'
#' @param trial_data Data frame with columns 'Arm', 'Outcome', and 'Batch'.
#' @param priors 2 x K matrix of alpha/beta parameters for each arm.
#' @param blocksize Numeric. Block size for the next block in the BRAR trial.
#' @param direction Character. "higher" or "lower".
#' @param postprobmethod Character. "simulation" or "exact".
#' @param multiarm_method Character. "top2" or "fixed".
#' @param tuning Numeric. Tuning exponent for allocation probabilities.
#' @param clipping Numeric or "adaptive". Minimum allocation probability per arm.
#' @param urn_alpha Numeric. Urn alpha parameter (if randmethod = "urn").
#' @param randmethod Character. "coin", "block", or "urn".
#' @param return Character. Either "allocations" (default) or "probabilities".
#' @return Either a vector of randomized allocations or a matrix of allocation probabilities.
#' @export
brar_randomization <- function(
    trial_data, priors, blocksize, direction = c("higher", "lower"),
    postprobmethod = c("simulation", "exact"), multiarm_method = c("top2", "fixed"),
    tuning = 1, clipping = 0, urn_alpha = NULL, randmethod = c("coin", "block", "urn"),
    return = c("allocations", "probabilities")) {


  direction = match.arg(direction)
  postprobmethod = match.arg(postprobmethod)
  multiarm_method = match.arg(multiarm_method)
  randmethod = match.arg(randmethod)
  return = match.arg(return)

  if (randmethod == "urn" && (urn_alpha <= 0 || !is.numeric(urn_alpha))) {
    stop("The 'urn_alpha' parameter must be a positive numeric value when randmethod='urn'.")
  }

  N = nrow(trial_data)
  arms = length(unique(trial_data$Arm))
  b = max(trial_data$Batch) + 1 # Used for adaptive clipping.

  allocations = numeric(N)
  allocation_probs_matrix = matrix(NA, nrow = 1, ncol = arms)

  # Initialize alpha/beta with priors + observed outcomes
  current_alpha = priors[1, ]
  current_beta = priors[2, ]
  for (k in 1:arms)
  {
    y_k = sum(trial_data$Outcome[trial_data$Arm == k])
    n_k = sum(trial_data$Arm == k)
    current_alpha[k] = current_alpha[k] + y_k
    current_beta[k] = current_beta[k] + n_k - y_k
  }

  # Posterior probabilities
  if (postprobmethod == "simulation") {
    alloc_probs = posterior_bin_sim(current_alpha, current_beta, direction)
  } else {
    alloc_probs = posterior_bin_exact(current_alpha, current_beta, direction)
  }

  # Multi-arm adjustments
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

  # Tuning
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



  # Clipping
  if (is.character(clipping) && clipping == "adaptive") {
    clipping_val = min(1/arms, (1/arms) * b^(-0.7))
  } else {
    clipping_val = max(0, clipping)
  }

  if (clipping_val > 0) {
    lower_bound_per_arm = clipping_val
    upper_bound_per_arm = 1 - (arms - 1) * clipping_val

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


  # Randomization
  if (randmethod == "coin") {
    allocations = sample(1:arms, blocksize, prob = alloc_probs_tuned, replace = TRUE)
  } else if (randmethod == "block") {
    counts = floor(alloc_probs_tuned * blocksize)
    remainder = blocksize - sum(counts)
    if (remainder > 0) {
      frac_parts = alloc_probs_tuned - counts / blocksize
      add_idx = sample(1:arms, remainder, prob = frac_parts, replace = FALSE)
      counts = counts + tabulate(add_idx, nbins = arms)
    }
    alloc_vector = rep(1:arms, counts)
    allocations = sample(alloc_vector, blocksize)
  } else if (randmethod == "urn") {
    count_in_block = rep(0, arms)
    for (i in 1:blocksize) {
      weights = pmax(urn_alpha * alloc_probs_tuned + (i-1) * alloc_probs_tuned - count_in_block, 1)
      prob_now = weights / sum(weights)
      draw = sample(1:arms, 1, prob = prob_now)
      allocations[i] = draw
      count_in_block[draw] = count_in_block[draw] + 1
    }
  }

  # Return the allocation probability for the first individual in the block.
  # The following will depend on the allocation of the previous, so it would
  # have to be different scenarios if all possibilities would be covered.
  allocation_probs_matrix = matrix(alloc_probs_tuned, ncol = arms, byrow = FALSE)

  if (return == "allocations") {
    return(allocations)
  } else {
    return(allocation_probs_matrix)
  }
}
