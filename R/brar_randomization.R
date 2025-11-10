#' Generate Randomization Sequences for a BRAR Trial (Full)
#'
#' Generates randomization sequences or allocation probabilities for a BRAR trial
#' given observed outcomes, priors, and allocation settings.
#'
#' @param trial_data Data frame with columns 'Arm', 'Outcome', 'Batch'.
#' @param priors 2 x K matrix of alpha/beta parameters for each arm.
#' @param blocksize Numeric. Block size used in the BRAR trial.
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
    tuning = 1, clipping = 0, urn_alpha = 0, randmethod = c("coin", "block", "urn"),
    return = c("allocations", "probabilities")
) {
  direction = match.arg(direction)
  postprobmethod = match.arg(postprobmethod)
  multiarm_method = match.arg(multiarm_method)
  randmethod = match.arg(randmethod)
  return = match.arg(return)

  N = nrow(trial_data)
  arms = length(unique(trial_data$Arm))

  allocations = numeric(N)
  allocation_probs_matrix = matrix(NA, nrow = N, ncol = arms)

  # Initialize alpha/beta with priors + observed outcomes
  current_alpha = priors[1, ]
  current_beta = priors[2, ]
  for (k in 1:arms) {
    y_k = sum(trial_data$Outcome[trial_data$Arm == k])
    n_k = sum(trial_data$Arm == k)
    current_alpha[k] = current_alpha[k] + y_k
    current_beta[k] = current_beta[k] + n_k - y_k
  }

  batches = unique(trial_data$Batch)

  for (b in batches) {
    block_idx = which(trial_data$Batch == b)
    block_size = length(block_idx)

    # Posterior probabilities
    if (postprobmethod == "simulation") {
      alloc_probs = posterior_bin_sim(current_alpha, current_beta, direction)
    } else {
      alloc_probs = posterior_bin_exact(current_alpha, current_beta, direction)
    }

    # Multi-arm adjustments
    if (arms > 2) {
      if (multiarm_method == "fixed") {
        if (alloc_probs[1] < 1/arms) {
          remaining_mass = 1 - 1/arms
          other_probs = alloc_probs[-1]
          other_probs = other_probs / sum(other_probs) * remaining_mass
          alloc_probs = c(1/arms, other_probs)
        }
      } else if (multiarm_method == "top2") {
        top2beta = 0.5
        alloc_probs = alloc_probs_T2TS(alloc_probs, top2beta)
      }
    }

    # Tuning
    if (tuning == 0) {
      alloc_probs_tuned = rep(1 / arms, arms)
    } else {
      alloc_probs_tuned = alloc_probs^tuning
      alloc_probs_tuned = alloc_probs_tuned / sum(alloc_probs_tuned)
    }

    # Clipping
    if (is.character(clipping) && clipping == "adaptive") {
      clipping_val = min(1/arms, (1/arms) * b^(-0.7))
    } else {
      clipping_val = max(0, clipping)
    }

    if (clipping_val > 0) {
      lower = clipping_val
      upper = 1 - (arms - 1) * clipping_val
      alloc_probs_tuned = pmax(alloc_probs_tuned, lower)
      alloc_probs_tuned = pmin(alloc_probs_tuned, upper)
      alloc_probs_tuned = alloc_probs_tuned / sum(alloc_probs_tuned)
    }

    # Randomization
    if (randmethod == "coin") {
      allocations[block_idx] = sample(1:arms, block_size, prob = alloc_probs_tuned, replace = TRUE)
    } else if (randmethod == "block") {
      counts = floor(alloc_probs_tuned * block_size)
      remainder = block_size - sum(counts)
      if (remainder > 0) {
        frac_parts = alloc_probs_tuned - counts / block_size
        add_idx = sample(1:arms, remainder, prob = frac_parts, replace = FALSE)
        counts = counts + tabulate(add_idx, nbins = arms)
      }
      alloc_vector = rep(1:arms, counts)
      allocations[block_idx] = sample(alloc_vector, block_size)
    } else if (randmethod == "urn") {
      count_in_block = rep(0, arms)
      for (i in 1:block_size) {
        weights = pmax(urn_alpha * alloc_probs_tuned + (i-1) * alloc_probs_tuned - count_in_block, 1)
        prob_now = weights / sum(weights)
        draw = sample(1:arms, 1, prob = prob_now)
        allocations[block_idx[i]] = draw
        count_in_block[draw] = count_in_block[draw] + 1
      }
    }

    allocation_probs_matrix[block_idx, ] = matrix(rep(alloc_probs_tuned, each = block_size), ncol = arms, byrow = FALSE)
  }

  if (return == "allocations") {
    return(allocations)
  } else {
    return(allocation_probs_matrix)
  }
}
