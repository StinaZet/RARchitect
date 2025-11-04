#' @title Simulate a Fixed Randomization (FR) Trial with Multiple Randomisation Methods
#'
#' @description
#' Simulate a multi-arm clinical trial using Fixed Randomisation (FR). Participants
#' are allocated to arms according to fixed allocation probabilities, using one of:
#' \code{"coin"} (independent draws), \code{"block"} (permuted block randomisation),
#' or \code{"urn"} (simple reinforcement urn / Pólya-like scheme).
#' Supports binary (Bernoulli), normal, and exponential outcomes and simulates trial
#' duration via a Poisson recruitment process (using your
#' \code{simulate_trial_duration_poisson_recruitment} function).
#'
#' @param outcome_type Character. `"binary"` or `"cont"`.
#' @param distribution Character. `"bernoulli"`, `"normal"`, or `"exponential"`.
#' @param arms Integer. Number of arms (>= 2).
#' @param N Integer. Total sample size.
#' @param direction Character. `"lower"` or `"higher"`.
#' @param known_var Logical. If `TRUE` treat normal variance as known (only applies to `distribution="normal"`).
#' @param modelpar Numeric vector or matrix. True data-generating parameters, where the first entry is for the control arm (see \code{simulate_brar_trial} documentation).
#' @param allocation_probs Numeric vector of length `arms` or character string. Fixed allocation probabilities (will be normalized).
#'   Defaults to equal allocation. If set to the string \code{"dunnett"} (and \code{arms >= 2}), the probabilities are
#'   set proportional to sqrt{k}:1:...:1, where k = arms-1 is the number of active arms (Arm 1 is control).
#' @param randmethod Character. Randomisation method: `"coin"` (default), `"block"`, or `"urn"`.
#' @param blocksize Integer. Block size for `"block"` randomisation. Ignored for `"coin"` and `"urn"`.
#' @param urn_alpha Numeric. Initial urn mass parameter for `"urn"` (controls initial ball counts). Default 3.
#' @param recruitment_rate Numeric. Poisson rate lambda for recruitment (per time unit). Default 100000.
#' @param observation_delay Numeric. Fixed delay from recruitment to outcome observation. Default 0.
#'
#' @return A data.frame with N rows and columns:
#' \itemize{
#'  \item `Participant`: 1..N
#'  \item `Recruitment time`
#'  \item `Outcome time`
#'  \item `Arm`: assigned arm (1..arms)
#'  \item `Outcome`: observed outcome (0/1 for binary, continuous for cont)
#'  \item `AP arm X`: fixed allocation probability for arm X, repeated for each participant
#' }
#' @export
#'
#' @examples
#' set.seed(301)
#' simulate_fr_trial(
#'   outcome_type = "binary", distribution = "bernoulli",
#'   arms = 2, N = 100, modelpar = c(0.6, 0.4),
#'   allocation_probs = c(0.5, 0.5),
#'   randmethod = "coin", recruitment_rate = 5, observation_delay = 10)
#'
#' set.seed(302)
#' simulate_fr_trial(
#'   outcome_type = "cont", distribution = "normal",
#'   arms = 3, N = 150,
#'   modelpar = matrix(c(10, 8, 6, 2, 2, 2), nrow = 2, byrow = TRUE),
#'   allocation_probs = c(0.5, 0.3, 0.2),
#'   randmethod = "block", blocksize = 10,
#'   recruitment_rate = 10, observation_delay = 5)
#'
#' set.seed(303)
#' simulate_fr_trial(
#'  outcome_type = "binary", distribution = "bernoulli",
#'  arms = 3, N = 120, modelpar = c(0.3, 0.5, 0.7),
#'  allocation_probs = "dunnett", # <-- NEW EXAMPLE
#'  randmethod = "urn", urn_alpha = 5,
#'  recruitment_rate = 8, observation_delay = 2)
#'
simulate_fr_trial <- function(outcome_type = c("binary", "cont"),
                              distribution = c("bernoulli", "normal", "exponential"),
                              arms, N, direction = c("lower", "higher"),
                              known_var = FALSE,
                              modelpar,
                              allocation_probs = NULL,
                              randmethod = c("coin", "block", "urn"),
                              blocksize,
                              urn_alpha = 3,
                              recruitment_rate = 100000,
                              observation_delay = 0) {

  # Argument matching and basic checks
  outcome_type = match.arg(outcome_type)
  distribution = match.arg(distribution)
  direction = match.arg(direction)
  randmethod = match.arg(randmethod)

  if (!is.numeric(arms) || arms %% 1 != 0 || arms <= 1) {
    stop("'arms' must be an integer greater than 1.")
  }
  if (!is.numeric(N) || N %% 1 != 0 || N <= 0) {
    stop("'N' must be a positive integer.")
  }
  if (recruitment_rate <= 0 || !is.numeric(recruitment_rate)) {
    stop("'recruitment_rate' must be a positive number.")
  }
  if (observation_delay < 0 || !is.numeric(observation_delay)) {
    stop("'observation_delay' must be a non-negative number.")
  }
  if (!is.numeric(urn_alpha) || urn_alpha <= 0) {
    stop("'urn_alpha' must be positive.")
  }

  # Validate distribution vs outcome_type
  if (outcome_type == "binary" && distribution != "bernoulli") {
    stop("For outcome_type = 'binary', distribution must be 'bernoulli'.")
  }
  if (outcome_type == "cont" && !(distribution %in% c("normal", "exponential"))) {
    stop("For outcome_type = 'cont', distribution must be 'normal' or 'exponential'.")
  }

  # Normalize or set allocation probabilities
  if (is.null(allocation_probs)) {
    # Default to equal allocation if NULL
    allocation_probs = rep(1 / arms, arms)
  } else if (is.character(allocation_probs) && allocation_probs == "dunnett") {
    # Implement Dunnett's square root rule
    if (arms < 2) stop("Dunnett allocation requires at least two arms (control + 1 active).")
    # Number of active arms
    n_active = arms - 1
    # Allocation for control arm 1 is sqrt(n_active), active arms are 1.
    raw_probs = c(sqrt(n_active), rep(1, n_active))
    # Normalize
    allocation_probs = raw_probs / sum(raw_probs)

  } else if (is.numeric(allocation_probs)) {
    # Handle custom numeric vector input.
    if (length(allocation_probs) != arms) {
      stop("'allocation_probs' must be a numeric vector with length equal to 'arms'.")
    }
    if (any(allocation_probs < 0)) stop("'allocation_probs' must be non-negative.")
    if (sum(allocation_probs) == 0) stop("'allocation_probs' must sum to a positive number.")
    allocation_probs = allocation_probs / sum(allocation_probs)
  } else {
    stop("'allocation_probs' must be NULL, the string \"dunnett\", or a numeric vector of length 'arms'.")
  }

  # --- Handle blocksize requirements ---
  if (randmethod == "block") {
    # Must be provided for block randomisation
    if (missing(blocksize)) stop("'blocksize' must be specified for block randomisation")
    if (!is.numeric(blocksize) || blocksize %% 1 != 0 || blocksize <= 0) {
      stop("'blocksize' must be a positive integer")
    }
    # If blocksize > N, reduce to single block
    if (blocksize > N) blocksize = N
  } else {
    # For coin or urn, ignore blocksize
    blocksize = NULL
  }

  # --- Generate Arm assignments according to randmethod ---
  Arm = integer(N)

  if (randmethod == "coin") {
    # Independent draws with fixed probabilities (allocation_probs)
    Arm = sample(1:arms, size = N, replace = TRUE, prob = allocation_probs)

  } else if (randmethod == "block") {
    # --- Permuted-block randomisation ---
    # Goal: Assign participants in "blocks" to approximately match allocation_probs
    # within each block, then shuffle to prevent predictability.
    n_blocks = ceiling(N / blocksize)
    idx = 1L
    for (b in seq_len(n_blocks)) {
      # Determine the size of this block (last block may be smaller)
      this_block_size = if (b < n_blocks) blocksize else (N - (b - 1) * blocksize)

      # Desired counts per arm based on allocation_probs
      raw_counts = allocation_probs * this_block_size
      base_counts = floor(raw_counts)              # integer part
      remainder = this_block_size - sum(base_counts)  # leftover slots

      # Distribute remainder probabilistically according to fractional parts
      if (remainder > 0) {
        fractional = raw_counts - base_counts
        # Safety: if all fractional parts are 0, distribute randomly
        if (all(fractional == 0)) fractional = rep(1/arms, arms)
        add_indices = sample(1:arms, size = remainder, prob = fractional, replace = TRUE)
        add_tab = tabulate(add_indices, nbins = arms)
        counts = base_counts + add_tab
      } else {
        counts = base_counts
      }

      # Construct block vector: repeat arm numbers according to counts
      block_vec = rep(seq_len(arms), times = counts)
      # Shuffle assignments within block to prevent predictability
      block_vec = sample(block_vec, size = length(block_vec), replace = FALSE)
      # Assign to main Arm vector
      Arm[idx:(idx + this_block_size - 1)] = block_vec
      idx = idx + this_block_size
    }

  } else if (randmethod == "urn") {
    # --- Simple Pólya-like urn randomisation ---
    # Goal: Use reinforced randomisation so frequently drawn arms get slightly higher chance
    # Start with initial integer "ball counts" proportional to allocation_probs * urn_alpha
    initial_raw = allocation_probs * urn_alpha
    balls = pmax(1, round(initial_raw))  # ensure at least 1 ball per arm

    # For each participant
    for (i in seq_len(N)) {
      # Current allocation probabilities are proportional to ball counts
      probs_now = balls / sum(balls)
      # Draw an arm based on current probabilities
      draw = sample(1:arms, size = 1, prob = probs_now)
      Arm[i] = draw
      # Reinforce: add 1 ball to the drawn arm
      balls[draw] = balls[draw] + 1L
    }

  } else {
    stop("Unsupported 'randmethod'.")
  }


  # --- Generate outcomes according to distribution & modelpar ---
  if (distribution == "bernoulli") {
    if (!is.numeric(modelpar) || length(modelpar) != arms) {
      stop("For 'bernoulli' distribution, 'modelpar' must be a numeric vector of length 'arms' containing success probabilities.")
    }
    Outcome = stats::rbinom(N, size = 1, prob = modelpar[Arm])

  } else if (distribution == "normal") {
    if (!is.matrix(modelpar) || nrow(modelpar) != 2 || ncol(modelpar) != arms) {
      stop("For 'normal' distribution, 'modelpar' must be a 2-row matrix: row1 = means, row2 = sds, with ncol == arms.")
    }
    means = modelpar[1, ]
    sds = modelpar[2, ]
    Outcome = stats::rnorm(N, mean = means[Arm], sd = sds[Arm])

  } else if (distribution == "exponential") {
    if (!is.numeric(modelpar) || length(modelpar) != arms) {
      stop("For 'exponential' distribution, 'modelpar' must be a numeric vector of length 'arms' containing rates.")
    }
    Outcome = stats::rexp(N, rate = modelpar[Arm])

  } else {
    stop("Unsupported 'distribution'.")
  }

  # --- Trial duration simulation using your existing helper ---
  # We'll treat each participant as its own block for the duration sim (this mirrors previous FR implementation).
  block_sizes = rep(1L, N)
  duration_results = simulate_trial_duration_poisson_recruitment(
    N = N,
    block_sizes = block_sizes,
    poisson_recruitment_rate_per_unit_time = recruitment_rate,
    outcome_observation_delay_per_patient = observation_delay,
    num_simulations = 1
  )

  # --- Build output dataframe ---
  df = data.frame(
    Participant = seq_len(N),
    `Recruitment time` = duration_results[, 2],
    `Outcome time` = duration_results[, 4],
    Arm = Arm,
    Outcome = Outcome,
    stringsAsFactors = FALSE
  )

  # Add fixed allocation probabilities as columns (AP arm X)
  ap_mat = matrix(rep(allocation_probs, each = N), nrow = N, byrow = FALSE)
  colnames(ap_mat) = paste0("AP arm ", seq_len(arms))
  df = cbind(df, as.data.frame(ap_mat, stringsAsFactors = FALSE))

  return(df)
}
