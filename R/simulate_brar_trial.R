#' @title Simulate a Response Adaptive Randomization (RAR) Trial for Binary or Normal Outcomes
#'
#' @description
#' This function simulates a multi-arm clinical trial using Bayesian Response-Adaptive
#' Randomization (BRAR)/Thompson Sampling with either a Beta-Bernoulli model for
#' binary outcomes or a Normal-Normal conjugate model for normal outcomes
#' (assuming known population variance). It also simulates trial
#' duration based on a Poisson recruitment process.
#'
#' @param outcome_type Character. Specifies the type of outcome to simulate.
#' Must be either `"binary"`, `"normal"`, or `"exponential"`.
#' @param arms Numeric. Number of arms in the trial.
#' @param N Numeric. Total sample size for the trial.
#' @param blocksize Numeric. (Fixed) size of each block of participants for the adaptive randomization.
#' @param known_var Logical. Indicates if the variance for the `"normal"` outcome type
#' is known (`TRUE`) or unknown (`FALSE`). Defaults to `FALSE`.
#' @param priors Matrix. Defines the prior parameters for the Bayesian model of each arm.
#' The structure of this matrix depends on `outcome_type` and `known_var`:
#' \itemize{
#' \item If `outcome_type = "binary"`: A 2-row matrix where the first row contains the alpha
#' parameters and the second row contains the beta parameters for the Beta
#' distributions of each arm. The number of columns must match `arms`.
#' e.g., `matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)` for two arms with uniform Beta(1,1) priors.
#' \item If `outcome_type = "normal"` and `known_var = TRUE`: A 2-row matrix
#' containing the Normal prior parameters for the mean of each arm . Each column
#' corresponds to an arm. Rows should be:
#' 1. Prior mean for the arm's mean.
#' 2. Prior standard deviation of the arm's mean.
#' Example for two arms with N(0,10) priors for means:
#' `matrix(c(0, 0, 10, 10), nrow = 2, byrow = TRUE)`
#' \item If `outcome_type = "normal"` and `known_var = FALSE`:  A 4-row matrix
#' containing the prior parameters for the mean of each arm . Each column
#' corresponds to an arm. Rows should be:
#' 1. mu0: prior mean for the arm's mean.
#' 2. kappa: "strength" of the prior mean (like a prior sample size).
#' 3. alpha: prior shape for the variance.
#' 4. beta: prior scale for the variance.
#' Example for two arms with mu0=0, kappa=0.1, alpha=1, beta=1:
#' `matrix(c(0, 0, 0.1, 0.1, 1, 1, 1, 1), nrow = 4, byrow = TRUE)`
#' \item If `outcome_type = "exponential"`: A 2-row matrix where the first row contains the shape
#' parameters and the second row contains the rate parameters for the Gamma
#' distributions of each arm. The number of columns must match `arms`.
#' }
#' @param modelpar Matrix or Numeric vector. Specifies the true parameters for each arm's
#' data generating process. The structure depends on `outcome_type`:
#' \itemize{
#' \item If `outcome_type = "binary"`: A numeric vector with the true success
#' probabilities for each Bernoulli arm. e.g., `c(0.6, 0.4)`.
#' \item If `outcome_type = "normal"`: A 2-row matrix specifying the true parameters
#' for each normal arm. Each column corresponds to an arm. Rows should be:
#' 1. True mean for the arm.
#' 2. True standard deviation for the arm. It is assumed to be
#' the known population standard deviation for the purpose of posterior updates.
#' Example for two arms: `matrix(c(10, 8, 2, 2), nrow = 2, byrow = TRUE)`
#' \item If `outcome_type = "exponential"`: A numeric vector with the rates for
#' each exponential arm. e.g., `c(0.2, 0.1)`.
#' }
#' @param tuning Numeric. A parameter (referred to as 'c' or 'gamma' in the
#' Wathen and Thall, 2017 paper) that shrinks allocation probabilities
#' towards equal randomization.
#' \itemize{
#' \item If `tuning = 1`, no shrinking occurs; original adaptive randomization probabilities are used.
#' \item As `tuning` approaches 0, probabilities shrink towards `1/arms` (equal randomization for `arms` arms).
#' \item A common value is 0.5. Defaults to 1.
#' }
#' @param clipping Character or Numeric. Forces a minimum and maximum probability
#' of selecting each arm. Defaults to 0.
#' \itemize{
#' \item If `clipping = 0` (Numeric), no clipping is applied.
#' \item If `clipping > 0` (Numeric), each arm's allocation probability is forced
#' to be at least `clipping` and at most `1 - (arms - 1) * clipping`.
#' Probabilities are then re-normalized to sum to 1.
#' It is **required** that `clipping * arms <= 1` for consistent bounds.
#' \item If `clipping = "adaptive"` (Character), adaptive clipping is applied as
#' `rho_min,t = 1/arms * (batch_number)^(-0.7)` (as in Hadad et al., 2021).
#' }
#' @param burnin Numeric. The number of initial participants allocated
#' before standard block processing begins. These participants are treated
#' as a single initial block, allocated using equal randomization. Defaults to 0.
#' @param ensure_all_arms_sampled Logical. If TRUE, ensures that at least
#' one participant is allocated to each arm within each block (if `blocksize >= arms`).
#' This is useful for ensuring initial exploration. Defaults to FALSE.
#' @param postprobmethod Character. (Applicable only when `outcome_type = "binary"`).
#' Method for calculating posterior probabilities. Can be `"simulation"` or `"exact"`. Defaults to `"simulation"`.
#' @param recruitment_rate Numeric. The rate parameter (lambda)
#' for the Poisson process that governs patient recruitment for trial duration simulation.
#' This is the average number of patients arriving per unit of time (e.g., patients per month).
#' Must be a positive value. Defaults to 100000.
#' @param observation_delay Numeric. The fixed time from a participant's
#' enrollment until their outcome is observed, used for trial duration simulation. Defaults to 0.
#'
#' @return A data frame with `N` rows and the following columns:
#' \itemize{
#' \item `Batch`: The block number for each participant.
#' \item `Recruitment time`: The time point at which each participant is recruited.
#' \item `Outcome time`: The time point at which each participant's outcome is observed.
#' \item `Arm`: The arm selected for the participant (1 to `arms`).
#' \item `Outcome`: The outcome associated with the selected arm (binary 0/1 for "binary",
#' continuous for "normal").
#' \item `AP arm X`: For each arm X, the allocation probability of Arm X for that block.
#' }
#' @export
#'
#' @examples
#' # Example 1: Simulate a Binary Outcome Trial
#' set.seed(101)
#' results_binary <- simulate_brar_trial(
#' outcome_type = "binary",
#' arms = 2, N = 100, blocksize = 10,
#' modelpar = c(0.6, 0.4),
#' priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
#' tuning = 1, clipping = 0, burnin = 0, ensure_all_arms_sampled = FALSE,
#' postprobmethod = "exact", # Or "simulation"
#' recruitment_rate = 5,
#' observation_delay = 30)
#' head(results_binary)
#'
#' # Example 2: Simulate a Normal Outcome Trial
#' set.seed(102)
#' model_params_norm <- matrix(c(10, 8, 2, 2), nrow = 2, byrow = TRUE)
#' prior_params_norm <- matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE)
#'
#' results_normal <- simulate_brar_trial(
#' outcome_type = "normal",
#' arms = 2, N = 200, blocksize = 10, known_var = TRUE,
#' priors = prior_params_norm,
#' modelpar = model_params_norm,
#' tuning = 1, clipping = 0, burnin = 10, ensure_all_arms_sampled = FALSE,
#' postprobmethod = "exact", # Or "simulation"
#' recruitment_rate = 10,
#' observation_delay = 15)
#' head(results_normal)
#'
#' # Example 3: Binary trial with tuning and clipping and uneven sample size
#' set.seed(103)
#' results_binary_tuned <- simulate_brar_trial(
#' outcome_type = "binary",
#' arms = 2, N = 149, burnin =15, blocksize = 15,
#' priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
#' modelpar = c(0.7, 0.5),
#' tuning = 0.5,
#' clipping = 0.05,
#' ensure_all_arms_sampled = TRUE,
#' postprobmethod = "exact",
#' recruitment_rate = 8,
#' observation_delay = 20)
#' tail(results_binary_tuned)
#'
#' # Example 4: Simulate a Normal Outcome Trial (unknown variance)
#' set.seed(104)
#' model_unvar <- matrix(c(10, 8, 2, 2), nrow = 2, byrow = TRUE)
#' # priors: mu0=0, kappa=0.1, alpha=1, beta=1 for each arm
#' prior_unvar <- matrix(c(0, 0,
#'                         0.1, 0.1,
#'                         1, 1,
#'                         1, 1),
#'                       nrow = 4, byrow = TRUE)
#'
#' results_normal_unvar <- simulate_brar_trial(
#'   outcome_type = "normal",
#'   arms = 2, N = 200, blocksize = 20, known_var = FALSE,
#'   priors = prior_unvar,
#'   modelpar = model_unvar,
#'   tuning = 1, clipping = 0, burnin = 10,
#'   ensure_all_arms_sampled = FALSE,
#'   postprobmethod = "simulation",
#'   recruitment_rate = 12,
#'   observation_delay = 10)
#' head(results_normal_unvar)
#'
#' #' # Example 5: Simulate an Exponential Outcome Trial
#' # True rates for arms: lambda1 = 0.2 (mean = 5), lambda2 = 0.1 (mean = 10)
#' set.seed(105)
#' prior_params_exp <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE) # Gamma(1,1) priors
#' true_rates_exp <- c(0.2, 0.1)
#'
#' results_exponential <- simulate_brar_trial(
#'   outcome_type = "exponential",
#'   arms = 2, N = 120, blocksize = 12,
#'   priors = prior_params_exp,
#'   modelpar = true_rates_exp,
#'   tuning = 1, clipping = 0, burnin = 0,
#'   ensure_all_arms_sampled = FALSE,
#'   postprobmethod = "simulation",
#'   recruitment_rate = 6,
#'   observation_delay = 5)
#' head(results_exponential)
#'
simulate_brar_trial <- function(outcome_type = c("binary", "normal"),
                                arms = 2, N, blocksize, known_var = FALSE,
                                priors, modelpar, tuning = 1, clipping = 0, burnin = 0,
                                ensure_all_arms_sampled = FALSE,
                                postprobmethod = "simulation",
                                recruitment_rate = 100000,
                                observation_delay = 0) {

  # Input validation for outcome_type
  outcome_type <- match.arg(outcome_type)

  # Validate common parameters
  if (N <= 0 || !is.numeric(N) || N %% 1 != 0) {
    stop("Parameter 'N' must be a positive integer.")
  }
  if (blocksize <= 0 || !is.numeric(blocksize) || blocksize %% 1 != 0) {
    stop("Parameter 'blocksize' must be a positive integer.")
  }
  if (arms <= 1 || !is.numeric(arms) || arms %% 1 != 0) {
    stop("Parameter 'arms' must be an integer greater than 1.")
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
  if (burnin < 0 || !is.numeric(burnin) || burnin %% 1 != 0) {
    stop("The 'burnin' parameter must be a non-negative integer.")
  }
  if (ensure_all_arms_sampled && blocksize < arms && blocksize > 0) {
    stop("If 'ensure_all_arms_sampled' is TRUE, 'blocksize' must be at least 'arms' (or 0 if burnin covers initial).")
  }
  if (ensure_all_arms_sampled && blocksize == 1 && arms > 1) {
    stop("You have a block size of 1 and more than one arm, but you asked to ensure each arm is sampled at least once in each block! This is impossible. Please change the block size or set 'ensure_all_arms_sampled = FALSE'.")
  }
  # Validation for renamed parameters
  if (recruitment_rate <= 0 || !is.numeric(recruitment_rate)) {
    stop("Parameter 'recruitment_rate' must be a positive number for the trial duration simulation (as it represents a Poisson rate).")
  }
  if (observation_delay < 0 || !is.numeric(observation_delay)) {
    stop("Parameter 'observation_delay' must be a positive number for the trial duration simulation.")
  }

  if (is.null(postprobmethod)) {
    stop("'postprobmethod' must be specified ('simulation' or 'exact').")
  }
  if (!(postprobmethod %in% c("simulation", "exact"))) {
    stop("Invalid 'postprobmethod'. Must be 'simulation' or 'exact'.")
  }

  # Delegate to specific simulation functions based on outcome_type
  if (outcome_type == "binary") {
    if (known_var == TRUE) {
      warning("You have chosen the variance to be known, this option is only available for normal outcome.")
    }
    #if (postprobmethod == "simulation") {
    #  warning("You have chosen the method 'simulation' for calculating posterior probabilities. It is possible to calculate the posterior probabilities exactly for this type of outcome variable.")
    #}

    results <- .simulate_brar_trial_binary(
      arms = arms, N = N, blocksize = blocksize,
      priors = priors, modelpar = modelpar, tuning = tuning,
      clipping = clipping, burnin = burnin,
      ensure_all_arms_sampled = ensure_all_arms_sampled,
      postprobmethod = postprobmethod
    )
  } else if (outcome_type == "normal") {

    #if (postprobmethod == "simulation") {
    #  warning("You have chosen the method 'simulation' for calculating posterior probabilities. It is possible to calculate the posterior probabilities exactly for this type of outcome variable.")
    #}

    results <- .simulate_brar_trial_normal(
      arms = arms, N = N, blocksize = blocksize, known_var = known_var,
      priors = priors, modelpar = modelpar, tuning = tuning,
      clipping = clipping, burnin = burnin,
      ensure_all_arms_sampled = ensure_all_arms_sampled,
      postprobmethod = postprobmethod
    )
  } else if (outcome_type == "exponential") {

    #if (postprobmethod == "simulation") {
    #  warning("You have chosen the method 'simulation' for calculating posterior probabilities. It is possible to calculate the posterior probabilities exactly for this type of outcome variable.")
    #}

    results <- .simulate_brar_trial_exp(
      arms = arms, N = N, blocksize = blocksize,
      priors = priors, modelpar = modelpar, tuning = tuning,
      clipping = clipping, burnin = burnin,
      ensure_all_arms_sampled = ensure_all_arms_sampled,
      postprobmethod = postprobmethod
    )
  } else {
    # This block should ideally not be reached due to match.arg, but as a safeguard.
    stop("Invalid 'outcome_type'. Must be 'binary', 'normal' or 'exponential'.")
  }

  # --- Call simulate_trial_duration_poisson_recruitment ---
  # Construct block_sizes for duration simulation based on burnin and blocksize
  remaining_N <- N - burnin

  if (remaining_N < 0) {
    stop("Burn-in period ('burnin') cannot be greater than total sample size ('N').")
  }

  #if (remaining_N > 0 && remaining_N %% blocksize != 0) {
  #  stop("The number of participants after the burn-in period (N - burnin) must be a multiple of 'blocksize' for consistent trial duration simulation. Please adjust N, burnin, or blocksize.")
  #}

  if (N > 0) { # Only run duration simulation if N is positive
    num_main_blocks <- remaining_N / blocksize

    if (burnin > 0) {
      # The first 'burnin' participants are treated as one initial block
      # The remaining participants are grouped into blocks of 'blocksize'
      block_sizes_for_duration_sim <- c(burnin, rep(blocksize, num_main_blocks))
    } else {
      # If no burnin, all participants are in blocks of 'blocksize'
      block_sizes_for_duration_sim <- rep(blocksize, N / blocksize)
    }

    # Call the duration simulation function with new parameter names
    duration_results <- simulate_trial_duration_poisson_recruitment(
      N = N,
      block_sizes = block_sizes_for_duration_sim,
      poisson_recruitment_rate_per_unit_time = recruitment_rate,
      outcome_observation_delay_per_patient = observation_delay,
      num_simulations = 1
    )

    # Merge timing results with the main simulation results
    results = cbind(results[,1], duration_results[,2], duration_results[,4], results[,-1])
    colnames(results) = c("Batch", "Recruitment time", "Outcome time", "Arm", "Outcome", paste0("AP arm ", 1:arms))

  } else { # N is 0, add empty columns for consistency if 'results' is not already empty with all expected columns
    # Assuming .simulate_brar_trial_binary and _normal handle N=0 correctly returning
    # an empty data frame with appropriate columns, we ensure the duration columns are added.
    # If results is already a data.frame with the expected RAR trial columns but 0 rows:
    if (is.data.frame(results) && nrow(results) == 0) {
      results$EstimatedRecruitmentTime <- numeric(0)
      results$StdDevRecruitmentTime <- numeric(0)
      results$EstimatedObservationTime <- numeric(0)
      results$StdDevObservationTime <- numeric(0)
    } else {
      # Fallback for unexpected empty 'results' structure
      results <- data.frame(
        Batch = numeric(0),
        Arm = numeric(0),
        Outcome = numeric(0),
        # Assuming AlloProb_ArmX will be handled by the specific outcome functions for N=0
        EstimatedRecruitmentTime = numeric(0),
        StdDevRecruitmentTime = numeric(0),
        EstimatedObservationTime = numeric(0),
        StdDevObservationTime = numeric(0)
      )
    }
  }
  # --- End Call simulate_trial_duration_poisson_recruitment ---

  return(results)
}
