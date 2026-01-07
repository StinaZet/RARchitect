#' @title Analyze a Single BRAR Trial
#'
#' @description
#' Wrapper function to analyze a single BRAR or Fixed Randomization (FR) trial replicate
#' with binary, normal, or exponential outcomes.
#'
#' Returns estimated treatment effects, p-values, and confidence intervals.
#' Also reports patient benefit (number of patients on the estimated best arm).
#'
#' @param outcome_type Character. Specifies the type of outcome to simulate.
#' Must be either `"binary"` or `"cont"`.
#' @param distribution Character. Specifies the distribution of the outcome.
#' Must be either `"bernoulli"`, `"normal"`, or `"exponential"`.
#' @param trial_data Data frame containing trial results.
#'   Must include columns 'Arm' and 'Outcome'.
#' @param N Total sample size.
#' @param arms Number of treatment arms.
#' @param direction `"lower"` or `"higher"` (indicates which outcome is favorable).
#' @param priors Matrix of prior parameters for Bayesian estimation (2 x K for binary).
#' @param known_var Logical, for normal outcomes with known variance.
#' @param estimation_method `"MLE"`, `"IPW"`, or `"post_mean"`.
#' @param test_method `"standard"`, `"exact"`, `"randomization"`, `"AP"`, or `"simulation"`.
#' @param CI_method `"standard"` or `"simulation"`.
#' @param effect_measure Currently only `"riskdifference"` supported.
#' @param multiple_tests Logical; if `TRUE`, test all experimental arms vs control with Bonferroni adjustment.
#' @param onesided Logical; if `TRUE`, perform one-sided test.
#' @param critval Numerical. if the critical value for the simulation-based test is known it can be specified here.
#' If it is unknown, leave it for `NULL` and the critical value is calculated in the function.
#' @param B Number of simulations/permutations for simulation/randomization-based tests.
#' @param alpha Significance level (used for CI and tests).
#' @param prob_threshold Posterior probability threshold (for Bayesian decision rules).
#' @param ... Additional arguments passed to test methods (e.g. `randmethod`, `blocksize`, `postprobmethod`, `multiarm_method`).
#'
#' @return A list with:
#' \describe{
#'   \item{test_results_df}{Data frame with estimated effects, confidence intervals, and p-values.}
#'   \item{patient_benefit}{Number of patients allocated to the best estimated arm.}
#'   \item{mean_outcome}{Overall mean outcome across patients.}
#' }
#'
#' @details
#' **Notes:**
#' - MLE estimation may be biased under adaptive randomization.
#' - Standard test and standard confidence intervals may not be valid for BRAR data.
#' - For adaptive designs, prefer "post_mean" or "IPW" estimation and "randomization" or "AP" testing.
#'
#' @examples
#' # --- Example setup: simulate a BRAR trial ---
#' set.seed(103)
#' results_binary_multiarm = simulate_brar_trial(
#'   outcome_type = "binary",
#'   distribution = "bernoulli",
#'   direction = "higher",
#'   arms = 3, N = 150, burnin = 15, blocksize = 1,
#'   priors = matrix(c(1, 1, 1, 1, 1, 1), nrow = 2, byrow = TRUE),
#'   modelpar = c(0.5, 0.7, 0.8),
#'   tuning = 1,
#'   clipping = 0,
#'   randmethod = "coin",
#'   postprobmethod = "simulation",
#'   recruitment_rate = 8,
#'   observation_delay = 20)
#'
#' # --- Example 1: Simulation-based test ---
#' analyze_brar_trial(outcome_type = "binary",
#'   trial_data = results_binary_multiarm,
#'   N = 150, arms = 3, direction = "higher",
#'   priors = matrix(c(1, 1, 1, 1, 1, 1), nrow = 2, byrow = TRUE),
#'   distribution = "bernoulli", burnin = 15,
#'   postprobmethod = "simulation",
#'   blocksize = 1, estimation_method = "IPW",
#'   test_method = "simulation", CI_method = "simulation",
#'   B = 10)
#'
#' # --- Example 2: standard test (not recommended for BRAR) ---
#' analyze_brar_trial(outcome_type = "binary",
#'   trial_data = results_binary_multiarm,
#'   N = 150, arms = 3, direction = "higher",
#'   priors = matrix(c(1, 1, 1, 1, 1, 1), nrow = 2, byrow = TRUE),
#'   distribution = "bernoulli", postprobmethod = "simulation",
#'   estimation_method = "MLE",
#'   test_method = "standard", CI_method = "standard")
#'
#' # --- Example 3: AP test (Allocation-Probability test) ---
#' analyze_brar_trial(outcome_type = "binary",
#'   trial_data = results_binary_multiarm,
#'   N = 150, arms = 3, direction = "higher",
#'   priors = matrix(c(1, 1, 1, 1, 1, 1), nrow = 2, byrow = TRUE),
#'   distribution = "bernoulli", postprobmethod = "simulation",
#'   estimation_method = "IPW",
#'   test_method = "AP", CI_method = "simulation",
#'   randmethod = "coin", blocksize = 1, B = 10,
#'   postprobmethod = "simulation")
#'
#' @export
analyze_brar_trial <- function(outcome_type = c("binary", "cont"),
                               distribution = c("bernoulli", "normal", "exponential"),
                               trial_data, N, arms, direction = c("lower", "higher"),
                               priors, known_var = NULL,
                               estimation_method = c("MLE", "IPW", "post_mean"),
                               test_method = c("standard", "exact", "randomization", "AP", "simulation"),
                               CI_method = c("standard", "simulation"),
                               effect_measure = "riskdifference",
                               multiple_tests = FALSE, onesided = TRUE,
                               critval = NULL, B = 10000,
                               alpha = 0.05, prob_threshold = NULL, ...) {

  # --- Validate arguments ---
  outcome_type = match.arg(outcome_type)
  distribution = match.arg(distribution)
  direction = match.arg(direction)
  estimation_method = match.arg(estimation_method)
  test_method = match.arg(test_method)
  CI_method = match.arg(CI_method)
  effect_measure = match.arg(effect_measure)

  # Error messages for the column names in the trial_data dataset.
  required_cols = c("Outcome", "Arm")
  missing_cols = setdiff(required_cols, names(trial_data))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The input dataset 'trial_data' is missing required column(s): ",
        paste(missing_cols, collapse = ", "),
        ".\nPlease ensure the dataset contains at least 'Outcome' and 'Arm' columns."
      )
    )
  }

  # AP test only works for one-sided alternative hypotheses.
  if (test_method == "AP" && onesided == FALSE) {
    stop("The AP test can only be used to test one-sided hypothesis. Please choose another hypothesis test or another alternative hypothesis.")
  }


  # Check that the columns with allocation probabilities are named correctly.
  if (estimation_method == "IPW" || test_method == "AP" || CI_method == "IPW")
  {
    # Check for correct 'AP arm X' column naming
    ap_cols = grep("^AP arm [0-9]+$", names(trial_data), value = TRUE)

    # Detect arms based on pattern
    if (length(ap_cols) == 0) {
      stop(
        "No columns matching the pattern 'AP arm X' were found in 'trial_data'.\n",
        "These columns are required for IPW estimation (allocation probabilities)."
      )
    }

    # Check if we have an AP arm column for *each* arm
    expected_ap_names = paste0("AP arm ", seq_len(arms))
    missing_ap = setdiff(expected_ap_names, ap_cols)
    if (length(missing_ap) > 0) {
      stop(
        paste0(
          "Missing allocation probability columns: ",
          paste(missing_ap, collapse = ", "),
          ".\nExpected columns: ",
          paste(expected_ap_names, collapse = ", "),
          ".\nPlease ensure the dataset contains one 'AP arm X' column for each arm."
        )
      )
    }
  }



  # --- Method warnings for BRAR data ---
  if (estimation_method == "MLE")
    warning("MLE estimation may be biased under BRAR. Consider 'IPW' or 'post_mean'.")
  if (test_method == "standard")
    warning("Standard test may not be appropriate for BRAR data. Consider 'randomization', 'simulation', or 'AP'.")
  if (CI_method == "standard")
    warning("Standard confidence intervals may not be appropriate for BRAR data. Consider 'simulation' CIs.")



  # --- Dispatch to correct helper based on distribution ---
  results = switch(
    distribution,
    "bernoulli" = .analyze_brar_trial_binary(
      trial_data = trial_data, priors = priors, N = N, arms = arms, direction = direction,
      estimation_method = estimation_method, test_method = test_method, CI_method = CI_method,
      effect_measure = effect_measure, alpha = alpha, multiple_tests = multiple_tests, B = B,
      onesided = onesided, ...
    ),
    "normal" = .analyze_brar_trial_normal(
      trial_data = trial_data, priors = priors, N = N, arms = arms, direction = direction,
      known_var = known_var,  B = B,
      estimation_method = estimation_method, test_method = test_method,
      alpha = alpha, multiple_tests = multiple_tests, onesided = onesided, ...
    ),
    "exponential" = .analyze_brar_trial_exp(
      trial_data = trial_data, priors = priors, N = N, arms = arms, direction = direction,
      estimation_method = estimation_method, test_method = test_method,  B = B,
      alpha = alpha, multiple_tests = multiple_tests, onesided = onesided, ...
    )
  )

  return(results)
}
