#' @title Analyze a Single BRAR Trial
#'
#' @description
#' Wrapper function to analyze a single BRAR or Fixed Randomization (FR) trial replicate
#' with binary, normal, or exponential outcomes.
#'
#' Returns estimated treatment effects, p-values, and confidence intervals.
#' Also reports patient benefit (number of patients on the estimated best arm).
#'
#' @param trial_data Data frame containing trial results.
#'   Must include columns 'Arm' and 'Outcome'.
#' @param N Total sample size.
#' @param arms Number of treatment arms.
#' @param direction "lower" or "higher" (indicates which outcome is favorable).
#' @param priors Matrix of prior parameters for Bayesian estimation (2 x K for binary).
#' @param distribution "bernoulli", "normal", or "exponential".
#' @param known_var Logical, for normal outcomes with known variance.
#' @param estimation_method "MLE", "IPW", or "post_mean" (for binary outcomes).
#' @param test_method "wald", "exact", "randomization", "AP", or "simulation".
#' @param CI_method "wald" or "simulation".
#' @param effect_measure Currently only "riskdifference" supported for binary outcomes.
#' @param multiple_tests Logical; if TRUE, test all experimental arms vs control with Bonferroni adjustment.
#' @param onesided Logical; if TRUE, perform one-sided test.
#' @param alpha Significance level (used for CI and tests).
#' @param prob_threshold Posterior probability threshold (for Bayesian decision rules).
#' @param ... Additional arguments passed to test methods (e.g. `B`, `policy_file`, `randmethod`, `blocksize`, `postprobmethod`, `multiarm_method`).
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
#' - Wald test and Wald confidence intervals may not be valid for BRAR data.
#' - For adaptive designs, prefer "post_mean" or "IPW" estimation and "randomization" or "AP" testing.
#'
#' @examples
#' # --- Example setup: simulate a BRAR trial with multiarm_method = "top2" ---
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
#'   observation_delay = 20
#' )
#'
#' # --- Example 1: Simulation-based test ---
#' analyze_brar_trial(
#'   trial_data = results_binary_multiarm,
#'   N = 150, arms = 3,
#'   direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow = 2),
#'   distribution = "bernoulli",
#'   estimation_method = "post_mean",
#'   test_method = "simulation", CI_method = "simulation",
#'   B = 2000
#' )
#'
#' # --- Example 2: Wald test (not recommended for BRAR) ---
#' analyze_brar_trial(
#'   trial_data = results_binary_multiarm,
#'   N = 150, arms = 3,
#'   direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow = 2),
#'   distribution = "bernoulli",
#'   estimation_method = "MLE",
#'   test_method = "wald", CI_method = "wald"
#' )
#'
#' # --- Example 3: AP test (Allocation-Probability test) ---
#' analyze_brar_trial(
#'   trial_data = results_binary_multiarm,
#'   N = 150, arms = 3,
#'   direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow = 2),
#'   distribution = "bernoulli",
#'   estimation_method = "IPW",
#'   test_method = "AP", CI_method = "simulation",
#'   randmethod = "coin", blocksize = 1,
#'   postprobmethod = "simulation", multiarm_method = "top2")
#'
#' @export
analyze_brar_trial <- function(
    trial_data, N, arms, direction, priors, distribution, known_var = NULL,
    estimation_method = c("MLE", "IPW", "post_mean"),
    test_method = c("wald", "exact", "randomization", "AP", "simulation"),
    CI_method = c("wald", "simulation"),
    effect_measure = "riskdifference",
    multiple_tests = FALSE, onesided = TRUE,
    alpha = 0.05, prob_threshold = NULL, ...
) {

  # --- Validate arguments ---
  distribution = match.arg(distribution, c("bernoulli", "normal", "exponential"))
  direction = match.arg(direction, c("lower", "higher"))
  estimation_method = match.arg(estimation_method, c("MLE", "IPW", "post_mean"))
  test_method = match.arg(test_method, c("wald", "exact", "randomization", "AP", "simulation"))
  CI_method = match.arg(CI_method, c("wald", "simulation"))

  if (distribution == "bernoulli") {
    effect_measure = match.arg(effect_measure, c("riskdifference"))
  }

  # --- Method warnings for BRAR data ---
  if (estimation_method == "MLE")
    warning("MLE estimation may be biased under BRAR. Consider 'IPW' or 'post_mean'.")
  if (test_method == "wald")
    warning("Wald test may not be appropriate for BRAR data. Consider 'randomization', 'simulation', or 'AP'.")
  if (CI_method == "wald")
    warning("Wald confidence intervals may not be appropriate for BRAR data. Consider 'simulation' CIs.")



  # --- Dispatch to correct helper based on distribution ---
  results = switch(
    distribution,
    "bernoulli" = .analyze_brar_trial_binary(
      trial_data = trial_data, priors = priors, N = N, arms = arms, direction = direction,
      estimation_method = estimation_method, test_method = test_method, CI_method = CI_method,
      effect_measure = effect_measure, alpha = alpha, multiple_tests = multiple_tests,
      onesided = onesided, ...
    ),
    "normal" = .analyze_brar_trial_normal(
      trial_data = trial_data, priors = priors, N = N, arms = arms, direction = direction,
      known_var = known_var,
      estimation_method = estimation_method, test_method = test_method,
      alpha = alpha, multiple_tests = multiple_tests, onesided = onesided, ...
    ),
    "exponential" = .analyze_brar_trial_exp(
      trial_data = trial_data, priors = priors, N = N, arms = arms, direction = direction,
      estimation_method = estimation_method, test_method = test_method,
      alpha = alpha, multiple_tests = multiple_tests, onesided = onesided, ...
    )
  )

  return(results)
}
