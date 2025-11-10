#' @title Analyze a Single BRAR/FR Trial Replicate
#'
#' @description
#' Wrapper function to analyze a single BRAR or Fixed Randomization trial replicate
#' with binary, normal, or exponential outcomes.
#' Returns effect estimates, p-values, and confidence intervals (no coverage).
#' Also reports patient benefit = number of patients on estimated best arm.
#'
#' @param trial_data Data frame with trial results.
#'   Must include columns 'Arm' and 'Outcome'.
#' @param N Total sample size.
#' @param arms Number of treatment arms.
#' @param direction "lower" or "higher" (indicates which outcome is favorable).
#' @param priors Matrix of prior parameters for Bayesian estimation (2 x K for binary).
#' @param distribution "bernoulli", "normal", or "exponential".
#' @param known_var Logical, used only for normal outcomes with known variance.
#' @param estimation_method "MLE", "IPW", or "post_mean" (for binary outcomes).
#' @param test_method "wald", "exact", "randomization", "AP", or "simulation".
#' @param CI_method "wald" or "simulation".
#' @param effect_measure Currently only "riskdifference" supported for binary outcomes.
#' @param multiple_tests Logical; if TRUE, test all experimental arms vs control with Bonferroni adjustment.
#' @param onesided Logical; if TRUE, perform one-sided test.
#' @param alpha Significance level (used for CI and tests).
#' @param prob_threshold Posterior probability threshold (for future Bayesian use).
#' @param ... Additional arguments passed to test methods (e.g., B, policy_file, randmethod, blocksize, postprobmethod, multiarm_method).
#'
#' @return List with components:
#'   - test_results_df: data frame with effect estimates, p-values, and CIs.
#'   - patient_benefit: number of patients on estimated best experimental arm.
#'   - mean_outcome: mean outcome across all patients.
#'
#' @details
#' Warnings are issued if methods may not be appropriate for BRAR data:
#'   - MLE estimation may be biased under adaptive randomization.
#'   - Wald test and Wald confidence intervals may not be valid under BRAR.
#'
#' @examples
#' # Simulation-based test
#' analyze_brar_trial(
#'   trial_data = trial_df, N = 100, arms = 3, direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow = 2),
#'   distribution = "bernoulli",
#'   test_method = "simulation", CI_method = "simulation",
#'   B = 5000
#' )
#'
#' # Wald test (with warnings)
#' analyze_brar_trial(
#'   trial_data = trial_df, N = 100, arms = 3, direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow = 2),
#'   distribution = "bernoulli",
#'   test_method = "wald", CI_method = "wald"
#' )
#'
#' # AP test example
#' analyze_brar_trial(
#'   trial_data = trial_df, N = 100, arms = 3, direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow = 2),
#'   distribution = "bernoulli",
#'   test_method = "AP", CI_method = "simulation",
#'   randmethod = "block", blocksize = 5, postprobmethod = "simulation"
#' )
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

  # --- Validate inputs ---
  distribution = match.arg(distribution, c("bernoulli", "normal", "exponential"))
  direction = match.arg(direction, c("lower", "higher"))
  estimation_method = match.arg(estimation_method, c("MLE", "IPW", "post_mean"))
  test_method = match.arg(test_method, c("wald", "exact", "randomization", "AP", "simulation"))
  CI_method = match.arg(CI_method, c("wald", "simulation"))

  if (distribution == "bernoulli") {
    effect_measure = match.arg(effect_measure, c("riskdifference"))
  }

  # --- Issue warnings for BRAR-specific considerations ---
  if (estimation_method == "MLE") {
    warning("MLE estimation may not be appropriate for BRAR data. Consider 'IPW' or 'post_mean'.")
  }
  if (test_method == "wald") {
    warning("Wald test may not be appropriate for BRAR data. Consider 'randomization', 'simulation', or 'AP'.")
  }
  if (CI_method == "wald") {
    warning("Wald confidence intervals may not be appropriate for BRAR data. Consider 'simulation' CIs.")
  }

  # --- Dispatch to internal helper functions ---
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
