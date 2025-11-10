#' @title Analyze a Single BRAR/FR Trial Replicate
#'
#' @description
#' Wrapper to analyze a trial replicate with real data. Assumes binary, normal, or exponential outcomes.
#' Reports effect estimates, p-values, and confidence intervals (no coverage).
#' Patient benefit = number of patients on estimated best arm.
#'
#' @param trial_data Data frame containing trial results.
#' @param N Total sample size.
#' @param arms Number of arms.
#' @param direction "lower" or "higher".
#' @param priors Matrix of prior parameters for Bayesian estimation.
#' @param distribution "bernoulli", "normal", or "exponential".
#' @param known_var Logical, for normal outcomes.
#' @param estimation_method "MLE", "IPW", or "post_mean".
#' @param test_method "wald", "exact", "randomization", "AP", or "simulation".
#' @param CI_method "wald" or "simulation".
#' @param effect_measure Currently only "riskdifference" for binary outcomes.
#' @param multiple_tests Logical, test all experimental arms vs control with Bonferroni adjustment.
#' @param onesided Logical, one-sided if TRUE.
#' @param alpha Significance level (used for CI).
#' @param prob_threshold Posterior probability threshold (future Bayesian use).
#' @param ... Additional arguments for test methods (e.g., B, policy_file, randmethod, blocksize).
#'
#' @return List with \code{test_results_df} and \code{patient_benefit}.
#' @export
#'
#' @examples
#' # Simulation-based test example
#' analyze_brar_trial(
#'   trial_data = trial_df, N = 100, arms = 3, direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow=2), distribution = "bernoulli",
#'   test_method = "simulation", CI_method = "simulation"
#' )
#'
#' # Wald test example (will give a warning)
#' analyze_brar_trial(
#'   trial_data = trial_df, N = 100, arms = 3, direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow=2), distribution = "bernoulli",
#'   test_method = "wald", CI_method = "wald"
#' )
#'
#' # AP test example
#' analyze_brar_trial(
#'   trial_data = trial_df, N = 100, arms = 3, direction = "higher",
#'   priors = matrix(c(1,1,1,1,1,1), nrow=2), distribution = "bernoulli",
#'   test_method = "AP", CI_method = "simulation",
#'   randmethod = "block", blocksize = 5, postprobmethod = "simulation"
#' )
analyze_brar_trial <- function(trial_data, N, arms, direction, priors, distribution,
                               estimation_method=c("MLE","post_mean"),
                               test_method=c("wald","randomization","simulation","AP"),
                               CI_method=c("wald","simulation"),
                               multiple_tests=FALSE, onesided=TRUE, alpha=0.05, ...) {

  # Validate inputs
  distribution = match.arg(distribution, c("bernoulli", "normal", "exponential"))
  estimation_method = match.arg(estimation_method)
  test_method = match.arg(test_method)
  CI_method = match.arg(CI_method)
  direction = match.arg(direction, c("higher","lower"))

  if (distribution == "bernoulli") {
    effect_measure = match.arg(effect_measure, c("riskdifference"))
  }

  # Warnings about methods not suitable for BRAR data
  if (estimation_method == "MLE") {
    warning("MLE estimation may not be appropriate for BRAR data. Consider 'IPW' or 'post_mean'.")
  }
  if (test_method == "wald") {
    warning("Wald test may not be appropriate for BRAR data. Consider 'randomization', 'simulation', or 'AP'.")
  }
  if (CI_method == "wald") {
    warning("Wald confidence intervals may not be appropriate for BRAR data. Consider 'simulation' CIs.")
  }

  # Dispatch to internal helper
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


