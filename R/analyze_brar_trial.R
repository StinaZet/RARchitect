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
#' @param ... Additional arguments for test methods (e.g., B, policy_file).
#'
#' @return List with \code{test_results_df} and \code{patient_benefit}.
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

  # Validate inputs
  distribution <- match.arg(distribution, c("bernoulli", "normal", "exponential"))
  direction <- match.arg(direction, c("lower", "higher"))
  estimation_method <- match.arg(estimation_method, c("MLE", "IPW", "post_mean"))
  test_method <- match.arg(test_method, c("wald", "exact", "randomization", "AP", "simulation"))
  CI_method <- match.arg(CI_method, c("wald", "simulation"))

  if (distribution == "bernoulli") {
    effect_measure <- match.arg(effect_measure, c("riskdifference"))
  }

  # Dispatch to internal helper
  results <- switch(
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
