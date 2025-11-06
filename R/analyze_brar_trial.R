#' @title Analyze a Single BRAR/FR Trial Replicate
#'
#' @description
#' Calculates the operating characteristics (estimation, hypothesis test results,
#' patient benefit, and allocation metrics) for a single simulated clinical trial
#' replicate.
#'
#' This function serves as a wrapper that validates inputs and dispatches
#' to the appropriate internal helper function based on the outcome's
#' \code{distribution} (e.g., \code{.analyze_brar_trial_binary}).
#'
#' @param trial_data A data frame containing the results of a single trial simulation
#'                   (output of \code{simulate_brar_trial} or \code{simulate_fr_trial}).
#' @param N Numeric. Total sample size for the trial.
#' @param arms Numeric. Number of arms in the trial.
#' @param direction Character. Specifies the direction of the desired treatment effect
#'                  (\code{"lower"} or \code{"higher"}).
#' @param priors Matrix. The prior parameters used in the simulation (needed for Bayesian estimation).
#' @param modelpar Matrix or Numeric vector. The true parameters for each arm (needed for benefit and power/T1E calculations).
#' @param distribution Character. Specifies the outcome distribution (\code{"bernoulli"}, \code{"normal"}, or \code{"exponential"}).
#' @param known_var Logical. If \code{TRUE}, the variance is known for Normal outcomes.
#' @param estimation_method Character. Method to estimate the effect for each arm (e.g., "MLE", "posterior_mean").
#' @param test_method Character. Hypothesis testing method (e.g., "wald", "exact_conditional", "bayesian_pp").
#' @param alpha Numeric. Significance level (frequentist Type I Error rate) or 1-CI level.
#' @param prob_threshold Numeric. The posterior probability threshold for Bayesian testing.
#' @param null_effect Numeric. Value representing the null hypothesis (typically 0).
#' @param ... Additional arguments to be passed to internal test methods (e.g., \code{B} for permutation test).
#'
#' @return A list containing the calculated metrics for the single replicate,
#'         as returned by the internal helper functions.
#' @export
analyze_brar_trial <- function(
    trial_data, N, arms, direction, priors, modelpar, distribution, known_var,
    estimation_method, test_method, alpha, prob_threshold, null_effect, ...
) {

  # 1. Input Validation
  # Basic checks before dispatching
  if (missing(trial_data) || missing(distribution) || missing(modelpar)) {
    stop("Arguments 'trial_data', 'distribution', and 'modelpar' must be provided.")
  }

  distribution <- match.arg(distribution, c("bernoulli", "normal", "exponential"))

  # 2. Dispatch to internal helper
  # Call the appropriate function based on the distribution

  results <- switch(
    distribution,

    "bernoulli" = .analyze_brar_trial_binary(
      trial_data = trial_data, N = N, arms = arms, direction = direction,
      priors = priors, modelpar = modelpar,
      estimation_method = estimation_method, test_method = test_method,
      alpha = alpha, prob_threshold = prob_threshold, null_effect = null_effect, ...
    ),

    "normal" = .analyze_brar_trial_normal(
      trial_data = trial_data, N = N, arms = arms, direction = direction,
      priors = priors, modelpar = modelpar, known_var = known_var,
      estimation_method = estimation_method, test_method = test_method,
      alpha = alpha, prob_threshold = prob_threshold, null_effect = null_effect, ...
    ),

    "exponential" = .analyze_brar_trial_exp(
      trial_data = trial_data, N = N, arms = arms, direction = direction,
      priors = priors, modelpar = modelpar,
      estimation_method = estimation_method, test_method = test_method,
      alpha = alpha, prob_threshold = prob_threshold, null_effect = null_effect, ...
    ),

    # Default case if no match (should be caught by match.arg, but good practice)
    stop(paste("No analysis function available for distribution:", distribution))
  )

  # 3. Return results from the helper
  return(results)
}
