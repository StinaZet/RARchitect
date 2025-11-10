#' @title Internal Analysis of a Single Binary BRAR Trial
#'
#' @description
#' Internal helper function to analyze a single BRAR trial replicate with binary outcomes.
#' Computes effect estimates, p-values, and confidence intervals for control vs. experimental arms.
#' Supports multiple experimental arms, one-sided or two-sided testing, and several test types
#' including Wald, randomization (permutation), Monte Carlo simulation, and AP test.
#'
#' @param trial_data Data frame with at least columns:
#'   - \code{Arm}: integer indicating treatment arm (1 = control, 2+ = experimental).
#'   - \code{Outcome}: binary outcome (0 or 1).
#' @param priors 2 x K numeric matrix of prior parameters (alpha/beta) for each arm.
#' @param N Integer. Total number of patients in the trial.
#' @param arms Integer. Number of treatment arms.
#' @param direction Character. Either "higher" or "lower" indicating which direction is considered better.
#' @param estimation_method Character. Method for estimating arm probabilities: "MLE" or "post_mean".
#' @param test_method Character. Test method: "wald", "randomization", "simulation", or "AP".
#' @param CI_method Character. Confidence interval method: "wald" or "simulation".
#' @param multiple_tests Logical. If TRUE, tests all experimental arms vs control using Bonferroni adjustment.
#' @param onesided Logical. If TRUE, performs a one-sided test in the specified direction.
#' @param alpha Numeric. Significance level for confidence intervals.
#' @param B Integer. Number of simulations/permutations for randomization, simulation, or AP test.
#' @param blocksize Integer. Block size for BRAR randomization sequences (used in AP/randomization tests).
#' @param ... Additional arguments passed to test methods (not used directly).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{test_results_df}}{Data frame with columns:
#'     - \code{ExperimentalArm}: tested experimental arm number.
#'     - \code{effect_estimate}: estimated risk difference vs control.
#'     - \code{p_value}: p-value from the chosen test.
#'     - \code{ci_low}: lower confidence interval bound.
#'     - \code{ci_high}: upper confidence interval bound.}
#'   \item{\code{patient_benefit}}{Number of patients assigned to the estimated best experimental arm.}
#'   \item{\code{mean_outcome}}{Mean outcome across all patients in the trial.}
#' }
#'
#' @details
#' This function is intended for internal use and is called by \code{analyze_brar_trial}.
#' It handles estimation for control and multiple experimental arms, calculates effect sizes,
#' and applies the specified statistical test and confidence interval method. For simulation-based
#' methods (AP or Monte Carlo), null distributions are generated according to the chosen method.
#'


.analyze_brar_trial_binary <- function(trial_data, priors, N, arms, direction,
                                       estimation_method="MLE", test_method="wald",
                                       CI_method="wald", multiple_tests=FALSE, onesided=TRUE,
                                       alpha=0.05, B=10000, blocksize=1, ...) {
  # Data summary
  arm_counts = table(trial_data$Arm)
  y_counts = tapply(trial_data$Outcome, trial_data$Arm, sum)
  arm_summary = lapply(1:arms, function(k) list(n=arm_counts[as.character(k)] %||% 0,
                                                 y=y_counts[as.character(k)] %||% 0))

  # Estimation
  arm_estimates = sapply(1:arms, function(k){
    n_k = arm_summary[[k]]$n
    y_k = arm_summary[[k]]$y
    if (n_k == 0) return(priors[1,k]/sum(priors[,k]))
    if (estimation_method=="MLE") return(y_k/n_k)
    if (estimation_method=="post_mean") return((priors[1,k]+y_k)/(sum(priors[,k])+n_k))
    stop("Only MLE or post_mean implemented")
  })

  # Determine arms to test
  if (multiple_tests) {
    arms_to_test = 2:arms
  } else {
    best_arm = if(direction=="higher") which.max(arm_estimates[-1])+1 else which.min(arm_estimates[-1])+1
    arms_to_test = best_arm
  }

  # Prepare result
  test_results_df = data.frame(
    ExperimentalArm = arms_to_test,
    effect_estimate = arm_estimates[arms_to_test]-arm_estimates[1],
    p_value = NA_real_, ci_low=NA_real_, ci_high=NA_real_
  )

  # Apply test
  n_c = arm_summary[[1]]$n
  y_c = arm_summary[[1]]$y
  y_e_list = lapply(arms_to_test, function(k) arm_summary[[k]]$y)
  n_e_list = lapply(arms_to_test, function(k) arm_summary[[k]]$n)

  if(test_method=="randomization") {
    test_results_df$p_value = randomization_test_brar(trial_data, priors, blocksize, arms_to_test)
  } else if(test_method=="simulation") {
    test_results_df$p_value = monte_carlo_test_brar(y_c, n_c, y_e_list, n_e_list)
  } else if(test_method=="AP") {
    test_results_df$p_value = ap_test_brar(trial_data, priors, blocksize)
  } else if(test_method=="wald") {
    for (i in seq_along(arms_to_test)) {
      p1 = y_e_list[[i]]/n_e_list[[i]]; p0 = y_c/n_c
      se = sqrt(p1*(1-p1)/n_e_list[[i]] + p0*(1-p0)/n_c)
      z = (p1-p0)/se
      test_results_df$p_value[i] = if(onesided) stats::pnorm(z, lower.tail = direction=="lower") else 2*stats::pnorm(abs(z), lower.tail=FALSE)
    }
  }

  # Confidence intervals
  for(i in seq_along(arms_to_test)) {
    if(CI_method=="wald") {
      diff = arm_estimates[arms_to_test[i]] - arm_estimates[1]
      se = sqrt((arm_estimates[arms_to_test[i]]*(1-arm_estimates[arms_to_test[i]])/n_e_list[[i]]) +
                   (arm_estimates[1]*(1-arm_estimates[1])/n_c))
      z = stats::qnorm(1-alpha)
      test_results_df[i,c("ci_low","ci_high")] = diff + c(-1,1)*z*se
    } else if(CI_method=="simulation") {
      test_results_df[i,c("ci_low","ci_high")] = perform_bootstrap_ci(y_c, n_c, y_e_list[[i]], n_e_list[[i]], alpha, onesided, direction)
    }
  }

  patient_benefit = sum(trial_data$Arm == arms_to_test[1])
  mean_outcome = mean(trial_data$Outcome)

  list(test_results_df=test_results_df, patient_benefit=patient_benefit, mean_outcome=mean_outcome)
}
