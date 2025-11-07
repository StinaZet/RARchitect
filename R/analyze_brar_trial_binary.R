#' @title Analyze a Single Binary (Bernoulli) Trial Replicate (Internal)
#'
#' @description
#' Internal helper function to calculate operating characteristics for a
#' single simulated clinical trial replicate with binary outcomes.
#'
#' @param trial_data A data frame from \code{simulate_brar_trial}.
#' @param effect_measure Character. The measure of effect. Currently only
#'   \code{"riskdifference"} is supported.
#' @param multiple_tests Logical. If \code{FALSE} (default), only test the
#'   *estimated* best experimental arm vs. control. If \code{TRUE}, test
#'   *all* experimental arms vs. control.
#' @param onesided Logical. If \code{TRUE} (default), performs a one-sided test
#'   based on \code{direction}. If \code{FALSE}, performs a two-sided test.
#' @param ... All other arguments passed from the main \code{analyze_brar_trial} wrapper.
#' @param B Numeric. Number of permutations for \code{test_method = "randomisation_based"}.
#'
#' @return A list containing the calculated metrics.
#' @keywords internal
.analyze_brar_trial_binary <- function(
    trial_data, N, arms, direction, priors, modelpar,
    estimation_method = c("MLE", "IPW", "post_mean"),
    test_method = c("wald", "exact", "randomization", "AP", "simulation"),
    CI_method = c("wald", "simulation"),
    effect_measure = "riskdifference",
    alpha, prob_threshold, null_effect,
    multiple_tests = FALSE, onesided = TRUE, B = 10000, ...
) {

  # 1. Data Processing and Summary Statistics
  # ----------------------------------------
  arm_counts = table(trial_data$Arm)
  observed_outcomes = tapply(trial_data$Outcome, trial_data$Arm, sum)

  arm_summary = list()
  for (k in 1:arms) {
    n_k = arm_counts[as.character(k)]
    y_k = observed_outcomes[as.character(k)]
    if (is.na(n_k)) n_k = 0
    if (is.na(y_k)) y_k = 0
    arm_summary[[k]] = list(n = n_k, y = y_k) # n = count, y = successes
  }

  # Determine the TRUE best arm (for allocation/benefit metrics)
  true_parameters = modelpar # modelpar is a vector of probabilities

  if (direction == "higher") {
    true_best_arm = which.max(true_parameters)
  } else {
    true_best_arm = which.min(true_parameters)
  }
  true_control_parameter_value = true_parameters[1]

  # 2. Estimation
  # ---------------------------------
  arm_estimates = numeric(arms)

  for (k in 1:arms) {
    n_k = arm_summary[[k]]$n
    y_k = arm_summary[[k]]$y

    if (n_k > 0) {
      if (estimation_method == "MLE") {
        arm_estimates[k] = y_k / n_k
      } else if (estimation_method == "IPW")
      {
        # --- Normalized IPW (Hajek) Estimator ---

        # --- Calculate weights for IPW ---
        # Get allocation probabilities from trial_data
        ap_matrix = as.matrix(trial_data[, paste0("AP arm ", 1:arms)])
        # Get the *actual* probability used to assign the patient
        pi_obs = ap_matrix[cbind(1:N, trial_data$Arm)]
        # Calculate inverse probability weights
        weights = 1 / pi_obs
        # --- END IPW WEIGHTS ---

        # Get indices for all patients assigned to arm k
        arm_k_indices = which(trial_data$Arm == k)

        # Numerator: Sum of weighted outcomes
        # (Outcome_i * Weight_i) for all patients in arm k
        numerator = sum(trial_data$Outcome[arm_k_indices] * weights[arm_k_indices])

        # Denominator: Sum of weights
        # (Weight_i) for all patients in arm k
        denominator = sum(weights[arm_k_indices])

        if (denominator > 0) {
          arm_estimates[k] = numerator / denominator
        } else {
          arm_estimates[k] = NA # Arm was never assigned, denominator is 0
        }
      } else if (estimation_method == "post_mean") {
        # Beta-Binomial posterior mean
        arm_estimates[k] = (priors[1, k] + y_k) / (priors[1, k] + priors[2, k] + n_k)
      } # Add other estimation methods (e.g., median) here
    } else {
      # Use prior mean if no patients are allocated
      arm_estimates[k] = priors[1, k] / (priors[1, k] + priors[2, k])
    }
  }


  # 2.5. Define which arms to test based on 'multiple_tests'
  # --------------------------------------------------------
  if (multiple_tests == TRUE) {
    # Test all experimental arms (2 to K) against control (Arm 1)
    arms_to_test = 2:arms # Vector of arm indices: 2, 3, ...
  } else {
    # Test only the *estimated* best experimental arm
    if (arms == 2) {
      estimated_best_exp_arm = 2 # Only one experimental arm
    } else {
      exp_arm_estimates = arm_estimates[-1] # Estimates for arms 2, 3, ...
      if (direction == "higher") {
        estimated_best_exp_arm = which.max(exp_arm_estimates) + 1
      } else {
        estimated_best_exp_arm = which.min(exp_arm_estimates) + 1
      }
    }
    arms_to_test = c(estimated_best_exp_arm) # Vector with one element
  }

  # 3. & 4. Hypothesis Testing and CI Coverage (Loop)
  # --------------------------------------------------

  # Initialize a data frame to store results for each comparison
  n_comparisons = length(arms_to_test)
  test_results_df = data.frame(
    ExperimentalArm = arms_to_test,
    effect_estimate = numeric(n_comparisons),
    p_value = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    test_outcome = 0, # Default: Fail to reject H0
    coverage = 0       # Default: No coverage
  )

  n_c = arm_summary[[1]]$n # Control count (same for all comparisons)
  y_c = arm_summary[[1]]$y # Control successes (same for all comparisons)

  # Define the alternative hypothesis string
  if (onesided) {
    alternative_str = if (direction == "higher") "greater" else "less"
    conf_level = 1 - alpha # e.g. alpha=0.05 -> conf.level=0.95 (one-sided)
  } else {
    alternative_str = "two.sided"
    conf_level = 1 - alpha # e.g. alpha=0.05 -> conf.level=0.95 (two-sided)
  }

  # --- Exact Test Validation --- (from previous update)
  if (test_method == "exact") {
    # ... (logic remains for validation) ...
    policy_file = list(...)$policy_file # Retrieve if available
    if (arms > 2) {
      stop("The 'exact' (CXS) test method is only supported for arms = 2.")
    }
    if (multiple_tests) {
      stop("The 'exact' (CXS) test method cannot be combined with 'multiple_tests = TRUE'.")
    }
    if (is.null(policy_file)) {
      stop("For 'test_method = \"exact\"', 'policy_file' must be specified.")
    }
  }

  # --- Load policy file ONCE if needed ---
  policy_data = NULL
  if (test_method == "exact") {
    # Use blocksize from parameters
    blocksize = list(...)$blocksize
    # This function loads from cache or disk
    policy_data = .load_cxs_policy(policy_file, N, blocksize)
  }
  # --- End Load ---


  for (i in 1:n_comparisons) {

    k_star = test_results_df$ExperimentalArm[i] # Current experimental arm to test
    n_e = arm_summary[[k_star]]$n # Experimental count
    y_e = arm_summary[[k_star]]$y # Experimental successes

    # --- Loop over effect_measure ---
    if (effect_measure == "riskdifference") {

      # --- Estimation (Risk Difference) ---
      test_results_df$effect_estimate[i] = arm_estimates[k_star] - arm_estimates[1]
      true_effect = true_parameters[k_star] - true_control_parameter_value

      # --- Hypothesis Testing (Risk Difference) ---
      p_value_temp = NA_real_
      ci_temp = c(NA_real_, NA_real_)

      # Run tests and CIs only if we have data in both arms
      if (n_c > 0 && n_e > 0)
      {
        counts = matrix(c(y_e, n_e - y_e, y_c, n_c - y_c), nrow = 2)

        # --- A. Test Method Switch ---
        if (test_method == "wald")
        {
          wald_test = try(stats::prop.test(counts, correct = TRUE), silent = TRUE)
          if (!inherits(wald_test, "try-error"))
          {
            z_stat = sign(y_e/n_e - y_c/n_c) * sqrt(wald_test$statistic)
            if (onesided)
            {
              p_value_temp = if (direction == "higher") stats::pnorm(z_stat, lower.tail = FALSE) else stats::pnorm(z_stat, lower.tail = TRUE)
            } else {
              p_value_temp = 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
            }
          }
        }
        else if (test_method == "exact")
        {
          p_value_temp = .perform_cxs_test(y_c = y_c, n_c = n_c, y_e = y_e, n_e = n_e, N = N, pvals_vector = policy_data$pvals_CXS)
          thr_bt = policy_data$thr_exact_BT
          if (!is.na(p_value_temp) && !is.na(thr_bt)) {
            if (p_value_temp <= thr_bt) {
              test_results_df$test_outcome[i] = 1
            }
          }
        }
        else if (test_method == "randomization")
        {
          p_value_temp = .perform_permutation_test(y_c, n_c, y_e, n_e,
                                                   alternative = alternative_str, B = B)
        }
        else if (test_method == "ap_test") {
          # p_value_temp = ...
        }
        # --- Monte Carlo (Parametric Bootstrap) Test ---
        else if (test_method == "simulation")
        {
          p_value_temp = perform_mc_test(y_c, n_c, y_e, n_e,
                                           alternative = alternative_str, B = B)
        }

        # --- B. CI Method Switch ---
        if (CI_method == "wald") {
          # Standard Wald/Normal CI (based on prop.test's method)
          prop_ci = try(stats::prop.test(counts, conf.level = conf_level)$conf.int, silent = TRUE)
          if (!inherits(prop_ci, "try-error")) {
            ci_temp = as.numeric(prop_ci)
          }
        } else if (CI_method == "simulation") {
          # Percentile Bootstrap CI
          ci_temp = perform_bootstrap_ci(y_c, n_c, y_e, n_e,
                                           alpha = alpha,
                                           onesided = onesided,
                                           direction = direction,
                                           B = B)
        }
        # --- End CI Method Switch ---


        # Set final frequentist test outcome (unless exact, handled above)
        if (test_method != "exact" && !is.na(p_value_temp)) {
          test_results_df$p_value[i] = p_value_temp
          if (p_value_temp < alpha) {
            test_results_df$test_outcome[i] = 1
          }
        }

        # Store CI results
        test_results_df$ci_low[i] = ci_temp[1]
        test_results_df$ci_high[i] = ci_temp[2]

        # --- Coverage (Risk Difference) ---
        # Check if the true effect lies within the calculated CI
        if (!is.na(ci_temp[1]) && !is.na(ci_temp[2])) {
          if (true_effect >= ci_temp[1] && true_effect <= ci_temp[2]) {
            test_results_df$coverage[i] = 1
          }
        }

      } # end if (n_c > 0 && n_e > 0)

    } else if (effect_measure == "riskratio") {

      # --- Placeholder for Risk Ratio ---
      # test_results_df$effect_estimate[i] = arm_estimates[k_star] / arm_estimates[1]
      # ... (Test logic for RR, e.g., log-transform, different variance)
      # ... (CI logic for RR)

    } else if (effect_measure == "oddsratio") {

      # --- Placeholder for Odds Ratio ---
      # p_e = arm_estimates[k_star]; p_c = arm_estimates[1]
      # test_results_df$effect_estimate[i] = (p_e/(1-p_e)) / (p_c/(1-p_c))
      # ... (Test logic for OR, e.g., logistic regression, Mantel-Haenszel)
      # ... (CI logic for OR)
    }
  } # --- End loop over arms_to_test ---

  # 5. Patient Benefit Calculation
  # ------------------------------
  trial_mean_outcome = sum(trial_data$Outcome) / N
  patient_benefit = trial_mean_outcome

  best_arm_alloc_count = sum(trial_data$Arm == true_best_arm)

  return(list(
    test_results_df = test_results_df, # Data frame of all test results
    patient_benefit = patient_benefit,
    best_arm_alloc_count = best_arm_alloc_count
  ))
}
