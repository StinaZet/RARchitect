#' @title Analyze a Single Binary (Bernoulli) Trial Replicate (Internal)
#'
#' @description
#' Internal helper function to calculate operating characteristics for a
#' single simulated clinical trial replicate with binary outcomes.
#'
#' @param trial_data A data frame from \code{simulate_brar_trial}.
#' @param ... All other arguments passed from the main \code{analyze_brar_trial} wrapper.
#' @param B Numeric. Number of permutations for \code{test_method = "randomisation_based"}.
#'
#' @return A list containing the calculated metrics.
#' @keywords internal
.analyze_brar_trial_binary <- function(
    trial_data, N, arms, direction, priors, modelpar,
    estimation_method, test_method, alpha, prob_threshold, null_effect, B = 1000, ...
) {

  # 1. Data Processing and Summary Statistics
  # ----------------------------------------
  arm_counts <- table(trial_data$Arm)
  observed_outcomes <- tapply(trial_data$Outcome, trial_data$Arm, sum)

  arm_summary <- list()
  for (k in 1:arms) {
    n_k <- arm_counts[as.character(k)]
    y_k <- observed_outcomes[as.character(k)]
    if (is.na(n_k)) n_k = 0
    if (is.na(y_k)) y_k = 0
    arm_summary[[k]] <- list(n = n_k, y = y_k) # n = count, y = successes
  }

  # Determine the TRUE best arm
  true_parameters <- modelpar # modelpar is a vector of probabilities

  if (direction == "higher") {
    true_best_arm <- which.max(true_parameters)
    true_best_parameter_value <- max(true_parameters)
  } else {
    true_best_arm <- which.min(true_parameters)
    true_best_parameter_value <- min(true_parameters)
  }
  true_control_parameter_value <- true_parameters[1]

  # 2. Estimation
  # ---------------------------------
  arm_estimates <- numeric(arms)
  for (k in 1:arms) {
    n_k <- arm_summary[[k]]$n
    y_k <- arm_summary[[k]]$y

    if (n_k > 0) {
      if (estimation_method == "MLE") {
        arm_estimates[k] <- y_k / n_k # P-hat
      } else if (estimation_method == "posterior_mean") {
        # Beta-Binomial posterior mean
        arm_estimates[k] <- (priors[1, k] + y_k) / (priors[1, k] + priors[2, k] + n_k)
      } # Add other estimation methods (e.g., median) here
    } else {
      # Use prior mean if no patients are allocated
      arm_estimates[k] <- priors[1, k] / (priors[1, k] + priors[2, k])
    }
  }

  # Identify the estimated best *experimental* arm (Arm > 1)
  exp_arm_estimates <- arm_estimates[-1]
  if (direction == "higher") {
    estimated_best_exp_arm <- which.max(exp_arm_estimates) + 1
  } else {
    estimated_best_exp_arm <- which.min(exp_arm_estimates) + 1
  }
  effect_estimate <- arm_estimates[estimated_best_exp_arm] - arm_estimates[1]

  # 3. Hypothesis Testing
  # -----------------------------------
  test_outcome <- 0 # Default: Fail to reject H0
  p_value <- NA     # Store p-value if frequentist

  k_star <- estimated_best_exp_arm
  n_c <- arm_summary[[1]]$n # Control count
  y_c <- arm_summary[[1]]$y # Control successes
  n_e <- arm_summary[[k_star]]$n # Experimental count
  y_e <- arm_summary[[k_star]]$y # Experimental successes

  if (n_c > 0 && n_e > 0) {

    # --- Wald (Normal Approximation) Test ---
    if (test_method == "wald") {
      counts <- matrix(c(y_e, n_e - y_e, y_c, n_c - y_c), nrow = 2)
      wald_test <- try(prop.test(counts, correct = TRUE), silent = TRUE)

      if (!inherits(wald_test, "try-error")) {
        z_stat <- sign(y_e/n_e - y_c/n_c) * sqrt(wald_test$statistic)
        p_value <- if (direction == "higher") pnorm(z_stat, lower.tail = FALSE) else pnorm(z_stat, lower.tail = TRUE)
      }
    }

    # --- Exact (Fisher's Conditional) Test ---
    else if (test_method == "exact_conditional") {
      counts <- matrix(c(y_e, n_e - y_e, y_c, n_c - y_c), nrow = 2)
      alternative <- if (direction == "higher") "greater" else "less"
      fisher_test <- fisher.test(counts, alternative = alternative)
      p_value <- fisher_test$p.value
    }

    # --- Randomization-based Test ---
    else if (test_method == "randomisation_based") {
      p_value <- .perform_permutation_test(y_c, n_c, y_e, n_e, direction, B = B)
    }

    # --- AP Test (Adaptive Permutation Test) (Placeholder) ---
    else if (test_method == "ap_test") {
      # p_value = calculate_ap_p_value(trial_data, k_star, ...)
    }

    # --- Bayesian Posterior Probability Test ---
    else if (test_method == "bayesian_pp") {
      # posterior_prob = ...
      # if (posterior_prob > prob_threshold) test_outcome <- 1
    }

    # Set final frequentist test outcome
    if (!is.na(p_value) && p_value < alpha) {
      test_outcome <- 1
    }
  }

  # 4. Confidence Interval / Credible Interval Coverage
  # --------------------------------------------------
  coverage <- 0
  true_effect <- true_parameters[true_best_arm] - true_control_parameter_value
  # ... CI/CRI calculation logic goes here ...

  # 5. Patient Benefit Calculation
  # ------------------------------
  trial_mean_outcome <- sum(trial_data$Outcome) / N
  patient_benefit <- trial_mean_outcome

  # 6. Allocation Metrics
  # ---------------------
  best_arm_alloc_count <- sum(trial_data$Arm == true_best_arm)

  return(list(
    effect_estimate = effect_estimate,
    test_outcome = test_outcome,
    p_value = p_value,
    patient_benefit = patient_benefit,
    best_arm_alloc_count = best_arm_alloc_count,
    coverage = coverage
  ))
}
