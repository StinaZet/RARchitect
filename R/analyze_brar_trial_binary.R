#' @title Analyze a Single Binary (Bernoulli) BRAR Trial Replicate (Internal)
#'
#' @description
#' Internal helper function to calculate inference for a single trial replicate
#' with binary outcomes.
#'
#' @param trial_data A data frame with trial results (must include columns
#'   'Arm' and 'Outcome'; optionally 'AP arm 1', 'AP arm 2', ... if using IPW).
#' @param priors Matrix of prior parameters for Bayesian estimation.
#' @param N Total number of patients.
#' @param arms Number of treatment arms.
#' @param direction "higher" or "lower".
#' @param estimation_method "MLE", "IPW", or "post_mean".
#' @param test_method "wald", "exact", "randomization", "AP", or "simulation".
#' @param CI_method "wald" or "simulation".
#' @param effect_measure Currently only "riskdifference".
#' @param alpha Significance level for CI calculation.
#' @param multiple_tests Logical; if TRUE, all experimental arms vs control.
#' @param onesided Logical; if TRUE, one-sided test.
#' @param B Number of simulations/permutations for simulation/randomization-based tests.
#' @param ... Additional arguments passed to test methods.
#'
#' @return List with:
#'   - \code{test_results_df}: Data frame with effect estimates, p-values, and CIs.
#'   - \code{patient_benefit}: Number of patients on estimated best experimental arm.
#'   - \code{mean_outcome}: Mean outcome across all patients.
#' @keywords internal
.analyze_brar_trial_binary <- function(
    trial_data, priors, N, arms, direction,
    estimation_method = c("MLE", "IPW", "post_mean"),
    test_method = c("wald", "exact", "randomization", "AP", "simulation"),
    CI_method = c("wald", "simulation"),
    effect_measure = "riskdifference",
    alpha = 0.05, multiple_tests = FALSE, onesided = TRUE, B = 10000, ...
) {

  # --- Argument matching ---
  estimation_method <- match.arg(estimation_method)
  test_method <- match.arg(test_method)
  CI_method <- match.arg(CI_method)
  effect_measure <- match.arg(effect_measure, c("riskdifference"))
  direction <- match.arg(direction, c("higher", "lower"))

  alpha_adj <- if (multiple_tests) alpha / (arms - 1) else alpha

  # --- Summarize data ---
  arm_counts <- table(trial_data$Arm)
  observed_outcomes <- tapply(trial_data$Outcome, trial_data$Arm, sum)
  arm_summary <- vector("list", arms)
  for (k in 1:arms) {
    n_k <- ifelse(is.na(arm_counts[as.character(k)]), 0, arm_counts[as.character(k)])
    y_k <- ifelse(is.na(observed_outcomes[as.character(k)]), 0, observed_outcomes[as.character(k)])
    arm_summary[[k]] <- list(n = n_k, y = y_k)
  }

  # --- Estimation ---
  arm_estimates <- numeric(arms)
  for (k in 1:arms) {
    n_k <- arm_summary[[k]]$n
    y_k <- arm_summary[[k]]$y
    if (n_k > 0) {
      if (estimation_method == "MLE") {
        arm_estimates[k] <- y_k / n_k
      } else if (estimation_method == "IPW") {
        ap_cols <- grep("^AP arm ", names(trial_data))
        ap_matrix <- as.matrix(trial_data[, ap_cols, drop = FALSE])
        pi_obs <- ap_matrix[cbind(seq_len(nrow(trial_data)), trial_data$Arm)]
        weights <- 1 / pi_obs
        arm_k_indices <- which(trial_data$Arm == k)
        numerator <- sum(trial_data$Outcome[arm_k_indices] * weights[arm_k_indices])
        denominator <- sum(weights[arm_k_indices])
        arm_estimates[k] <- ifelse(denominator > 0, numerator / denominator, NA)
      } else if (estimation_method == "post_mean") {
        arm_estimates[k] <- (priors[1, k] + y_k) / (priors[1, k] + priors[2, k] + n_k)
      }
    } else {
      arm_estimates[k] <- priors[1, k] / (priors[1, k] + priors[2, k])
    }
  }

  # --- Determine arms to test ---
  if (multiple_tests) {
    arms_to_test <- 2:arms
  } else {
    exp_arm_estimates <- arm_estimates[-1]
    estimated_best_exp_arm <- if (direction == "higher") which.max(exp_arm_estimates) + 1 else which.min(exp_arm_estimates) + 1
    arms_to_test <- c(estimated_best_exp_arm)
  }

  # --- Prepare result dataframe ---
  test_results_df <- data.frame(
    ExperimentalArm = arms_to_test,
    effect_estimate = arm_estimates[arms_to_test] - arm_estimates[1],
    p_value = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_
  )

  n_c <- arm_summary[[1]]$n
  y_c <- arm_summary[[1]]$y

  alternative_str <- if (onesided) {
    if (direction == "higher") "greater" else "less"
  } else "two.sided"

  conf_level <- if (onesided) 1 - alpha_adj else 1 - alpha_adj / 2

  # --- Simulation-based tests for all arms at once (AP, simulation, randomization) ---
  if (test_method %in% c("AP", "simulation", "randomization")) {
    # Collect parameters for null simulation
    sim_args <- list(
      trial_data = trial_data, priors = priors, arms = arms,
      direction = direction, B = B, multiple_tests = multiple_tests, ...
    )

    if (test_method == "AP") {
      test_results_df$p_value <- ap_test_brar(
        n_c = n_c, y_c = y_c, arm_summary = arm_summary,
        arms_to_test = arms_to_test, direction = direction,
        B = B, trial_data = trial_data, priors = priors, ...
      )
    } else if (test_method == "simulation") {
      test_results_df$p_value <- simulation_test_brar(
        n_c = n_c, y_c = y_c, arm_summary = arm_summary,
        arms_to_test = arms_to_test, direction = direction,
        B = B, trial_data = trial_data, priors = priors, ...
      )
    } else if (test_method == "randomization") {
      test_results_df$p_value <- randomization_test_brar(
        n_c = n_c, y_c = y_c, arm_summary = arm_summary,
        arms_to_test = arms_to_test, direction = direction,
        B = B, trial_data = trial_data, priors = priors, ...
      )
    }

  } else {
    # Wald and exact tests: loop over each arm individually
    for (i in seq_along(arms_to_test)) {
      k_star <- arms_to_test[i]
      n_e <- arm_summary[[k_star]]$n
      y_e <- arm_summary[[k_star]]$y

      if (n_c == 0 || n_e == 0) {
        test_results_df[i, c("effect_estimate", "p_value", "ci_low", "ci_high")] <- NA
        next
      }

      # Wald test
      if (test_method == "wald") {
        p1 <- y_e / n_e; p0 <- y_c / n_c
        se <- sqrt(p1 * (1 - p1) / n_e + p0 * (1 - p0) / n_c)
        z_stat <- (p1 - p0) / se
        if (onesided) {
          test_results_df$p_value[i] <- if (direction == "higher") stats::pnorm(z_stat, lower.tail = FALSE) else stats::pnorm(z_stat, lower.tail = TRUE)
        } else {
          test_results_df$p_value[i] <- 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
        }
      }

      # Exact test
      if (test_method == "exact") {
        policy_file <- list(...)$policy_file
        blocksize <- list(...)$blocksize
        policy_data <- .load_cxs_policy(policy_file, N, blocksize)
        test_results_df$p_value[i] <- .perform_cxs_test(y_c, n_c, y_e, n_e, N, policy_data$pvals_CXS)
      }
    }
  }

  # --- Confidence intervals ---
  for (i in seq_along(arms_to_test)) {
    k_star <- arms_to_test[i]
    n_e <- arm_summary[[k_star]]$n
    y_e <- arm_summary[[k_star]]$y
    if (n_c == 0 || n_e == 0) next

    if (CI_method == "wald") {
      diff_est <- y_e / n_e - y_c / n_c
      se <- sqrt((y_e / n_e * (1 - y_e / n_e) / n_e) + (y_c / n_c * (1 - y_c / n_c) / n_c))
      z_crit <- stats::qnorm(conf_level)
      test_results_df[i, c("ci_low", "ci_high")] <- diff_est + c(-1, 1) * z_crit * se
    } else if (CI_method == "simulation") {
      test_results_df[i, c("ci_low", "ci_high")] <- perform_bootstrap_ci(
        y_c, n_c, y_e, n_e, alpha = alpha_adj, onesided = onesided, direction = direction, B = B
      )
    }
  }

  # --- Patient benefit ---
  best_arm <- if (multiple_tests) {
    exp_arm_estimates <- arm_estimates[-1]
    if (direction == "higher") which.max(exp_arm_estimates) + 1 else which.min(exp_arm_estimates) + 1
  } else {
    arms_to_test[1]
  }
  patient_benefit <- sum(trial_data$Arm == best_arm)
  mean_outcome <- mean(trial_data$Outcome)

  return(list(
    test_results_df = test_results_df,
    patient_benefit = patient_benefit,
    mean_outcome = mean_outcome
  ))
}
