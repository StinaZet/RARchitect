# Calculate a simulation-based confidence interval based on Monte Carlo.
perform_bootstrap_ci <- function(y_c, n_c, y_e, n_e,
                                  alpha, onesided, direction, B = 10000) {

  # 1. Calculate observed proportions
  p_e_obs = y_e / n_e
  p_c_obs = y_c / n_c

  # 2. Bootstrap simulation (Parametric Bootstrap from observed)
  # Simulate the *number of successes* from the *observed* proportions
  y_e_sim = stats::rbinom(B, n_e, p_e_obs)
  y_c_sim = stats::rbinom(B, n_c, p_c_obs)

  # 3. Calculate simulated statistics
  T_sim = (y_e_sim / n_e) - (y_c_sim / n_c)

  # 4. Compute CI bounds based on percentiles
  if (onesided) {
    if (direction == "higher") {
      # We want a (1-alpha) lower bound, e.g., 5th percentile
      ci_low = stats::quantile(T_sim, alpha, na.rm = TRUE)
      ci_high = Inf
    } else { # "less"
      # We want a (1-alpha) upper bound, e.g., 95th percentile
      ci_low = -Inf
      ci_high = stats::quantile(T_sim, 1 - alpha, na.rm = TRUE)
    }
  } else { # two.sided
    # Standard two-sided (1-alpha) percentile CI
    ci_low = stats::quantile(T_sim, alpha / 2, na.rm = TRUE)
    ci_high = stats::quantile(T_sim, 1 - (alpha / 2), na.rm = TRUE)
  }

  return(c(ci_low, ci_high))
}

