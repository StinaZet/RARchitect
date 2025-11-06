# Randomization-based hypothesis test.
#' @title Perform a Randomization-Based Permutation Test

perform_permutation_test <- function(y_c, n_c, y_e, n_e,
                                      alternative = "two.sided", B = 1000) {

  # 1. Calculate the observed test statistic (Difference in Proportions)
  # Handle division by zero if n_c or n_e is 0 (though analyze_brar_trial checks this)
  if (n_e == 0 || n_c == 0) return(NA_real_)
  p_obs <- y_e/n_e - y_c/n_c

  # 2. Pool the data (Outcomes under H0)
  n_total <- n_c + n_e
  y_total <- y_c + y_e

  outcomes <- c(rep(1, y_total), rep(0, n_total - y_total))

  perm_stats <- numeric(B)

  # 3. Permutation Loop
  for (i in 1:B) {
    perm_outcomes <- sample(outcomes)
    y_e_perm <- sum(perm_outcomes[1:n_e])
    y_c_perm <- sum(perm_outcomes[(n_e + 1):n_total])

    perm_stats[i] <- y_e_perm/n_e - y_c_perm/n_c
  }

  # 4. Calculate p-value
  if (alternative == "greater") {
    # One-sided test for p_e > p_c (positive difference)
    p_value <- (sum(perm_stats >= p_obs) + 1) / (B + 1) # Add +1 for obs
  } else if (alternative == "less") {
    # One-sided test for p_e < p_c (negative difference)
    p_value <- (sum(perm_stats <= p_obs) + 1) / (B + 1) # Add +1 for obs
  } else {
    # Two-sided p-value: P(|T_perm| >= |T_obs|)
    p_value <- (sum(abs(perm_stats) >= abs(p_obs)) + 1) / (B + 1) # Add +1 for obs
  }

  return(p_value)
}
