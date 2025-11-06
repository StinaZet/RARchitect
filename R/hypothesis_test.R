# Randomization-based hypothesis test.
perform_permutation_test <- function(y_c, n_c, y_e, n_e, direction, B = 1000) {

  # 1. Calculate the observed test statistic (Difference in Proportions)
  p_obs <- y_e/n_e - y_c/n_c

  # 2. Pool the data (Outcomes under H0)
  n_total <- n_c + n_e
  y_total <- y_c + y_e

  # Create the vector of all observed outcomes (1 for success, 0 for failure)
  outcomes <- c(rep(1, y_total), rep(0, n_total - y_total))

  perm_stats <- numeric(B)

  # 3. Permutation Loop
  for (i in 1:B) {
    # Randomly shuffle the pooled outcomes
    perm_outcomes <- sample(outcomes)

    # Assign the first n_e outcomes to the Experimental arm and the rest to Control
    y_e_perm <- sum(perm_outcomes[1:n_e])
    y_c_perm <- sum(perm_outcomes[(n_e + 1):n_total])

    # Calculate the permuted difference in proportions
    perm_stats[i] <- y_e_perm/n_e - y_c_perm/n_c
  }

  # 4. Calculate p-value
  # P-value is the proportion of permuted statistics as extreme or more extreme than the observed.
  if (direction == "higher") {
    # One-sided test for p_e > p_c (positive difference)
    p_value <- sum(perm_stats >= p_obs) / B
  } else {
    # One-sided test for p_e < p_c (negative difference)
    p_value <- sum(perm_stats <= p_obs) / B
  }

  return(p_value)
}
