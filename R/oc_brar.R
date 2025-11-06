oc_brar <- function(
    # --- 1. Monte Carlo Simulation Arguments ---
  M,                      # Numeric: Number of replicates to run.
  # --- 2. Analysis & Testing Arguments ---
  estimation_method = c("posterior_mean", "MLE", "bayesian_median"), # Character: Method to estimate the effect for each arm.
  test_method = c("bayesian_pp", "frequentist_z", "frequentist_t"), # Character: Hypothesis testing method.
  alpha = 0.05,           # Numeric: Significance level (frequentist Type I Error rate) or confidence level (Bayesian threshold).
  prob_threshold = 0.95,  # Numeric: The posterior probability threshold for Bayesian testing (e.g., P(best | Data) > 0.95).
  null_effect = 0,        # Numeric: Value representing the null hypothesis (e.g., H0: mu_A - mu_C = 0).

  # --- 3. All Arguments from simulate_XXX_trial ---
  # (These will be passed directly to the single-replicate function)
  outcome_type = c("binary", "cont"),
  distribution = c("bernoulli", "normal", "exponential"),
  arms, N, blocksize, direction = c("lower", "higher"),
  known_var = FALSE, priors, modelpar, tuning = 1,
  clipping = 0, burnin = 0, randmethod = "coin", urn_alpha = 3,
  postprobmethod = "simulation", multiarm_method = c("fixed", "top2"),
  recruitment_rate = 100000, observation_delay = 0,
  ...
) {
  # ... Function implementation ...
}
