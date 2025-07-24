test_that("tuning 0 gives ER", {
 a = simulate_brar_trial(
     outcome_type = "binary",
     arms = 2, N = 100, blocksize = 10,
     modelpar = c(0.6, 0.4),
     priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
     tuning = 0, clipping = 0, burnin = 0, ensure_all_arms_sampled = FALSE,
     postprobmethod = "simulation",
     recruitment_rate = 5,
     observation_delay = 30)[,6]

 b = rep(0.5, 100)

 expect_equal(a, b)

})
