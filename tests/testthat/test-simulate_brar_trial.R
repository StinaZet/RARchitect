test_that("tuning 0 gives ER", {
 a = suppressWarnings(simulate_brar_trial(
     outcome_type = "binary",
     distribution = "bernoulli",
     arms = 2, N = 100, blocksize = 10,
     modelpar = c(0.6, 0.4),
     priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
     tuning = 0, clipping = 0, burnin = 0, randmethod = "coin",
     postprobmethod = "simulation",
     recruitment_rate = 5,
     observation_delay = 30)[,6])

 b = rep(0.5, 100)

 expect_equal(a, b)

})

test_that("simulate_brar_trial runs for all outcome types", {
  set.seed(1)

  expect_silent(simulate_brar_trial(
    outcome_type = "binary",
    distribution = "bernoulli",
    arms = 2, N = 40, blocksize = 1,
    modelpar = c(0.6, 0.4),
    priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
    tuning = 0, clipping = 0, burnin = 0, randmethod = "coin",
    postprobmethod = "exact", direction = "higher",
    recruitment_rate = 5,
    observation_delay = 30))

  expect_silent(simulate_brar_trial(
    outcome_type = "cont",
    distribution = "normal",
    arms = 2, N = 40, blocksize = 1,
    modelpar = matrix(c(0, 1, 1, 1), nrow = 2, byrow = TRUE),
    priors = matrix(c(0, 0, 0.1, 0.1, 1, 1, 1, 1), nrow = 4, byrow = TRUE),
    tuning = 1, clipping = 0, burnin = 0, randmethod = "coin",
    postprobmethod = "simulation",  direction = "higher",
    recruitment_rate = 5,
    observation_delay = 30))

  expect_silent(simulate_brar_trial(
    outcome_type = "cont",
    distribution = "exponential",
    arms = 2, N = 40, blocksize = 1,
    modelpar = c(1, 2),
    priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),tuning = 1,
    clipping = 0, burnin = 0, randmethod = "coin",
    postprobmethod = "simulation",  direction = "higher",
    recruitment_rate = 5,
    observation_delay = 30))
})

test_that("simulate_brar_trial is reproducible", {
  set.seed(123)
  x = simulate_brar_trial(
    outcome_type = "binary",
    distribution = "bernoulli",
    arms = 2, N = 40, blocksize = 1,
    modelpar = c(0.6, 0.4),
    priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
    tuning = 1, clipping = 0, burnin = 0, randmethod = "coin",
    postprobmethod = "exact", direction = "higher",
    recruitment_rate = 5,
    observation_delay = 30)

  set.seed(123)
  y = simulate_brar_trial(
    outcome_type = "binary",
    distribution = "bernoulli",
    arms = 2, N = 40, blocksize = 1,
    modelpar = c(0.6, 0.4),
    priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
    tuning = 1, clipping = 0, burnin = 0, randmethod = "coin",
    postprobmethod = "exact", direction = "higher",
    recruitment_rate = 5,
    observation_delay = 30)

  expect_equal(x$Arm, y$Arm)
  expect_equal(x$Outcome, y$Outcome)
})

library(testthat)

# --------------------------------------------------------------------------------
# Test simulate_brar_trial for basic functionality, edge cases, and output validity
# --------------------------------------------------------------------------------

test_that("simulate_brar_trial requires positive N", {
  expect_error(
    simulate_brar_trial(
      outcome_type = "binary",
      distribution = "bernoulli",
      direction = "higher",
      arms = 2,
      N = 0,                       # invalid input
      blocksize = 5,
      priors = matrix(c(1, 1, 1, 1), nrow = 2),
      modelpar = c(0.5, 0.5),
      tuning = 1,
      clipping = 0,
      burnin = 0,
      postprobmethod = "simulation",
      randmethod = "coin"
    ),
    "Parameter 'N' must be a positive integer."
  )
})


test_that("simulate_brar_trial produces correct number of rows and columns", {
  N = 20
  arms = 2
  blocksize = 1
  priors = matrix(c(1, 1, 1, 1), nrow=2)
  modelpar = c(0.5, 0.6)

  res <- suppressWarnings(simulate_brar_trial(
    outcome_type = "binary",
    distribution = "bernoulli",
    direction = "higher",
    arms = arms,
    N = N,
    blocksize = blocksize,
    priors = priors,
    modelpar = modelpar,
    tuning = 1,
    clipping = 0,
    burnin = 0,
    postprobmethod = "simulation",
    randmethod = "coin",
    urn_alpha = 1,
    multiarm_method = "fixed"
  ))

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), N)
  expect_equal(ncol(res), 5 + arms) # Batch, Arm, Outcome, time in, time out + allocation_probs
  expect_true(all(res$Arm %in% 1:arms))
})

test_that("simulate_brar_trial respects burn-in period", {
  N = 10
  burnin = 4
  arms = 2
  blocksize = 2
  priors = matrix(c(1, 1, 1, 1), nrow=2)
  modelpar = c(0.5, 0.6)

  res <- simulate_brar_trial(
    outcome_type = "binary",
    distribution = "bernoulli",
    direction = "higher",
    arms = arms,
    N = N,
    blocksize = blocksize,
    priors = priors,
    modelpar = modelpar,
    tuning = 1,
    clipping = 0,
    burnin = burnin,
    postprobmethod = "simulation",
    randmethod = "urn",
    urn_alpha = 1
  )

  # First burnin block should have roughly equal allocation probabilities
  expect_true(all(res$Batch[1:burnin] == 1))
  expect_true(all(res$EstimatedRecruitmentTime[1:burnin] >= 0)) # for Poisson
})

test_that("simulate_brar_trial supports multi-arm > 2 with top2 method", {
  N = 15
  arms = 3
  blocksize = 3
  priors = matrix(1, nrow=2, ncol=arms)
  modelpar = c(0.3, 0.5, 0.7)

  res <- simulate_brar_trial(
    outcome_type = "binary",
    distribution = "bernoulli",
    direction = "higher",
    arms = arms,
    N = N,
    blocksize = blocksize,
    priors = priors,
    modelpar = modelpar,
    tuning = 1,
    clipping = 0,
    burnin = 0,
    postprobmethod = "simulation",
    randmethod = "urn",
    urn_alpha = 1,
    multiarm_method = "top2"
  )

  expect_true(all(res$Arm %in% 1:arms))
  expect_equal(nrow(res), N)
  expect_equal(ncol(res), 5 + arms)
})

test_that("simulate_trial_duration_poisson_recruitment computes expected times", {
  N = 20
  block_sizes = rep(5, 4)
  poisson_rate = 5
  obs_delay = 10
  num_sim = 100

  res <- simulate_trial_duration_poisson_recruitment(
    N = N,
    block_sizes = block_sizes,
    poisson_recruitment_rate_per_unit_time = poisson_rate,
    outcome_observation_delay_per_patient = obs_delay,
    num_simulations = num_sim
  )

  expect_equal(nrow(res), N)
  expect_true(all(res$EstimatedRecruitmentTime >= 0))
  expect_true(all(res$EstimatedObservationTime >= res$EstimatedRecruitmentTime))
  expect_true(all(res$StdDevRecruitmentTime >= 0))
  expect_true(all(res$StdDevObservationTime >= 0))
})



# --------------------------------------------------------------------------------
# Stress Test: Large trial to check adaptive allocation
# --------------------------------------------------------------------------------
test_that("simulate_brar_trial adapts allocation to better performing arms in large trial", {
  set.seed(123)  # for reproducibility

  burnin = 20
  N = 500
  blocksize = 10
  arms = 2

  res <- simulate_brar_trial(
    outcome_type = "cont",
    distribution = "normal",
    direction = "higher",
    arms = arms,
    N = N,
    blocksize = blocksize,
    priors =  matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE),
    modelpar = matrix(c(0, 1, 1, 1), nrow = 2, byrow = TRUE),
    tuning = 1,
    clipping = 0,
    burnin = burnin,
    postprobmethod = "simulation",
    randmethod = "block",
    known_var = TRUE
  )

  # Compute proportion of allocations after burn-in
  allocation_after_burnin <- res$Arm[(burnin + 1):N]
  prop_arm2 <- mean(allocation_after_burnin == 2)
  prop_arm1 <- mean(allocation_after_burnin == 1)

  # Arm 2 is better, so should have a higher allocation proportion
  expect_gt(prop_arm2, prop_arm1)

  # Sanity check: all arm assignments valid
  expect_true(all(res$Arm %in% 1:arms))
})

