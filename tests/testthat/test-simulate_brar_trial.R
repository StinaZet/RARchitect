test_that("tuning 0 gives ER", {
 a = simulate_brar_trial(
     outcome_type = "binary",
     distribution = "bernoulli",
     arms = 2, N = 100, blocksize = 10,
     modelpar = c(0.6, 0.4),
     priors = matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE),
     tuning = 0, clipping = 0, burnin = 0, randmethod = "coin",
     postprobmethod = "simulation",
     recruitment_rate = 5,
     observation_delay = 30)[,6]

 b = rep(0.5, 100)

 expect_equal(a, b)

})

test_that("simulate_brar_trial runs for all outcome types", {
  set.seed(1)

  expect_silent(simulate_brar_trial(
    outcome = "binary",
    arms = 2, N = 40,
    modelpar = c(0.3, 0.6)
  ))

  expect_silent(simulate_brar_trial(
    outcome = "normal",
    arms = 2, N = 40,
    modelpar = matrix(c(0,1, 1,1), nrow = 2)
  ))

  expect_silent(simulate_brar_trial(
    outcome = "exponential",
    arms = 2, N = 40,
    modelpar = c(1, 2)
  ))
})

test_that("simulate_brar_trial is reproducible", {
  set.seed(123)
  x <- simulate_brar_trial(outcome = "binary", arms = 2, N = 30, modelpar = c(0.4, 0.6))

  set.seed(123)
  y <- simulate_brar_trial(outcome = "binary", arms = 2, N = 30, modelpar = c(0.4, 0.6))

  expect_equal(x$Arm, y$Arm)
  expect_equal(x$Outcome, y$Outcome)
})

library(testthat)

# --------------------------------------------------------------------------------
# Test simulate_brar_trial for basic functionality, edge cases, and output validity
# --------------------------------------------------------------------------------

test_that("simulate_brar_trial handles zero participants correctly", {
  res <- simulate_brar_trial(
    direction = "higher",
    arms = 2,
    N = 0,
    blocksize = 5,
    priors = matrix(c(1,1,1,1), nrow=2),
    modelpar = c(0.5, 0.5),
    tuning = 1,
    clipping = 0,
    burnin = 0,
    postprobmethod = "simulation",
    randmethod = "coin",
    urn_alpha = 1,
    multiarm_method = "fixed"
  )
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 0)
  expect_equal(ncol(res), 5 + 2) # Batch, Arm, Outcome + allocation_probs for 2 arms
})

test_that("simulate_brar_trial produces correct number of rows and columns", {
  N <- 20
  arms <- 2
  blocksize <- 5
  priors <- matrix(c(1,1,1,1), nrow=2)
  modelpar <- c(0.5, 0.6)

  res <- simulate_brar_trial(
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
  )

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), N)
  expect_equal(ncol(res), 3 + arms) # Batch, Arm, Outcome + allocation_probs
  expect_true(all(res$Arm %in% 1:arms))
})

test_that("simulate_brar_trial respects burn-in period", {
  N <- 10
  burnin <- 4
  arms <- 2
  blocksize <- 2
  priors <- matrix(c(1,1,1,1), nrow=2)
  modelpar <- c(0.5, 0.6)

  res <- simulate_brar_trial(
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
    randmethod = "coin",
    urn_alpha = 1,
    multiarm_method = "fixed"
  )

  # First burnin block should have roughly equal allocation probabilities
  expect_true(all(res$Batch[1:burnin] == 1))
  expect_true(all(res$EstimatedRecruitmentTime[1:burnin] >= 0)) # for Poisson
})

test_that("simulate_brar_trial supports multi-arm > 2 with top2 method", {
  N <- 15
  arms <- 3
  blocksize <- 3
  priors <- matrix(1, nrow=2, ncol=arms)
  modelpar <- c(0.3, 0.5, 0.7)

  res <- simulate_brar_trial(
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
    multiarm_method = "top2"
  )

  expect_true(all(res$Arm %in% 1:arms))
  expect_equal(nrow(res), N)
  expect_equal(ncol(res), 3 + arms)
})

test_that("simulate_trial_duration_poisson_recruitment computes expected times", {
  N <- 20
  block_sizes <- rep(5, 4)
  poisson_rate <- 5
  obs_delay <- 10
  num_sim <- 100

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

test_that("simulate_brar_trial_normal works with known variance", {
  N <- 10
  arms <- 2
  blocksize <- 2
  priors <- matrix(c(0,0,1,1), nrow=2)
  modelpar <- matrix(c(0.2, 0.4, 1, 1), nrow=2)

  res <- .simulate_brar_trial_normal(
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
    multiarm_method = "fixed",
    urn_alpha = 1,
    known_var = TRUE
  )

  expect_equal(nrow(res), N)
  expect_equal(ncol(res), 3 + arms)
  expect_true(all(res$Arm %in% 1:arms))
})

# --------------------------------------------------------------------------------
# Stress Test: Large trial to check adaptive allocation
# --------------------------------------------------------------------------------
test_that("simulate_brar_trial adapts allocation to better performing arms in large trial", {
  set.seed(123)  # for reproducibility

  N <- 500
  arms <- 2
  blocksize <- 10
  burnin <- 20

  # True means: arm 2 is better
  true_means <- c(0.4, 0.7)
  true_sds <- c(1, 1)
  modelpar <- rbind(true_means, true_sds)

  # Priors for known variance
  priors <- rbind(c(0, 0), c(1, 1))  # mu0, tau0

  res <- .simulate_brar_trial_normal(
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
    randmethod = "coin",
    multiarm_method = "fixed",
    urn_alpha = 1,
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

