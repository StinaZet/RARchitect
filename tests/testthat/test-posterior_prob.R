### Unit tests for posterior probability updates for binary outcome. ###

test_that("posterior_bin_sim returns valid probabilities", {
  set.seed(1)
  ap <- posterior_bin_sim(M = 5000, alphas = c(5,5), betas = c(5,5), direction = "higher")

  expect_length(ap, 2)
  expect_equal(sum(ap), 1, tolerance = 0.01)
  expect_true(all(ap >= 0))
})

test_that("posterior_bin_exact symmetric case gives 0.5", {
  ap <- posterior_bin_exact(c(5,5), c(5,5), direction = "higher")
  expect_equal(ap, c(0.5, 0.5))
})

test_that("posterior_bin_sim favors better arm", {
  set.seed(2)
  ap <- posterior_bin_sim(M = 5000, alphas = c(20,5), betas = c(5,20), direction = "higher")
  expect_gt(ap[1], ap[2])
})

test_that("direction = 'lower' flips probabilities", {
  ap1 <- posterior_bin_exact(c(10,2), c(2,10), direction = "higher")
  ap2 <- posterior_bin_exact(c(10,2), c(2,10), direction = "lower")
  expect_equal(ap1, rev(ap2))
})

### Stop binary ###

### Unit tests for posterior probability updates for normal outcome. ###

test_that("posterior_norm_exact symmetric case", {
  ap <- posterior_norm_exact(c(0,0), c(1,1), direction = "higher")
  expect_equal(ap, c(0.5, 0.5))
})

test_that("posterior_norm_sim returns valid probabilities", {
  set.seed(1)
  ap <- posterior_norm_sim(5000, means = c(0,1), sds = c(1,1), direction = "higher")
  expect_equal(sum(ap), 1, tolerance = 0.02)
  expect_gt(ap[2], ap[1])
})

test_that("posterior_norm_unknownvar_sim works for 3 arms", {
  set.seed(1)
  ap <- posterior_norm_unknownvar_sim(
    M = 3000,
    mu_n = c(0, 0.5, 1),
    kappa_n = c(5,5,5),
    alpha_n = c(5,5,5),
    beta_n = c(1,1,1), direction = "higher"
  )

  expect_length(ap, 3)
  expect_equal(sum(ap), 1, tolerance = 0.02)
})

### Stop normal ###

### Unit tests for posterior probability updates for exponential outcome. ###

test_that("posterior_exp_exact symmetric case", {
  ap <- posterior_exp_exact(c(5,5), c(2,2), direction = "higher")
  expect_equal(ap, c(0.5, 0.5))
})

test_that("posterior_exp_sim favors larger mean", {
  set.seed(1)
  ap <- posterior_exp_sim(5000, shapes = c(10,2), rates = c(1,1), direction = "higher")
  expect_gt(ap[1], ap[2])
})

### Stop exponential ###

### Unit tests for Top 2 TS. ###

test_that("alloc_probs_T2TS returns valid probabilities", {
  raw <- c(0.7, 0.2, 0.1)
  ap <- alloc_probs_T2TS(raw, beta = 0.5)

  expect_length(ap, 3)
  expect_equal(sum(ap), 1)
  expect_true(all(ap >= 0))
})

test_that("alloc_probs_T2TS preserves symmetry", {
  raw <- c(1/3, 1/3, 1/3)
  ap <- alloc_probs_T2TS(raw, beta = 0.5)

  expect_equal(ap, raw)
})

### Stop T2TS ###
