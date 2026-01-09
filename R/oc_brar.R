#' @title Operating Characteristics for a Bayesian BRAR Design
#'
#' @description
#' Simulates and analyzes \code{M} independent Bayesian response-adaptive randomization
#' (BRAR) trials and computes operating characteristics including treatment effect
#' performance, power/type I error, confidence interval coverage, and allocation behavior.
#'
#' @param M Integer. Number of Monte Carlo trial replicates.
#' @param true_effect Numeric vector. True treatment effect(s) for experimental arms
#'   (on the scale returned by the analyzer, e.g., risk difference).
#' @param seed Optional integer. Random seed for reproducibility.
#' @param plot Logical. If TRUE, produces a boxplot of sample sizes per arm.
#' @param ... Additional arguments passed to both
#'   \code{\link{simulate_brar_trial}} and \code{\link{analyze_brar_trial}}.
#'
#' @return A list with components:
#' \describe{
#'   \item{effect}{Mean estimates and Monte Carlo error.}
#'   \item{testing}{Power or type I error.}
#'   \item{coverage}{Confidence interval coverage.}
#'   \item{allocation}{Mean and distribution of sample sizes per arm.}
#'   \item{raw}{Raw Monte Carlo outputs.}
#' }
#'
#' @export
oc_brar <- function(M,
                    true_effect,
                    seed = NULL,
                    plot = TRUE,
                    ...) {

  if (M <= 0 || M %% 1 != 0) stop("'M' must be a positive integer.")
  if (!is.null(seed)) set.seed(seed)

  args = list(...)
  alpha = ifelse(is.null(args$alpha), 0.05, args$alpha)

  # --- storage ---
  eff_mat   = NULL
  pval_mat  = NULL
  cover_mat = NULL
  n_alloc   = NULL

  # --- main loop ---
  for (m in seq_len(M)) {

    trial_data = do.call(simulate_brar_trial, args)

    analysis = do.call(analyze_brar_trial,
                        c(list(trial_data = trial_data), args))

    res = analysis$test_results_df

    # --- treatment effects ---
    eff_mat = rbind(eff_mat, res$effect_estimate)

    # --- hypothesis tests ---
    pval_mat = rbind(pval_mat, res$p_value)

    # --- CI coverage ---
    cover_mat = rbind(
      cover_mat,
      (true_effect >= res$ci_low) & (true_effect <= res$ci_high)
    )

    # --- allocation ---
    alloc_tab = table(factor(trial_data$Arm, levels = seq_len(args$arms)))
    n_alloc = rbind(n_alloc, as.numeric(alloc_tab))
  }

  colnames(eff_mat)   = paste0("Arm", seq_len(ncol(eff_mat)) + 1)
  colnames(pval_mat)  = colnames(eff_mat)
  colnames(cover_mat) = colnames(eff_mat)
  colnames(n_alloc)   = paste0("Arm", seq_len(ncol(n_alloc)))

  # --- operating characteristics ---

  mean_eff = colMeans(eff_mat, na.rm = TRUE)
  mc_error = apply(eff_mat, 2, sd, na.rm = TRUE) / sqrt(M)

  rejection_prob = colMeans(pval_mat < alpha, na.rm = TRUE)
  coverage_prob  = colMeans(cover_mat, na.rm = TRUE)
  mean_alloc     = colMeans(n_alloc)

  oc = list(

    effect = data.frame(
      arm = names(mean_eff),
      mean_estimate = mean_eff,
      mc_error = mc_error,
      row.names = NULL
    ),

    testing = data.frame(
      arm = names(rejection_prob),
      rejection_prob = rejection_prob,
      row.names = NULL
    ),

    coverage = data.frame(
      arm = names(coverage_prob),
      coverage = coverage_prob,
      row.names = NULL
    ),

    allocation = list(
      mean_n = mean_alloc,
      all_n  = n_alloc
    )
  )

  # --- boxplot ---
  if (plot) {
    graphics::boxplot(
      n_alloc,
      names = colnames(n_alloc),
      ylab = "Number of patients",
      main = "Distribution of sample size per arm"
    )
  }

  return(list(
    oc = oc,
    raw = list(
      effect_estimates = eff_mat,
      p_values = pval_mat,
      coverage = cover_mat,
      allocations = n_alloc
    ),
    settings = args
  ))
}
