#' @title Operating Characteristics for a Bayesian BRAR Design
#'
#' @description
#' Simulates and analyzes \code{M} independent Bayesian response-adaptive randomization
#' (BRAR) trials using \code{\link{simulate_brar_trial}} and
#' \code{\link{analyze_brar_trial}}, and summarizes operating characteristics such as
#' rejection probabilities, allocation proportions, patient benefit, and trial duration.
#'
#' @param M Integer. Number of Monte Carlo trial replicates.
#' @param seed Optional integer. Random seed for reproducibility.
#' @param ... Additional arguments passed to both
#'   \code{\link{simulate_brar_trial}} and \code{\link{analyze_brar_trial}}.
#'
#' @return A list with components:
#' \describe{
#'   \item{raw_results}{List containing per-replicate trial data and analysis results.}
#'   \item{oc}{Named list of summarized operating characteristics.}
#'   \item{settings}{List of design and analysis parameters used.}
#' }
#'
#' @details
#' For each replicate, the function:
#' \enumerate{
#'   \item Simulates a BRAR trial using \code{simulate_brar_trial()}.
#'   \item Analyzes the trial using \code{analyze_brar_trial()}.
#'   \item Extracts quantities of interest (e.g., p-values, allocation, patient benefit).
#' }
#'
#' By default, the following operating characteristics are reported:
#' \itemize{
#'   \item Rejection probability (per tested arm)
#'   \item Mean patient benefit
#'   \item Mean overall outcome
#'   \item Mean and SD of allocation proportions
#'   \item Mean trial duration (if outcome times are present)
#' }
#'
#' This function is designed to be extensible. Users can easily add new operating
#' characteristics such as power, type I error, probability of correct selection,
#' expected sample size, or bias and RMSE of estimators.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' oc <- oc_brar(
#'   M = 100,
#'   seed = 1,
#'   outcome_type = "binary",
#'   distribution = "bernoulli",
#'   direction = "higher",
#'   arms = 2, N = 100, blocksize = 10,
#'   modelpar = c(0.5, 0.65),
#'   priors = matrix(c(1, 1, 1, 1), nrow = 2),
#'   tuning = 1, clipping = 0, burnin = 0,
#'   randmethod = "coin",
#'   postprobmethod = "exact",
#'   estimation_method = "IPW",
#'   test_method = "simulation",
#'   CI_method = "simulation",
#'   B = 500, alpha = 0.05
#' )
#'
#' oc$oc
#' }
oc_brar <- function(M,
                    seed = NULL,
                    ...) {

  if (missing(M) || M <= 0 || M %% 1 != 0)
    stop("'M' must be a positive integer.")

  if (!is.null(seed)) set.seed(seed)

  args = list(...)

  # --- storage ---
  trial_data_list  = vector("list", M)
  analysis_list    = vector("list", M)

  reject_mat      = NULL
  patient_benefit = numeric(M)
  mean_outcome    = numeric(M)
  alloc_props     = NULL
  trial_duration  = rep(NA_real_, M)

  # --- main loop ---
  for (m in seq_len(M)) {

    trial_data = do.call(simulate_brar_trial, args)

    analysis = do.call(analyze_brar_trial,
                        c(list(trial_data = trial_data), args))

    trial_data_list[[m]] = trial_data
    analysis_list[[m]]   = analysis

    # --- rejection indicators ---
    if (!is.null(analysis$test_results_df$p_value)) {
      reject_mat = rbind(reject_mat,
                          analysis$test_results_df$p_value < args$alpha)
    }

    # --- patient benefit and outcome ---
    patient_benefit[m] = analysis$patient_benefit
    mean_outcome[m]    = analysis$mean_outcome

    # --- allocation proportions ---
    alloc_tab = table(trial_data$Arm)
    alloc_props = rbind(alloc_props, alloc_tab / sum(alloc_tab))

    # --- trial duration ---
    if ("Outcome time" %in% names(trial_data)) {
      trial_duration[m] = max(trial_data$`Outcome time`)
    }
  }

  # --- operating characteristics ---
  oc = list(
    rejection_prob = colMeans(reject_mat, na.rm = TRUE),
    mean_patient_benefit = mean(patient_benefit, na.rm = TRUE),
    sd_patient_benefit   = stats::sd(patient_benefit, na.rm = TRUE),
    mean_outcome = mean(mean_outcome, na.rm = TRUE),
    mean_allocation = colMeans(alloc_props, na.rm = TRUE),
    sd_allocation   = apply(alloc_props, 2, stats::sd, na.rm = TRUE),
    mean_trial_duration = if (all(is.na(trial_duration)))
      NULL else mean(trial_duration, na.rm = TRUE)
  )

  return(list(
    raw_results = list(
      trial_data = trial_data_list,
      analysis   = analysis_list
    ),
    oc = oc,
    settings = args
  ))
}
