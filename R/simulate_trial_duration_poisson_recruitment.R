#' @title Simulate Trial Durations with Poisson Recruitment and Variable Block Sizes
#' @description
#' Simulates the total duration of a clinical trial multiple times, accounting for
#' patient recruitment following a Poisson process (meaning inter-arrival times
#' are exponentially distributed) and a fixed time for outcome observation.
#' Recruitment occurs in sequential blocks of potentially different sizes,
#' and a new block can only begin enrollment after all participants from the
#' previous block have been fully enrolled and their outcomes observed.
#' This function returns the estimated (average) recruitment and observation times
#' for each participant across all simulations, along with their standard deviations.
#'
#' @param N Numeric. Total sample size for the trial. Must be equal to `sum(block_sizes)`.
#' @param block_sizes Numeric vector. A vector where each element specifies the
#' size of a sequential recruitment block. The sum of all elements in `block_sizes`
#' must equal `N`. All elements must be positive integers.
#' @param poisson_recruitment_rate_per_unit_time Numeric. The rate parameter (lambda)
#' for the Poisson process that governs patient recruitment. This is the average
#' number of patients arriving per unit of time (e.g., patients per month).
#' Inter-arrival times for patient recruitment are drawn from an exponential
#' distribution with this rate.
#' @param outcome_observation_delay_per_patient Numeric. The fixed time from a participant's
#' enrollment until their outcome is observed. For example, if 30, it takes
#' 30 days for an outcome to be observed after a participant is enrolled.
#' @param num_simulations Integer. The number of times to run the simulation.
#' A higher number provides a more stable estimate of the average recruitment
#' and observation times for each participant and their variability.
#'
#' @return A data frame with `N` rows and the following columns:
#' \itemize{
#' \item `SubjectID`: The unique identifier for each participant (1 to N).
#' \item `EstimatedRecruitmentTime`: The average time point at which each participant
#' is recruited, estimated across all simulation runs.
#' \item `StdDevRecruitmentTime`: The standard deviation of recruitment times
#' for each participant across all simulation runs. If `num_simulations` is 1,
#' this will be 0.
#' \item `EstimatedObservationTime`: The average time point at which each participant's
#' outcome is observed, estimated across all simulation runs.
#' \item `StdDevObservationTime`: The standard deviation of observation times
#' for each participant across all simulation runs. If `num_simulations` is 1,
#' this will be 0.
#' }
#' If `N` is 0, an empty data frame with the specified columns is returned.
#' The estimated total trial duration can be obtained by taking `max(EstimatedObservationTime)`
#' from the returned data frame.
#' @export
#'
#' @examples
#' # Example 1: Fixed block size (same as previous behavior, but using vector)
#' set.seed(123) # For reproducibility
#' estimated_times_1 <- simulate_trial_duration_poisson_recruitment(
#' N = 200,
#' block_sizes = rep(10, 20), # 20 blocks of 10 patients each
#' poisson_recruitment_rate_per_unit_time = 5,
#' outcome_observation_delay_per_patient = 30,
#' num_simulations = 1000
#' )
#' estimated_total_duration_1 <- max(estimated_times_1$EstimatedObservationTime)
#' print(estimated_total_duration_1)
#' head(estimated_times_1)
#' summary(estimated_times_1$StdDevRecruitmentTime) # Check variability
#'
#' # Example 2: Variable block sizes
#' set.seed(456)
#' estimated_times_2 <- simulate_trial_duration_poisson_recruitment(
#' N = 100,
#' block_sizes = c(10, 15, 20, 25, 30), # Sum is 100
#' poisson_recruitment_rate_per_unit_time = 10,
#' outcome_observation_delay_per_patient = 15,
#' num_simulations = 500
#' )
#' estimated_total_duration_2 <- max(estimated_times_2$EstimatedObservationTime)
#' print(estimated_total_duration_2)
#' tail(estimated_times_2)
#' summary(estimated_times_2$StdDevObservationTime) # Check variability
#'
#' # Example 3: Single patient trial, single simulation (StdDev will be 0)
#' set.seed(789)
#' estimated_times_3 <- simulate_trial_duration_poisson_recruitment(
#' N = 1,
#' block_sizes = c(1), # Single block of 1 patient
#' poisson_recruitment_rate_per_unit_time = 2,
#' outcome_observation_delay_per_patient = 10,
#' num_simulations = 1
#' )
#' estimated_total_duration_3 <- max(estimated_times_3$EstimatedObservationTime)
#' print(estimated_total_duration_3)
#' print(estimated_times_3)
simulate_trial_duration_poisson_recruitment <- function(N,
                                                        block_sizes,
                                                        poisson_recruitment_rate_per_unit_time,
                                                        outcome_observation_delay_per_patient,
                                                        num_simulations) {

  # Input Validation
  if (N < 0) {
    stop("Parameter 'N' cannot be negative.")
  }
  if (!is.numeric(block_sizes) || any(block_sizes <= 0) || any(block_sizes %% 1 != 0)) {
    stop("Parameter 'block_sizes' must be a numeric vector of positive integers.")
  }
  if (sum(block_sizes) != N) {
    stop("The sum of 'block_sizes' must be equal to 'N' (total sample size).")
  }
  if (poisson_recruitment_rate_per_unit_time <= 0) {
    stop("Parameter 'poisson_recruitment_rate_per_unit_time' must be positive.")
  }
  if (outcome_observation_delay_per_patient < 0) {
    stop("Parameter 'outcome_observation_delay_per_patient' cannot be negative.")
  }
  if (num_simulations <= 0 || num_simulations %% 1 != 0) {
    stop("Parameter 'num_simulations' must be a positive integer.")
  }

  # Handle edge case: no patients
  if (N == 0) {
    return(data.frame(
      SubjectID = numeric(0),
      EstimatedRecruitmentTime = numeric(0),
      StdDevRecruitmentTime = numeric(0),
      EstimatedObservationTime = numeric(0),
      StdDevObservationTime = numeric(0)
    ))
  }

  # Initialize matrices to store times from each simulation run for calculating variability
  all_recruitment_times <- matrix(NA, nrow = N, ncol = num_simulations)
  all_observation_times <- matrix(NA, nrow = N, ncol = num_simulations)

  num_batches <- length(block_sizes)

  for (s in 1:num_simulations) {
    # Initialize for each simulation run
    current_global_time <- 0
    simulated_recruitment_times_for_run <- numeric(N)
    simulated_observation_times_for_run <- numeric(N)

    # Track the number of patients already assigned to previous blocks
    patients_processed_so_far <- 0

    # Loop through each batch within the current simulation run
    for (b_idx in 1:num_batches) {
      current_batch_size <- block_sizes[b_idx]

      start_index <- patients_processed_so_far + 1
      end_index <- patients_processed_so_far + current_batch_size

      # Simulate inter-arrival times for patients within this batch
      inter_arrival_times_in_batch <- stats::rexp(current_batch_size,
                                           rate = poisson_recruitment_rate_per_unit_time)

      # Calculate cumulative recruitment times within this batch relative to its start
      cumulative_arrival_times_relative_to_batch_start <- cumsum(inter_arrival_times_in_batch)

      # Adjust to global time for the current simulation run
      batch_recruitment_times_global <- current_global_time +
        cumulative_arrival_times_relative_to_batch_start

      # Store recruitment times for the current run
      simulated_recruitment_times_for_run[start_index:end_index] <- batch_recruitment_times_global

      # Calculate observation times for subjects in this batch
      batch_observation_times_global <- batch_recruitment_times_global +
        outcome_observation_delay_per_patient

      # Store observation times for the current run
      simulated_observation_times_for_run[start_index:end_index] <- batch_observation_times_global

      # Update 'current_global_time' for the next batch.
      # This enforces the blocked accrual: next batch starts after *all* outcomes
      # from the current batch are observed.
      current_global_time <- max(batch_observation_times_global)

      # Update count of processed patients for the next batch
      patients_processed_so_far <- patients_processed_so_far + current_batch_size
    }

    # Store the results of the current simulation run in the matrices
    all_recruitment_times[, s] <- simulated_recruitment_times_for_run
    all_observation_times[, s] <- simulated_observation_times_for_run
  }

  # Calculate the average recruitment and observation times across all simulations
  estimated_recruitment_times <- rowMeans(all_recruitment_times)
  estimated_observation_times <- rowMeans(all_observation_times)

  # Calculate standard deviations. Handle num_simulations = 1 gracefully (sd will be 0).
  if (num_simulations > 1) {
    std_dev_recruitment_times <- apply(all_recruitment_times, 1, stats::sd)
    std_dev_observation_times <- apply(all_observation_times, 1, stats::sd)
  } else {
    std_dev_recruitment_times <- rep(0, N)
    std_dev_observation_times <- rep(0, N)
  }

  # Round recruitment time and outcome time to nearest integer above, except if the recruitment is very fast, then we round down (emulating instant recruitment).
  estimated_recruitment_times = ifelse(max(estimated_recruitment_times < 0.1), 0, ceiling(estimated_recruitment_times))
  estimated_observation_times = ifelse(max(estimated_recruitment_times < 0.1), floor(estimated_observation_times), ceiling(estimated_observation_times))

  # Create the final data frame
  estimated_subject_details_df <- data.frame(
    SubjectID = 1:N,
    EstimatedRecruitmentTime = estimated_recruitment_times,
    StdDevRecruitmentTime = std_dev_recruitment_times,
    EstimatedObservationTime = estimated_observation_times,
    StdDevObservationTime = std_dev_observation_times
  )
  colnames(estimated_subject_details_df) = c("SubjectID", "EstimatedRecruitmentTime", "StdDevRecruitmentTime", "EstimatedObservationTime", "StdDevObservationTime")
  return(estimated_subject_details_df)
  #return( data.frame(
  #  SubjectID = 1:N,
  #  EstimatedRecruitmentTime = estimated_recruitment_times,
  #  StdDevRecruitmentTime = std_dev_recruitment_times,
  #  EstimatedObservationTime = estimated_observation_times,
  #  StdDevObservationTime = std_dev_observation_times
  #))
}
