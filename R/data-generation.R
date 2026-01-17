# data-generation.R
# Functions for generating synthetic longitudinal eGFR trajectories
# for testing and demonstration purposes

#' Generate a single synthetic eGFR trajectory
#'
#' Internal helper function to create one realistic-looking eGFR time series
#' with specified pattern.
#'
#' @param id Identifier for the trajectory (character or numeric)
#' @param pattern One of: "stable", "slow_linear", "rapid_breakpoint", "insufficient"
#' @param years Numeric vector of time points (default: 0:26)
#' @param baseline_egfr Starting eGFR value (default: 125)
#' @param noise_sd Standard deviation of Gaussian noise (default: 3)
#'
#' @return A tibble with columns: id, year, eGFR, true_pattern
#' @keywords internal
generate_one_trajectory <- function(
    id,
    pattern,
    years = 0:26,
    baseline_egfr = 125,
    noise_sd = 3
) {
  e <- rnorm(length(years), sd = noise_sd)
  
  if (pattern == "stable") {
    egfr <- baseline_egfr + 0.2 * years + e
  } else if (pattern == "slow_linear") {
    egfr <- baseline_egfr - 1.2 * years + e
  } else if (pattern == "rapid_breakpoint") {
    # Random breakpoint between year 10 and 18 for variety
    bp <- sample(10:18, 1)
    before_idx <- years < bp
    after_idx  <- years >= bp
    
    before <- baseline_egfr - 0.8 * years[before_idx]
    after  <- before[length(before)] - 5.5 * (years[after_idx] - bp)
    
    egfr <- numeric(length(years))
    egfr[before_idx] <- before
    egfr[after_idx]  <- after
    egfr <- egfr + e
  } else if (pattern == "insufficient") {
    # Simulate missing data / short follow-up
    egfr <- rep(NA_real_, length(years))
    n_keep <- sample(3:6, 1)
    egfr[seq_len(n_keep)] <- baseline_egfr - 1.5 * years[seq_len(n_keep)] +
                             rnorm(n_keep, sd = 4)
  } else {
    stop("Unknown pattern. Choose from: stable, slow_linear, rapid_breakpoint, insufficient")
  }
  
  tibble::tibble(
    id          = id,
    year        = years,
    eGFR        = egfr,
    true_pattern = pattern
  )
}


#' Generate synthetic eGFR dataset for testing
#'
#' Creates a balanced dataset containing multiple trajectories of each type.
#' Useful for testing classification algorithms and creating vignettes/examples.
#'
#' @param n_per_pattern Number of trajectories per pattern type (default: 4)
#' @param years Time points (default: 0:26 → 27 annual measurements)
#' @param baseline_egfr Starting eGFR value (default: 125)
#' @param noise_sd Standard deviation of measurement noise (default: 3)
#' @param seed Random seed for reproducibility (default: NULL = no fixed seed)
#'
#' @return A tibble with columns: id, year, eGFR, true_pattern
#' @export
#' @examples
#' set.seed(20250116)
#' sim_data <- generate_synthetic_egfr(n_per_pattern = 5)
#' table(sim_data$true_pattern)
generate_synthetic_egfr <- function(
    n_per_pattern = 4,
    years = 0:26,
    baseline_egfr = 125,
    noise_sd = 3,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  patterns <- c("stable", "slow_linear", "rapid_breakpoint", "insufficient")
  
  ids <- seq_len(n_per_pattern * length(patterns))
  pattern_vec <- rep(patterns, each = n_per_pattern)
  
  purrr::map2_dfr(
    .x = ids,
    .y = pattern_vec,
    .f = ~ generate_one_trajectory(
      id = .x,
      pattern = .y,
      years = years,
      baseline_egfr = baseline_egfr,
      noise_sd = noise_sd
    )
  ) %>%
    dplyr::mutate(id = factor(id, levels = unique(id)))
}


#' Small example synthetic dataset (for package data)
#'
#' A small pre-generated dataset for immediate use in examples and tests.
#'
#' @format A tibble with 108 rows and 4 columns:
#' \describe{
#'   \item{id}{Trajectory identifier}
#'   \item{year}{Time in years (0–26)}
#'   \item{eGFR}{Estimated GFR value}
#'   \item{true_pattern}{True underlying pattern}
#' }
#' @source Generated with default parameters
#' @keywords datasets
"example_egfr_trajectories"
