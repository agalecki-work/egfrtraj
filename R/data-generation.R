# data-generation.R
# Functions for generating synthetic longitudinal eGFR trajectories
# for testing and demonstration purposes

#' Generate a single synthetic eGFR trajectory
#'
#' Internal helper function to create one realistic-looking eGFR time series
#' with specified pattern.
#'
#' @param id Identifier for the trajectory (character or numeric)
#' @param pattern One of: "linear", "quadratic", "with_breakpoint", "insufficient_data"
#' @param years Numeric vector of time points (default: 0:26)
#' @param baseline_egfr Starting eGFR value (default: 125)
#' @param noise_sd Standard deviation of Gaussian noise (default: 2)
#'
#' @return A tibble with columns: id, time, egfr, true_pattern
#' @keywords internal
generate_one_trajectory <- function(
    id,
    pattern,
    years = 0:26,
    baseline_egfr = 125,
    noise_sd = 2
) {
  e <- rnorm(length(years), sd = noise_sd)
  
  if (pattern == "linear") {
    egfr <- baseline_egfr - 1.1 * years + e
  } else if (pattern == "quadratic") {
    egfr <- baseline_egfr - 0.8 * years - 0.08 * years^2 + e
  } else if (pattern == "with_breakpoint") {
    bp <- sample(10:18, 1)
    before_idx <- years < bp
    after_idx  <- years >= bp
    
    before <- baseline_egfr - 0.8 * years[before_idx]
    after  <- before[length(before)] - 8.5 * (years[after_idx] - bp)
    
    egfr <- numeric(length(years))
    egfr[before_idx] <- before
    egfr[after_idx]  <- after
    egfr <- egfr + e
  } else if (pattern == "insufficient_data") {
    egfr <- rep(NA_real_, length(years))
    n_keep <- sample(2:4, 1)
    egfr[seq_len(n_keep)] <- baseline_egfr - 1.5 * years[seq_len(n_keep)] +
                             rnorm(n_keep, sd = noise_sd * 1.3)
  } else {
    stop("Unknown pattern. Choose from: linear, quadratic, with_breakpoint, insufficient_data")
  }
  
  tibble::tibble(
    id          = id,
    time        = years,
    egfr        = egfr,
    true_pattern = pattern
  )
}


#' Generate synthetic eGFR dataset for testing
#'
#' Creates a balanced dataset with one trajectory per pattern (total 108 rows).
#'
#' @param n_per_pattern Number of trajectories per pattern type (default: 1)
#' @param years Time points (default: 0:26)
#' @param baseline_egfr Starting eGFR value (default: 125)
#' @param noise_sd Standard deviation of measurement noise (default: 2)
#' @param seed Random seed for reproducibility (default: NULL = no fixed seed)
#'
#' @return A tibble with columns: id, time, egfr, true_pattern
#' @export
#' @examples
#' set.seed(20250116)
#' sim_data <- generate_synthetic_egfr()
#' table(sim_data$true_pattern)
generate_synthetic_egfr <- function(
    n_per_pattern = 1,
    years = 0:26,
    baseline_egfr = 125,
    noise_sd = 0.1,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  patterns <- c("linear", "quadratic", "with_breakpoint", "insufficient_data")
  
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
#' A small pre-generated dataset with one trajectory per pattern (total 108 rows).
#'
#' @format A tibble with 108 rows and 4 columns:
#' \describe{
#'   \item{id}{Trajectory identifier}
#'   \item{time}{Time in years (0–26)}
#'   \item{egfr}{Estimated GFR value}
#'   \item{true_pattern}{True underlying pattern ("linear", "quadratic", "with_breakpoint", "insufficient_data")}
#' }
#' @source Generated with default parameters (noise_sd = 2, n_per_pattern = 1)
#' @keywords datasets
"example_egfr_trajectories"
