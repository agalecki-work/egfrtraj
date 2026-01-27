# data-generation.R
# Functions for generating synthetic longitudinal eGFR trajectories
# for testing and demonstration purposes

#' Generate a single synthetic eGFR trajectory
#'
#' Internal helper function to create one realistic-looking eGFR time series
#' with specified pattern.
#'
#' @param id Identifier for the trajectory (character or numeric)
#' @param pattern One of: "linear", "quadratic", "with_breakpoint", MYNOTE: Do not create data with "insufficient_data" pattern
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
    after  <- before[length(before)] - 6.5 * (years[after_idx] - bp)
    
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


# data-generation.R
# Functions for generating synthetic longitudinal eGFR trajectories

#' Generate a single synthetic eGFR trajectory
#'
#' Internal helper function to create one realistic-looking eGFR time series
#' with specified pattern.
#'
#' @param id Identifier for the trajectory (character or numeric)
#' @param pattern One of: "linear", "quadratic", "with_breakpoint"
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
    after  <- before[length(before)] - 6.5 * (years[after_idx] - bp)
    
    egfr <- numeric(length(years))
    egfr[before_idx] <- before
    egfr[after_idx]  <- after
    egfr <- egfr + e
  } else {
    stop("Unknown pattern. Choose from: linear, quadratic, with_breakpoint")
  }
  
  tibble::tibble(
    id          = as.character(id),   # ensure character
    time        = years,
    egfr        = egfr,
    true_pattern = pattern            # kept as character
  )
}


#' Generate synthetic eGFR dataset for testing
#'
#' Creates a balanced dataset with multiple eGFR trajectories per pattern (linear, quadratic, or segmented),
#' including controlled proportions of missing values.
#'
#' @param pct_miss Numeric vector specifying the percentage of randomly assigned
#'   missing values for each subject (e.g. `c(0, 10, 50, 90)`). The length of
#'   this vector determines how many trajectories/subjects are generated **per pattern**.
#' @param years Time points (default: 0:26)
#' @param baseline_egfr Starting eGFR value (default: 125)
#' @param noise_sd Standard deviation of measurement noise (default: 2)
#' @param seed Random seed for reproducibility (default: NULL = no fixed seed)
#'
#' @return A tibble with columns: id (character), time, egfr, true_pattern (character)
#' @export
#'
#' @examples
#' set.seed(20250116)
#' sim_data <- generate_synthetic_egfr(pct_miss = c(0, 10, 50, 90))
#' str(sim_data)           # id and true_pattern are character
#' table(sim_data$id)      # Number of rows per Id
#' sim_data |> filter(id == 35)
generate_synthetic_egfr <- function(
    pct_miss = c(0, 10, 50, 90),
    years = 0:26,
    baseline_egfr = 125,
    noise_sd = 2,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  patterns <- c("linear", "quadratic", "with_breakpoint")

  n_per_pattern <- length(pct_miss)

  combo <- expand.grid(
    pattern_idx = seq_along(patterns),
    miss_idx    = seq_along(pct_miss),
    stringsAsFactors = FALSE
  )

  # ID format: pattern number (1–3) + miss level / 10
  combo$id_base <- sprintf("%d%d", combo$pattern_idx, pct_miss[combo$miss_idx] / 10)

  trajectories <- purrr::pmap_dfr(
    .l = list(
      id       = combo$id_base,
      pattern  = patterns[combo$pattern_idx],
      pct_miss = pct_miss[combo$miss_idx]
    ),
    .f = function(id, pattern, pct_miss) {

      traj <- generate_one_trajectory(
        id          = id,
        pattern     = pattern,
        years       = years,
        baseline_egfr = baseline_egfr,
        noise_sd    = noise_sd
      )

      # Apply missingness
      n <- nrow(traj)
      if (pct_miss > 0 && pct_miss <= 100) {
        n_miss <- round(n * pct_miss / 100)
        if (n_miss > 0) {
          miss_idx <- sample(seq_len(n), size = n_miss, replace = FALSE)
          traj$egfr[miss_idx] <- NA_real_
        }
      }

      # If ALL values became NA → put one realistic value back
      if (all(is.na(traj$egfr))) {
        one_time <- sample(years, size = 1)
        one_idx  <- which(traj$time == one_time)
        traj$egfr[one_idx] <- baseline_egfr - 1.5 * one_time + rnorm(1, sd = noise_sd * 1.3)
      }

      traj
    }
  )

  # Final touch-ups – id and true_pattern stay character
  trajectories |>
    dplyr::arrange(id, time)
}