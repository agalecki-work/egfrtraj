# utils.R
# Internal utility functions for the egfr-trajectories package

#' Prepare longitudinal eGFR data
#'
#' Standardizes column names, removes missing values, and filters by minimum observations.
#' Used internally by classify_egfr_trajectories().
#'
#' @param data A data frame or tibble in long format
#' @param id_col Name of the ID column
#' @param time_col Name of the time/year column
#' @param egfr_col Name of the eGFR column
#' @param min_obs Minimum number of valid observations per individual
#'
#' @return A prepared tibble ready for classification
#' @keywords internal
prepare_egfr_data <- function(data, id_col, time_col, egfr_col, min_obs) {
  data %>%
    dplyr::rename(
      id   = dplyr::all_of(id_col),
      time = dplyr::all_of(time_col),
      egfr = dplyr::all_of(egfr_col)
    ) %>%
    dplyr::select(id, time, egfr) %>%
    dplyr::filter(!is.na(egfr), !is.na(time)) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::n() >= min_obs) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = as.factor(id))
}


#' Validate input parameters for classification
#'
#' Checks that user-provided thresholds and parameters are sensible.
#' Throws informative errors if invalid.
#'
#' @param min_obs Minimum observations per ID
#' @param slope_p Slope p-value threshold
#' @param quad_p Quadratic p-value threshold
#' @param delta_r2 Minimum R² improvement
#' @param max_knots Maximum number of knots
#'
#' @return Invisible TRUE if all valid
#' @keywords internal
validate_classification_params <- function(min_obs, slope_p, quad_p, delta_r2, max_knots) {
  if (!is.numeric(min_obs) || min_obs < 3 || min_obs > 50) {
    stop("min_obs must be an integer between 3 and 50", call. = FALSE)
  }
  
  if (!is.numeric(slope_p) || slope_p <= 0 || slope_p > 1) {
    stop("slope_p must be a number between 0 and 1", call. = FALSE)
  }
  
  if (!is.numeric(quad_p) || quad_p <= 0 || quad_p > 1) {
    stop("quad_p must be a number between 0 and 1", call. = FALSE)
  }
  
  if (!is.numeric(delta_r2) || delta_r2 < 0 || delta_r2 > 1) {
    stop("delta_r2 must be a non-negative number ≤ 1", call. = FALSE)
  }
  
  if (!is.numeric(max_knots) || max_knots < 1 || max_knots > 5) {
    stop("max_knots must be an integer between 1 and 5", call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Summarize pattern distribution with proportions
#'
#' Helper for print/summary methods.
#'
#' @param patterns Factor or character vector of assigned patterns
#'
#' @return A tibble with pattern, count, and proportion
#' @keywords internal
summarize_patterns <- function(patterns) {
  tibble::tibble(
    pattern = names(table(patterns, useNA = "ifany")),
    count   = as.integer(table(patterns, useNA = "ifany")),
    percent = round(prop.table(table(patterns, useNA = "ifany")) * 100, 1)
  ) %>%
    dplyr::mutate(percent = paste0(percent, "%"))
}


#' Check if segmented package is available and loaded
#'
#' @return Logical
#' @keywords internal
has_segmented <- function() {
  requireNamespace("segmented", quietly = TRUE)
}


#' Safe segmented fit with fallback
#'
#' Attempts segmented regression and returns NULL on failure.
#'
#' @param lm_fit Linear model object
#' @param max_knots Maximum knots to try
#'
#' @return segmented object or NULL
#' @keywords internal
safe_segmented <- function(lm_fit, max_knots = 2) {
  if (!has_segmented()) {
    return(NULL)
  }
  
  tryCatch(
    segmented::segmented(
      lm_fit,
      seg.Z = ~ time,
      psi = list(time = median(lm_fit$model$time)),
      control = segmented::seg.control(K = max_knots, stop.if = 2)
    ),
    error = function(e) NULL
  )
}
