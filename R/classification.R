# classification.R
# Core functions for classifying longitudinal eGFR trajectories
# Implements the full-cohort, unbiased classification approach

#' Classify a single eGFR trajectory
#'
#' Internal function that applies the decision tree to one individual's time series.
#'
#' @param df A data frame with columns `time` (numeric) and `egfr` (numeric)
#' @param slope_p Significance threshold for overall and final slope (default: 0.05)
#' @param quad_p Significance threshold for quadratic term (default: 0.1)
#' @param delta_r2 Minimum R² improvement required for non-linearity (default: 0.05)
#' @param max_knots Maximum number of knots in segmented model (default: 2)
#'
#' @return A one-row tibble with classification results
#' @keywords internal
classify_single_trajectory <- function(
    df,
    slope_p = 0.05,
    quad_p = 0.1,
    delta_r2 = 0.05,
    max_knots = 2
) {
  # Early exit if too few points (should already be filtered)
  if (nrow(df) < 5) {
    return(tibble::tibble(
      pattern = "insufficient_data",
      breakpoint_time = NA_real_,
      n_knots = NA_integer_
    ))
  }

  # 1. Fit linear model & check overall trend
  lm_lin <- lm(egfr ~ time, data = df)
  coef_sum <- coef(summary(lm_lin))
  
  overall_slope <- coef_sum["time", "Estimate"]
  overall_p     <- coef_sum["time", "Pr(>|t|)"]
  
  if (is.na(overall_p) || overall_p > slope_p || overall_slope >= 0) {
    return(tibble::tibble(
      pattern = "stable",
      breakpoint_time = NA_real_,
      n_knots = 0L
    ))
  }

  # 2. Test for non-linearity (quadratic)
  lm_quad <- try(lm(egfr ~ time + I(time^2), data = df), silent = TRUE)
  if (inherits(lm_quad, "try-error")) {
    return(tibble::tibble(
      pattern = "slow_linear",
      breakpoint_time = NA_real_,
      n_knots = 0L
    ))
  }
  
  quad_anova <- try(anova(lm_lin, lm_quad), silent = TRUE)
  quad_p_val <- if (inherits(quad_anova, "try-error")) NA else quad_anova$`Pr(>F)`[2]
  
  r2_lin  <- summary(lm_lin)$r.squared
  r2_quad <- summary(lm_quad)$r.squared
  delta   <- r2_quad - r2_lin
  
  if (is.na(quad_p_val) || quad_p_val > quad_p || delta < delta_r2) {
    return(tibble::tibble(
      pattern = "slow_linear",
      breakpoint_time = NA_real_,
      n_knots = 0L
    ))
  }

  # 3. Segmented regression (breakpoint detection)
  seg_fit <- try(
    segmented::segmented(
      lm_lin,
      seg.Z = ~ time,
      psi = list(time = median(df$time)),
      control = segmented::seg.control(K = max_knots, stop.if = 2)
    ),
    silent = TRUE
  )
  
  if (inherits(seg_fit, "try-error") || is.null(seg_fit$psi)) {
    return(tibble::tibble(
      pattern = "slow_linear",
      breakpoint_time = NA_real_,
      n_knots = 0L
    ))
  }
  
  n_knots <- length(seg_fit$psi[,1])
  slopes <- try(segmented::slope(seg_fit)$time, silent = TRUE)
  
  if (inherits(slopes, "try-error") || nrow(slopes) == 0) {
    return(tibble::tibble(
      pattern = "slow_linear",
      breakpoint_time = NA_real_,
      n_knots = 0L
    ))
  }
  
  # Check final segment
  last_slope_est <- slopes[nrow(slopes), "Estimate"]
  last_slope_p   <- slopes[nrow(slopes), "p-value"]
  
  if (n_knots == 0 || is.na(last_slope_p) || 
      last_slope_p > slope_p || last_slope_est >= 0) {
    return(tibble::tibble(
      pattern = "slow_linear",
      breakpoint_time = NA_real_,
      n_knots = 0L
    ))
  }
  
  # Valid breakpoint found
  first_bp <- seg_fit$psi[1, "Est."]
  
  tibble::tibble(
    pattern         = "rapid_with_breakpoint",
    breakpoint_time = round(first_bp, 1),
    n_knots         = as.integer(n_knots)
  )
}


#' Classify eGFR trajectories for multiple individuals
#'
#' Main function of the package: classifies longitudinal eGFR data
#' for all individuals in the dataset.
#'
#' @param data A data frame or tibble in long format
#' @param id_col Name of the ID column (default: "id")
#' @param time_col Name of the time/year column (default: "year")
#' @param egfr_col Name of the eGFR column (default: "eGFR")
#' @param min_obs Minimum number of valid observations per individual (default: 5)
#' @param slope_p Significance level for slopes (default: 0.05)
#' @param quad_p Significance level for quadratic term (default: 0.1)
#' @param delta_r2 Minimum R² improvement for non-linearity (default: 0.05)
#' @param max_knots Maximum number of breakpoints to allow (default: 2)
#' @param verbose Logical. Print progress/summary messages? (default: TRUE)
#'
#' @return An object of class `egfr_trajectories` with components:
#'   - `results`: tibble with classification for each individual
#'   - `call`: the function call
#'   - `params`: parameters used
#' @export
#' @examples
#' \dontrun{
#'   # Simulated data
#'   sim <- generate_synthetic_egfr(n_per_pattern = 3)
#'   result <- classify_egfr_trajectories(sim)
#'   print(result)
#' }
classify_egfr_trajectories <- function(
    data,
    id_col     = "id",
    time_col   = "year",
    egfr_col   = "eGFR",
    min_obs    = 5,
    slope_p    = 0.05,
    quad_p     = 0.1,
    delta_r2   = 0.05,
    max_knots  = 2,
    verbose    = TRUE
) {
  # Prepare data
  prepared <- prepare_egfr_data(
    data,
    id_col   = id_col,
    time_col = time_col,
    egfr_col = egfr_col,
    min_obs  = min_obs
  )
  
  if (verbose) {
    message("Prepared data for ", n_distinct(prepared$id), " individuals")
  }
  
  # Classify each trajectory
  results <- prepared %>%
    tidyr::nest(data = c(time, egfr)) %>%
    dplyr::mutate(
      classification = purrr::map(
        data,
        ~ classify_single_trajectory(
          .x,
          slope_p    = slope_p,
          quad_p     = quad_p,
          delta_r2   = delta_r2,
          max_knots  = max_knots
        )
      )
    ) %>%
    tidyr::unnest(classification) %>%
    dplyr::select(id, pattern, breakpoint_time, n_knots)
  
  if (verbose) {
    message("Classification complete. Pattern distribution:")
    print(table(results$pattern))
  }
  
  structure(
    list(
      results = results,
      call    = match.call(),
      params  = list(
        min_obs    = min_obs,
        slope_p    = slope_p,
        quad_p     = quad_p,
        delta_r2   = delta_r2,
        max_knots  = max_knots
      )
    ),
    class = "egfr_trajectories"
  )
}


#' Prepare longitudinal eGFR data
#'
#' Standardizes column names and applies basic filtering.
#'
#' @inheritParams classify_egfr_trajectories
#'
#' @return A tibble ready for classification
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