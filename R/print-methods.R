# print-methods.R
# Custom print and summary methods for egfr_trajectories objects

#' Print method for egfr_trajectories objects
#'
#' Displays a concise, human-readable overview of the classification results.
#'
#' @param x An object of class `"egfr_trajectories"`
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns `x`
#' @export
#' @method print egfr_trajectories
print.egfr_trajectories <- function(x, ...) {
  cat("eGFR Trajectory Classification Results\n")
  cat("────────────────────────────────────\n")
  
  # Basic info
  n_total <- nrow(x$results)
  n_classified <- sum(!x$results$pattern %in% c("insufficient_data", NA))
  cat("Total individuals processed: ", n_total, "\n")
  cat("Successfully classified:     ", n_classified, " (", 
      round(n_classified / n_total * 100, 1), "%)\n\n")
  
  # Pattern distribution
  cat("Pattern Distribution:\n")
  pat_tab <- table(x$results$pattern, useNA = "ifany")
  pat_prop <- round(prop.table(pat_tab) * 100, 1)
  
  print_df <- data.frame(
    Pattern = names(pat_tab),
    Count   = as.integer(pat_tab),
    Percent = paste0(pat_prop, "%")
  )
  print(print_df, row.names = FALSE, right = FALSE)
  cat("\n")
  
  # Breakpoint summary (only for rapid_with_breakpoint)
  rapid <- x$results %>% dplyr::filter(pattern == "rapid_with_breakpoint")
  if (nrow(rapid) > 0) {
    cat("Breakpoint Timing Summary (rapid_with_breakpoint):\n")
    cat("  - Number of cases: ", nrow(rapid), "\n")
    cat("  - Mean breakpoint time: ", 
        round(mean(rapid$breakpoint_time, na.rm = TRUE), 1), " years\n")
    cat("  - Median breakpoint time: ", 
        round(median(rapid$breakpoint_time, na.rm = TRUE), 1), " years\n")
    cat("  - Range: ", 
        round(min(rapid$breakpoint_time, na.rm = TRUE), 1), "–",
        round(max(rapid$breakpoint_time, na.rm = TRUE), 1), " years\n\n")
  }
  
  cat("Use `summary()` for more details or `plot()` for visualization.\n")
  invisible(x)
}


#' Summary method for egfr_trajectories objects
#'
#' Provides a more detailed summary than the default print method.
#'
#' @param object An object of class `"egfr_trajectories"`
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns a list with summary statistics
#' @export
#' @method summary egfr_trajectories
summary.egfr_trajectories <- function(object, ...) {
  res <- object$results
  
  # Pattern counts
  pat_counts <- table(res$pattern, useNA = "ifany")
  
  # Breakpoint statistics
  rapid <- res %>% dplyr::filter(pattern == "rapid_with_breakpoint")
  bp_stats <- if (nrow(rapid) > 0) {
    list(
      n = nrow(rapid),
      mean_bp = mean(rapid$breakpoint_time, na.rm = TRUE),
      median_bp = median(rapid$breakpoint_time, na.rm = TRUE),
      min_bp = min(rapid$breakpoint_time, na.rm = TRUE),
      max_bp = max(rapid$breakpoint_time, na.rm = TRUE),
      n_knots_mean = mean(rapid$n_knots, na.rm = TRUE)
    )
  } else {
    list(n = 0, mean_bp = NA, median_bp = NA, min_bp = NA, max_bp = NA, n_knots_mean = NA)
  }
  
  # Used parameters
  params <- object$params
  
  structure(
    list(
      n_total = nrow(res),
      pattern_counts = pat_counts,
      breakpoint_stats = bp_stats,
      parameters = params,
      call = object$call
    ),
    class = "summary.egfr_trajectories"
  )
}


#' Print method for summary.egfr_trajectories
#'
#' @param x A summary object from `summary.egfr_trajectories()`
#' @param digits Number of decimal places (default: 2)
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method print summary.egfr_trajectories
print.summary.egfr_trajectories <- function(x, digits = 2, ...) {
  cat("Summary of eGFR Trajectory Classification\n")
  cat("────────────────────────────────────────\n\n")
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  cat("Parameters used:\n")
  cat("  Min observations per ID: ", x$parameters$min_obs, "\n")
  cat("  Slope significance level: ", x$parameters$slope_p, "\n")
  cat("  Quadratic p-value threshold: ", x$parameters$quad_p, "\n")
  cat("  Min ΔR² for non-linearity: ", x$parameters$delta_r2, "\n")
  cat("  Max knots allowed: ", x$parameters$max_knots, "\n\n")
  
  cat("Classification Results:\n")
  cat("  Total trajectories: ", x$n_total, "\n\n")
  
  cat("Pattern Counts:\n")
  print(as.data.frame(x$pattern_counts), right = FALSE)
  cat("\n")
  
  cat("Breakpoint Statistics (rapid_with_breakpoint):\n")
  if (x$breakpoint_stats$n > 0) {
    cat("  Number of rapid cases: ", x$breakpoint_stats$n, "\n")
    cat("  Mean breakpoint time:  ", round(x$breakpoint_stats$mean_bp, digits), " years\n")
    cat("  Median breakpoint time:", round(x$breakpoint_stats$median_bp, digits), " years\n")
    cat("  Range:                 ", round(x$breakpoint_stats$min_bp, digits), "–",
        round(x$breakpoint_stats$max_bp, digits), " years\n")
    cat("  Average knots per case:", round(x$breakpoint_stats$n_knots_mean, digits), "\n")
  } else {
    cat("  No rapid_with_breakpoint trajectories detected.\n")
  }
  
  invisible(x)
}
