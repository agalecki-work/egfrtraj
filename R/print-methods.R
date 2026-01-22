# print-methods.R
# Custom print and summary methods for egfr_traj objects
#' @export
print.egfr_traj <- function(x, ...) {
  cat("<egfr_traj>\n")
  cat("  ID:          ", x$id, "\n")
  cat("  Pattern:     ", x$models$info[[1]]$pattern, "\n")
  cat("  Observations:", nrow(x$trajectory), "\n")
  if (!is.null(x$add_cols)) {
    cat("  Add. columns:", paste(setdiff(names(x$add_cols), "id"), collapse = ", "), "\n")
  }
  invisible(x)
}


#' Summarise a single eGFR trajectory object
#'
#' @param object An `egfr_traj` object
#' @param ... Additional arguments (currently ignored)
#'
#' @return A one-row tibble with key summary statistics
#' @export
summary.egfr_traj <- function(object, ...) {
  info <- object$models$info[[1]]
  traj <- object$trajectory

  # Safe time range
  time_min <- if (nrow(traj) > 0 && any(!is.na(traj$time))) min(traj$time, na.rm = TRUE) else NA_integer_
  time_max <- if (nrow(traj) > 0 && any(!is.na(traj$time))) max(traj$time, na.rm = TRUE) else NA_integer_
  time_str <- if (!is.na(time_min) && !is.na(time_max)) {
    sprintf("%d – %d", time_min, time_max)
  } else "NA – NA"

  # Safe eGFR values
  egfr_start  <- if (nrow(traj) > 0 && !is.na(traj$egfr[1])) round(traj$egfr[1], 1) else NA_real_
  egfr_end    <- if (nrow(traj) > 0 && !is.na(traj$egfr[nrow(traj)])) round(traj$egfr[nrow(traj)], 1) else NA_real_
  egfr_change <- if (nrow(traj) > 1 && !any(is.na(traj$egfr[c(1, nrow(traj))]))) {
                    round(egfr_start - egfr_end, 1)
                  } else NA_real_

  tibble::tibble(
    id              = object$id,
    pattern         = info$pattern,
    n_valid         = info$nobs,
    n_missing       = info$nmiss,
    time_range      = time_str,
    egfr_start      = egfr_start,
    egfr_end        = egfr_end,
    egfr_change     = egfr_change,
    breakpoint_time = if (!is.na(info$breakpoint_time)) round(info$breakpoint_time, 1) else NA_real_,
    breakpoint_se   = if (!is.na(info$bp_se)) round(info$bp_se, 2) else NA_real_
  )
}


#' Summarise results from classify_multiple_trajectories
#'
#' @param object Named list of `egfr_traj` objects
#' @param ... Additional arguments (currently ignored)
#'
#' @return A tibble with one row per subject
#' @export
summary.egfr_traj_list <- function(object, ...) {
  purrr::map_dfr(object, summary.egfr_traj)
}


