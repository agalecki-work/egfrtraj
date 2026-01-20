# plot-methods.R
# Visualization methods for egfr_trajectories objects

#' @importFrom ggplot2 autoplot
NULL


#' Plot method for egfrtraj_single objects
#'
#' Plots the observed eGFR trajectory, fitted values from the selected model,
#' 95% prediction interval ribbon, and breakpoint line (if applicable).
#'
#' @param x An object of class "egfrtraj_single" (from classify_single_trajectory_ic())
#' @param data Ignored (uses x$data internally)
#' @param ncol Number of columns in facet (default: 1, single subject)
#' @param scales Y-scale behavior ("fixed" or "free_y", default: "fixed")
#' @param show_breakpoints Logical. Draw vertical line at breakpoint? (default: TRUE)
#' @param point_size Size of observed points (default: 2)
#' @param line_width Width of fitted line (default: 1.2)
#' @param ribbon_alpha Transparency of prediction interval ribbon (default: 0.2)
#' @param ... Additional arguments passed to ggplot2::ggplot()
#'
#' @return A ggplot object (invisibly returned)
#' @export
#' @method plot egfrtraj_single
plot.egfrtraj_single <- function(
    x,
    data = NULL,  # ignored - uses x$data
    ncol = 1,
    scales = "fixed",
    show_breakpoints = TRUE,
    point_size = 2,
    line_width = 1.2,
    ribbon_alpha = 0.2,
    ...
) {
  # Extract data with predictions
  plot_data <- x$data

  # Base plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = time, y = egfr)
  ) +
    ggplot2::geom_point(size = point_size, color = "darkblue", alpha = 0.7) +
    ggplot2::geom_line(
      ggplot2::aes(y = fitted),
      color = "blue",
      linewidth = line_width
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
      alpha = ribbon_alpha,
      fill = "blue"
    ) +
    ggplot2::labs(
      title    = paste("eGFR Trajectory - Selected Model:", x$selected_pattern),
      subtitle = paste("BIC winner:", x$selected_model_name),
      x        = "Time (years)",
      y        = "eGFR (mL/min/1.73 m²)"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  # Add breakpoint line if present and selected
  if (show_breakpoints && !is.na(x$breakpoint_time) && x$selected_model_name == "segmented") {
    p <- p + ggplot2::geom_vline(
      xintercept = x$breakpoint_time,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    )
  }

  # Apply scales & facets (single subject → ncol=1 usually)
  p <- p + ggplot2::facet_wrap(~ id, scales = scales, ncol = ncol)

  print(p)
  invisible(p)
}

#' Autoplot method for egfrtraj_single objects
#'
#' Returns a ggplot object of the eGFR trajectory with selected model fit,
#' prediction intervals, and breakpoint line (if applicable).
#'
#' @param object An egfrtraj_single object
#' @param ... Additional arguments passed to ggplot2::ggplot()
#'
#' @return A ggplot object
#' @export
#' @method autoplot egfrtraj_single
autoplot.egfrtraj_single <- function(object, ...) {
  plot_data <- object$data

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = time, y = egfr)
  ) +
    ggplot2::geom_point(size = 2, color = "darkblue", alpha = 0.7) +
    ggplot2::geom_line(
      ggplot2::aes(y = fitted),
      color = "blue",
      linewidth = 1.2
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
      alpha = 0.2,
      fill = "blue"
    ) +
    ggplot2::labs(
      title    = paste("eGFR Trajectory - Selected:", object$selected_pattern),
      subtitle = paste("Model:", object$selected_model_name),
      x        = "Time (years)",
      y        = "eGFR (mL/min/1.73 m²)"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  # Add breakpoint line if present
  if (!is.na(object$breakpoint_time) && object$selected_model_name == "segmented") {
    p <- p + ggplot2::geom_vline(
      xintercept = object$breakpoint_time,
      linetype = "dashed",
      color = "red",
      linewidth = 1
    )
  }

  p
}


#' Plot method for egfrtraj_single objects
#'
#' Convenience wrapper that calls autoplot() and prints the result.
#'
#' @param x An egfrtraj_single object
#' @param ... Passed to autoplot.egfrtraj_single()
#'
#' @export
#' @method plot egfrtraj_single
plot.egfrtraj_single <- function(x, ...) {
  p <- autoplot.egfrtraj_single(x, ...)
  print(p)
  invisible(p)
}

