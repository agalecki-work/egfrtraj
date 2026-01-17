# plot-methods.R
# Visualization methods for egfr_trajectories objects

#' Plot method for egfr_trajectories objects
#'
#' Creates a faceted line plot of all individual eGFR trajectories,
#' colored by assigned pattern.
#'
#' @param x An object of class `"egfr_trajectories"`
#' @param data Optional original long-format data frame (if not stored in object)
#' @param ncol Number of columns in facet grid (default: 5)
#' @param scales Should facet scales be "free_y", "fixed", etc.? (default: "free_y")
#' @param show_breakpoints Logical. Draw vertical lines at detected breakpoints? (default: TRUE)
#' @param alpha Transparency of lines (default: 0.7)
#' @param point_size Size of individual points (default: 1.2)
#' @param ... Additional arguments passed to ggplot2::ggplot()
#'
#' @return A ggplot object (invisibly returned)
#' @export
#' @method plot egfr_trajectories
plot.egfr_trajectories <- function(
    x,
    data = NULL,
    ncol = 5,
    scales = "free_y",
    show_breakpoints = TRUE,
    alpha = 0.7,
    point_size = 1.2,
    ...
) {
  # Use stored results; fall back to provided data if needed
  if (is.null(data)) {
    if (!exists("data", envir = parent.frame())) {
      stop("No original data found. Provide 'data' argument or ensure data is available.")
    }
    data <- get("data", envir = parent.frame())  # fallback
  }
  
  # Join classification with original data
  plot_data <- data %>%
    dplyr::left_join(x$results, by = "id") %>%
    dplyr::mutate(
      pattern = forcats::fct_relevel(
        pattern,
        "stable", "slow_linear", "rapid_with_breakpoint", "insufficient_data"
      )
    )
  
  # Base plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = time, y = egfr, group = id, color = pattern)
  ) +
    ggplot2::geom_line(alpha = alpha, linewidth = 0.6) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::scale_color_manual(
      values = c(
        "stable"               = "#2ca02c",        # green
        "slow_linear"          = "#1f77b4",        # blue
        "rapid_with_breakpoint" = "#d62728",       # red
        "insufficient_data"    = "#7f7f7f",        # gray
        "NA"                   = "#9467bd"         # purple for any NA
      ),
      na.value = "#9467bd",
      drop = FALSE
    ) +
    ggplot2::labs(
      title    = "Individual eGFR Trajectories by Classified Pattern",
      x        = "Time (years)",
      y        = "eGFR (mL/min/1.73 m²)",
      color    = "Pattern"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(face = "bold"),
      strip.text      = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::facet_wrap(~id, ncol = ncol, scales = scales)
  
  # Optional: add vertical lines for breakpoints
  if (show_breakpoints) {
    bp_data <- x$results %>%
      dplyr::filter(pattern == "rapid_with_breakpoint" & !is.na(breakpoint_time))
    
    if (nrow(bp_data) > 0) {
      p <- p + ggplot2::geom_vline(
        data = bp_data,
        ggplot2::aes(xintercept = breakpoint_time),
        linetype = "dashed", color = "#d62728", linewidth = 0.8,
        alpha = 0.6
      )
    }
  }
  
  print(p)
  invisible(p)
}


#' Autoplot method (ggplot2 compatible)
#'
#' Allows ggplot2-style modification after plotting.
#'
#' @param object An `egfr_trajectories` object
#' @param ... Arguments passed to `plot.egfr_trajectories()`
#'
#' @return A ggplot object
#' @export
#' @method autoplot egfr_trajectories
autoplot.egfr_trajectories <- function(object, ...) {
  plot.egfr_trajectories(object, ...)
}
