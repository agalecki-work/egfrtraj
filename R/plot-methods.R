# plot-methods-tmp.R
# Visualization methods for egfr_trajectories objects
#' Plot an eGFR trajectory
#'
#' @param x An `egfr_traj` object
#' @param fitted Character string: `"best"` (default), `"lin"`, `"quad"`, `"seg"`, or `NULL` (observed only).
#' @param show_ci Logical. Show CI ribbon for the fitted model? Default `TRUE`.
#' @param ... Passed to ggplot2
#'
#' @return ggplot object
#' @export
plot.egfr_traj <- function(x, 
                           fitted = "best",
                           show_ci = TRUE,
                           ...) {
  
  info <- x$models$info[[1]]
  traj <- x$trajectory
  
  # ── Insufficient data ──────────────────────────────────────────────────────
  if (info$pattern == "insufficient_data") {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 6, color = "grey50") +
        theme_minimal() +
        labs(title = sprintf("ID = %s (n = %d)", x$id, info$nobs)) +
        coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
    )
  }
  
  # ── Selected/best model prefix & info ──────────────────────────────────────
  best_prefix <- switch(
    info$pattern,
    linear    = "lin",
    quadratic = "quad",
    segmented = "seg",
    "lin"
  )
  
  # Best model info for subtitle
  best_info <- if (best_prefix == "lin") {
    slope_val <- x$models$linear[[1]]$slope
    slope_se  <- x$models$linear[[1]]$slope_se
    sprintf("Slope = %.2f%s", 
            slope_val, 
            if (!is.na(slope_se)) sprintf(" (SE = %.2f)", slope_se) else "")
  } else if (best_prefix == "quad") {
    "quadratic model"
  } else if (best_prefix == "seg") {
    bp_val <- info$breakpoint_time
    bp_se  <- info$bp_se
    sprintf("Breakpoint = %.1f%s", 
            bp_val, 
            if (!is.na(bp_se)) sprintf(" (SE = %.2f)", bp_se) else "")
  } else {
    "selected model"
  }
  
  best_subtitle <- sprintf("Best model: %s. %s", info$pattern, best_info)
  
  # Delta for observed-only
  delta_txt <- if (nrow(traj) > 1 && !any(is.na(traj$egfr[c(1, nrow(traj))]))) {
    delta <- traj$egfr[nrow(traj)] - traj$egfr[1]
    sprintf("Δ = %.1f", delta)
  } else "Δ = NA"
  
  # ── Determine plotted model ────────────────────────────────────────────────
  if (is.null(fitted) || fitted == "") {
    plotted_prefix <- NULL
  } else if (fitted == "best") {
    plotted_prefix <- best_prefix
  } else {
    plotted_prefix <- match.arg(fitted, c("lin", "quad", "seg"))
  }
  
  # ── Fitted model info for subtitle (if different from best) ─────────────────
  fitted_info <- if (!is.null(plotted_prefix) && plotted_prefix != best_prefix) {
    if (plotted_prefix == "lin") {
      slope_val <- x$models$linear[[1]]$slope
      slope_se  <- x$models$linear[[1]]$slope_se
      sprintf("Slope = %.2f%s", 
              slope_val, 
              if (!is.na(slope_se)) sprintf(" (SE = %.2f)", slope_se) else "")
    } else if (plotted_prefix == "quad") {
      "quadratic model"
    } else if (plotted_prefix == "seg") {
      bp_val <- info$breakpoint_time
      bp_se  <- info$bp_se
      sprintf("Breakpoint = %.1f%s", 
              bp_val, 
              if (!is.na(bp_se)) sprintf(" (SE = %.2f)", bp_se) else "")
    } else "selected model"
  } else {
    NULL
  }
  
  # ── Base plot ──────────────────────────────────────────────────────────────
  p <- ggplot(traj, aes(x = time, y = egfr)) +
    geom_point(color = "black", size = 2.2, alpha = 0.85, shape = 16) +  # always black dots
    labs(
      title    = sprintf("ID = %s (n = %d)", x$id, info$nobs),
      subtitle = if (!is.null(plotted_prefix)) {
                   if (!is.null(fitted_info)) {
                     paste0(best_subtitle, "\nFitted line shown: ", plotted_prefix, ". ", fitted_info)
                   } else {
                     best_subtitle
                   }
                 } else {
                   delta_txt
                 },
      x = "Time (in years)",
      y = "eGFR (mL/min/1.73 m²)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 11)
    )
  
  # ── Add fitted line ────────────────────────────────────────────────────────
  if (!is.null(plotted_prefix)) {
    fitted_col <- paste0(plotted_prefix, "_fitted")
    low_col    <- paste0(plotted_prefix, "_conf.low")
    high_col   <- paste0(plotted_prefix, "_conf.high")
    
    if (fitted_col %in% names(traj)) {
      line_color <- switch(
        plotted_prefix,
        "lin" = "black",
        "quad" = "#7B3294",
        "seg"  = "#CA0020"
      )
      
      is_best <- plotted_prefix == best_prefix
      
      p <- p + 
        geom_line(aes(y = .data[[fitted_col]]), 
                  color = line_color,
                  linewidth = if(is_best) 1.35 else 0.95,
                  linetype = if(is_best) "solid" else "dashed")
      
      # CI ribbon
      if (show_ci && all(c(low_col, high_col) %in% names(traj))) {
        p <- p + 
          geom_ribbon(aes(ymin = .data[[low_col]], ymax = .data[[high_col]]),
                      alpha = 0.12, fill = line_color)
      }
    }
  }
  
  # ── Breakpoint line (if segmented is plotted) ──────────────────────────────
  if (plotted_prefix == "seg" && !is.na(info$breakpoint_time)) {
    p <- p + 
      geom_vline(xintercept = info$breakpoint_time, 
                 linetype = "dotted", color = "darkgreen", linewidth = 0.9) +
      annotate("text", x = info$breakpoint_time + 1, 
               y = max(traj$egfr, na.rm = TRUE) * 0.92,
               label = sprintf("Breakpoint: %.1f", info$breakpoint_time),
               color = "darkgreen", hjust = 0, size = 3.8)
  }
  
  p
}
