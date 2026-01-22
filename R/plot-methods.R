# plot-methods.R
# Visualization methods for egfr_trajectories objects

#' @importFrom ggplot2  # autoplot
NULL


#' Plot an eGFR trajectory
#'
#' @param x An `egfr_traj` object (from `classify_single_trajectory()`)
#' @param plot_model Character vector specifying **additional** fitted models to show.
#'   Options: `"lin"`, `"quad"`, `"seg"`. Default: `c("lin", "quad", "seg")` (all).
#'   The **selected model** is **always plotted**, even if not in `plot_model`.
#'   Use `character(0)` or `NULL` to show **only the selected model**.
#' @param show_ci Logical. Show confidence bands for the selected model? Default `TRUE`.
#' @param ... Passed to ggplot2
#'
#' @return A ggplot object
#' @export
plot.egfr_traj <- function(x, 
                           plot_model = c("lin", "quad", "seg"),
                           show_ci = TRUE,
                           ...) {
  info <- x$models$info[[1]]
  traj <- x$trajectory

  # ── Handle insufficient data ───────────────────────────────────────────────
  if (info$pattern == "insufficient_data") {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "Insufficient data\n(no valid eGFR points)", 
               size = 6, color = "grey50") +
      theme_minimal() +
      labs(title = sprintf("Subject %s – Insufficient Data", x$id)) +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), 
            axis.title = element_blank())
    return(p)
  }

  # ── Selected model prefix ──────────────────────────────────────────────────
  selected_prefix <- switch(info$pattern,
                            "linear"    = "lin",
                            "quadratic" = "quad",
                            "segmented" = "seg",
                            "lin")

  # ── Models to plot (always include selected) ───────────────────────────────
  if (length(plot_model) == 0 || is.null(plot_model)) {
    all_models <- selected_prefix
  } else {
    plot_model <- match.arg(plot_model, choices = c("lin", "quad", "seg"), several.ok = TRUE)
    all_models <- unique(c(plot_model, selected_prefix))
  }

  # ── Base plot ──────────────────────────────────────────────────────────────
  p <- ggplot(traj, aes(x = time, y = egfr)) +
    geom_point(color = "steelblue", size = 2, alpha = 0.7) +
    labs(
      title    = sprintf("Subject %s – %s trajectory", x$id, info$pattern),
      subtitle = sprintf("Observations: %d   eGFR change: %.1f", 
                         info$nobs, traj$egfr[1] - traj$egfr[nrow(traj)]),
      x = "Time", y = "eGFR"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # ── Model aesthetics ───────────────────────────────────────────────────────
  model_colors <- c("lin" = "grey60", "quad" = "purple", "seg" = "tomato")
  model_labels <- c("lin" = "Linear", "quad" = "Quadratic", "seg" = "Segmented (selected)")

  # ── Add lines ──────────────────────────────────────────────────────────────
  for (m in all_models) {
    fitted_col <- paste0(m, "_fitted")
    if (fitted_col %in% names(traj)) {
      is_selected <- m == selected_prefix
      p <- p + 
        geom_line(aes(y = .data[[fitted_col]]), 
                  color = model_colors[m],
                  linewidth = if (is_selected) 1.4 else 0.9,
                  linetype = if (is_selected) "solid" else "dashed")
    }
  }

  # ── Unified legend ─────────────────────────────────────────────────────────
  used_colors <- model_colors[match(all_models, names(model_colors))]
  used_labels <- model_labels[match(all_models, names(model_labels))]

  p <- p + 
    scale_color_manual(values = used_colors,
                       labels = used_labels,
                       name = "Model")

  # ── Confidence bands (selected model only) ─────────────────────────────────
  if (show_ci) {
    low_col  <- paste0(selected_prefix, "_conf.low")
    high_col <- paste0(selected_prefix, "_conf.high")
    if (all(c(low_col, high_col) %in% names(traj))) {
      p <- p + 
        geom_ribbon(aes(ymin = .data[[low_col]], ymax = .data[[high_col]]),
                    alpha = 0.15, fill = model_colors[selected_prefix])
    }
  }

  # ── Breakpoint line ────────────────────────────────────────────────────────
  if (info$pattern == "segmented" && !is.na(info$breakpoint_time)) {
    p <- p + 
      geom_vline(xintercept = info$breakpoint_time, 
                 linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
      annotate("text", x = info$breakpoint_time + 1, 
               y = max(traj$egfr, na.rm = TRUE) * 0.9,
               label = sprintf("Breakpoint: %.1f", info$breakpoint_time),
               color = "darkgreen", hjust = 0)
  }

  p
}
