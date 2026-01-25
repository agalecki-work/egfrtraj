
# ──────────────────────────────────────────────────────────────────────────────
# Helper functions (all internal)
# ──────────────────────────────────────────────────────────────────────────────

.plot_insufficient_data <- function(x, info) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 6, color = "grey50") +
    theme_minimal() +
    labs(title = sprintf("ID = %s (n = %d)", x$id, info$nobs)) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
}

.get_selected_prefix <- function(pattern) {
  switch(
    pattern,
    linear    = "lin",
    quadratic = "quad",
    segmented = "seg",
    "lin"
  )
}

.extract_slope_info <- function(models, selected_prefix) {
  if (selected_prefix == "lin") {
    val <- models$linear[[1]]$slope
    se  <- models$linear[[1]]$slope_se
  } else if (selected_prefix == "quad") {
    val <- models$quadratic[[1]]$linear
    se  <- models$quadratic[[1]]$linear_se
  } else if (selected_prefix == "seg") {
    val <- models$segmented[[1]]$slope_post
    se  <- models$segmented[[1]]$slope_post_se
  } else {
    val <- NA_real_
    se  <- NA_real_
  }
  
  txt <- if (!is.na(val)) {
    if (!is.na(se)) sprintf("%.2f (SE = %.2f)", val, se) else sprintf("%.2f", val)
  } else "NA"
  
  list(val = val, se = se, txt = txt)
}

.process_fitted_argument <- function(fitted, selected_prefix) {
  if (is.null(fitted) || length(fitted) == 0) {
    return(character(0))
  }
  
  fitted <- gsub("^(fitted|best)$", selected_prefix, fitted)
  models_to_plot <- unique(fitted)
  
  invalid <- setdiff(models_to_plot, c("lin", "quad", "seg"))
  if (length(invalid) > 0) {
    warning("Ignoring unknown model names: ", paste(invalid, collapse = ", "))
    models_to_plot <- intersect(models_to_plot, c("lin", "quad", "seg"))
  }
  
  models_to_plot
}

.determine_models_to_plot <- function(models_to_plot, selected_prefix) {
  if (length(models_to_plot) > 0) {
    unique(c(models_to_plot, selected_prefix))
  } else {
    character(0)
  }
}

.base_ggplot <- function(traj, id, nobs) {
  ggplot(traj, aes(x = time, y = egfr)) +
    geom_point(aes(color = "Observed"), size = 2.2, alpha = 0.85, shape = 16) +
    labs(
      title = sprintf("ID = %s (n = %d)", id, nobs),
      x = "Time", y = "eGFR"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 11)
    )
}

.get_model_color <- function(model_prefix) {
  model_colors <- c(
    lin  = "grey50",
    quad = "#7B3294",
    seg  = "#CA0020"
  )
  model_colors[[model_prefix]] %||% "black"
}

.build_legend_breaks <- function(all_to_plot) {
  c("Observed", all_to_plot)
}

.build_legend_labels <- function(all_to_plot) {
  labels = c("Observed", c(lin = "Linear", quad = "Quadratic", seg = "Segmented")[all_to_plot])
  unname(labels)
}

.add_model_lines <- function(p, traj, all_to_plot, selected_prefix) {
  for (m in all_to_plot) {
    fitted_col <- paste0(m, "_fitted")
    if (fitted_col %in% names(traj)) {
      is_selected <- m == selected_prefix
      line_color <- .get_model_color(m)
      p <- p + 
        geom_line(aes(y = .data[[fitted_col]]), 
                  color = line_color,
                  linewidth = if(is_selected) 1.35 else 0.95,
                  linetype = if(is_selected) "solid" else "longdash")
    }
  }
  p
}

.add_legend <- function(p, all_to_plot) {
  breaks <- .build_legend_breaks(all_to_plot)
  labels <- .build_legend_labels(all_to_plot)
  
  p + 
    scale_color_identity(
      name = NULL,
      breaks = breaks,
      labels = labels,
      guide = guide_legend(override.aes = list(
        linetype = c("blank", rep("solid", length(all_to_plot))),
        shape = c(16, rep(NA, length(all_to_plot)))
      ))
    )
}

.add_ci_ribbon <- function(p, traj, first_model) {
  low  <- paste0(first_model, "_conf.low")
  high <- paste0(first_model, "_conf.high")
  
  if (all(c(low, high) %in% names(traj))) {
    fill_color <- .get_model_color(first_model)
    p + 
      geom_ribbon(aes(ymin = .data[[low]], ymax = .data[[high]]),
                  alpha = 0.12, fill = fill_color)
  } else {
    p
  }
}

.add_breakpoint_annotation <- function(p, bp_time, max_egfr) {
  p + 
    geom_vline(xintercept = bp_time, linetype = "dotted", color = "darkgreen", linewidth = 0.9) +
    annotate("text", x = bp_time + 1, 
             y = max_egfr * 0.92,
             label = sprintf("Breakpoint: %.1f", bp_time),
             color = "darkgreen", hjust = 0, size = 3.8)
}

