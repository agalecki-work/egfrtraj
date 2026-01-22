# classify_single_trajectory
#' Classify a single eGFR trajectory
#'
#' @param df Data frame with data for exactly one subject
#' @param id_col,time_col,egfr_col Column names
#' @param bic_tie_threshold Default 4
#' @param ci_level Default 0.95
#'
#' @return List with id, trajectory, models
#' @export
classify_single_trajectory <- function(
    df,
    id_col             = "id",
    time_col           = "time",
    egfr_col           = "egfr",
    bic_tie_threshold  = 4,
    ci_level           = 0.95
) {
  # ── Early filtering & renaming ─────────────────────────────────────────────
  df_clean <- df %>%
    dplyr::select(
      id   = all_of(id_col),
      time = all_of(time_col),
      egfr = all_of(egfr_col)
    ) %>%
    dplyr::filter(!is.na(time), !is.na(egfr))

  n_valid    <- nrow(df_clean)
  n_miss     <- nrow(df) - n_valid

  # Safe subject_id even if df_clean is empty
  subject_id <- if (nrow(df) > 0) unique(df[[id_col]])[1] else NA_character_

  # ── EARLY EXIT ─ must be here, before any other logic ──────────────────────
 

 if (n_valid < 5L) {
  return(list(
    id         = subject_id,
    trajectory = tibble::tibble(
      id   = subject_id,
      time = NA_integer_,
      egfr = NA_real_
    ),
    models     = tibble::tibble(
      info = list(list(
        pattern         = "insufficient_data",
        nobs            = n_valid,
        nmiss           = n_miss,
        breakpoint_time = NA_real_,
        bp_se           = NA_real_
      )),
      linear    = list(list(
        intercept = NA_real_,
        slope     = NA_real_,
        note      = "Insufficient data - no model fitted"
      )),
      quadratic = list(list(
        intercept = NA_real_,
        linear    = NA_real_,
        quadratic = NA_real_,
        note      = "Insufficient data - no model fitted"
      )),
      segmented = list(list(
        intercept_pre = NA_real_,
        slope_pre     = NA_real_,
        slope_post    = NA_real_,
        note          = "Insufficient data - no model fitted"
      ))
    )
  ) |> structure(class = "egfr_traj")
)
}

  # ── Only reached if n_valid >= 5 ───────────────────────────────────────────
  trajectory <- df_clean   # already has id, time, egfr
    nested_info <- list(
      pattern         = "insufficient_data",   # temporary
      nobs            = n_valid,
      nmiss           = n_miss,
      breakpoint_time = NA_real_,
      bp_se           = NA_real_
    )

  
  # ── Model fitting ──────────────────────────────────────────────────────────
  lm_lin <- lm(egfr ~ time, data = df_clean)

  lm_quad <- try(lm(egfr ~ time + I(time^2), data = df_clean), silent = TRUE)

  seg_fit <- try(
    segmented::segmented(
      lm_lin,
      seg.Z = ~ time,
      psi   = median(df_clean$time),
      control = segmented::seg.control(K = 1L, n.boot = 10L, it.max = 30L)
    ),
    silent = TRUE
  )

  models_list <- list(
    linear    = lm_lin,
    quadratic = if (!inherits(lm_quad, "try-error")) lm_quad else NULL,
    segmented = if (!inherits(seg_fit, "try-error") && !is.null(seg_fit$psi)) seg_fit else NULL
  )

  aic_vals <- purrr::map_dbl(models_list, ~ if (!is.null(.x)) AIC(.x) else NA_real_)
  bic_vals <- purrr::map_dbl(models_list, ~ if (!is.null(.x)) BIC(.x) else NA_real_)

  # Select best
  valid_models <- names(bic_vals)[!is.na(bic_vals)]
  best_pattern <- if (length(valid_models) == 0) "linear" else {
    min_bic <- min(bic_vals[valid_models])
    winner <- valid_models[which.min(bic_vals[valid_models])]
    if ("segmented" %in% valid_models && bic_vals["segmented"] <= min_bic + bic_tie_threshold) {
      "segmented"
    } else winner
  }

  bp_time <- NA_real_
  bp_se   <- NA_real_
  if (best_pattern == "segmented" && !is.null(models_list$segmented$psi)) {
    bp_time <- round(models_list$segmented$psi[1, "Est."], 1)
    bp_se   <- models_list$segmented$psi[1, "St.Err"]
  }

  # Update info
  nested_info <- list(list(
    pattern         = best_pattern,
    nobs            = n_valid,
    nmiss           = n_miss,
    breakpoint_time = bp_time,
    bp_se           = bp_se
  ))

  # Predictions
  pred_lin <- suppressWarnings(predict(lm_lin, interval = "prediction", level = ci_level))
  pred_quad <- if (!is.null(models_list$quadratic))
    suppressWarnings(predict(models_list$quadratic, interval = "prediction", level = ci_level)) else NULL
  pred_seg  <- if (!is.null(models_list$segmented))
    suppressWarnings(predict(models_list$segmented, interval = "prediction", level = ci_level)) else NULL

  trajectory <- trajectory |>
    dplyr::mutate(
      lin_fitted     = pred_lin[, "fit"],
      lin_conf.low   = pred_lin[, "lwr"],
      lin_conf.high  = pred_lin[, "upr"],
      quad_fitted    = if (!is.null(pred_quad)) pred_quad[, "fit"]     else NA_real_,
      quad_conf.low  = if (!is.null(pred_quad)) pred_quad[, "lwr"]     else NA_real_,
      quad_conf.high = if (!is.null(pred_quad)) pred_quad[, "upr"]     else NA_real_,
      seg_fitted     = if (!is.null(pred_seg))  pred_seg[, "fit"]      else NA_real_,
      seg_conf.low   = if (!is.null(pred_seg))  pred_seg[, "lwr"]      else NA_real_,
      seg_conf.high  = if (!is.null(pred_seg))  pred_seg[, "upr"]      else NA_real_
    )

  # Model coefficients
  extract_coefs <- function(m, name) {
    if (is.null(m)) return(list())
    s <- summary(m)
    cf <- coef(m)
    if (name == "linear") {
      list(intercept = cf[1], slope = cf[2],
           intercept_se = s$coefficients[1,2], slope_se = s$coefficients[2,2],
           AIC = aic_vals[[name]], BIC = bic_vals[[name]])
    } else if (name == "quadratic") {
      list(intercept = cf[1], linear = cf[2], quadratic = cf[3],
           intercept_se = s$coefficients[1,2], linear_se = s$coefficients[2,2],
           quadratic_se = s$coefficients[3,2],
           AIC = aic_vals[[name]], BIC = bic_vals[[name]])
    } else if (name == "segmented") {
      if (is.null(m$psi) || nrow(m$psi) == 0) return(list())
      list(intercept_pre = cf[1], slope_pre = cf[2], slope_post = cf[3],
           intercept_se = s$coefficients[1,2], slope_pre_se = s$coefficients[2,2],
           slope_post_se = s$coefficients[3,2],
           AIC = aic_vals[[name]], BIC = bic_vals[[name]])
    } else list()
  }

  models <- tibble::tibble(
    info      = nested_info,
    linear    = list(extract_coefs(models_list$linear,    "linear")),
    quadratic = list(extract_coefs(models_list$quadratic, "quadratic")),
    segmented = list(extract_coefs(models_list$segmented, "segmented"))
  )

  list(
    id         = subject_id,
    trajectory = trajectory,
    models     = models 
  ) |> structure(class = "egfr_traj")
}

# =============================================================================
# classify_multiple_trajectories
# =============================================================================

#' Classify multiple eGFR trajectories
#'
#' @inheritParams classify_single_trajectory
#' @param data Data frame with multiple (or single) subjects
#' @param best_model_only If `TRUE` (default), keep only selected model results
#' @param estimates_only If `TRUE`, drop standard error columns
#' @param add_cols Character vector of additional columns to return separately
#'
#' @return Named list of `egfr_traj` objects (one per subject)
#' @export
classify_multiple_trajectories <- function(
    data,
    id_col             = "id",
    time_col           = "time",
    egfr_col           = "egfr",
    add_cols           = NULL,
    bic_tie_threshold  = 4,
    ci_level           = 0.95,
    best_model_only    = TRUE,
    estimates_only     = FALSE
) {
  subjects <- unique(data[[id_col]])
  if (length(subjects) == 0) stop("No subjects found")

  result <- setNames(vector("list", length(subjects)), as.character(subjects))

  for (subj in subjects) {
    subj_data <- data[data[[id_col]] == subj, , drop = FALSE]

    single <- classify_single_trajectory(
      df                = subj_data,
      id_col            = id_col,
      time_col          = time_col,
      egfr_col          = egfr_col,
      bic_tie_threshold = bic_tie_threshold,
      ci_level          = ci_level
    )
    
   # cat("Subject:", subj, "| info length:", length(single$models$info), "\n")
   # if (length(single$models$info) > 0) {
   #  cat("  info[[1]] exists:", !is.null(single$models$info[[1]]), "\n")
   #}


    # pattern <- single$models$info[[1]]$pattern
   pattern <- if (length(single$models$info) > 0 && length(single$models$info[[1]]) > 0) {
      single$models$info[[1]]$pattern
   } else {
      "insufficient_data"
   }

    # ── Trajectory: keep only selected model ───────────────────────────────
    
   if (pattern == "insufficient_data") {
     # No fitted values exist – keep only core columns
     traj <- single$trajectory |>
    dplyr::select(id, time, egfr)
    } else {
    # Normal case: select and rename the selected model's columns
    prefix <- switch(pattern,
      "linear"    = "lin",
      "quadratic" = "quad",
      "segmented" = "seg",
      "lin"  # fallback
    )

    keep_cols <- c("id", "time", "egfr",
                   paste0(prefix, c("_fitted", "_conf.low", "_conf.high")))

    traj <- single$trajectory |>
      dplyr::select(any_of(keep_cols)) |>   # any_of() is safer
      dplyr::rename_with(
        ~ gsub(paste0("^", prefix, "_"), "", .x),
        .cols = starts_with(prefix)           # now safe inside rename_with
      ) |>
      dplyr::rename(fitted_value = fitted)
  } # if (pattern 

    # ── Models: simplify if requested ──────────────────────────────────────
    models_s <- single$models

    if (best_model_only) {
      # Select info + columns starting with the prefix
      models_s <- models_s |>
        dplyr::select(
          info,
          dplyr::starts_with(paste0(prefix, "_"))
        ) |>
        # Rename the selected columns by removing the prefix
        dplyr::rename_with(
          ~ gsub(paste0("^", prefix, "_"), "", .x),
          .cols = starts_with(paste0(prefix, "_"))
        )
    }

    if (estimates_only) {
      models_s <- models_s |>
        dplyr::select(-dplyr::ends_with("_se"))
    }

    out <- list(
      id         = subj,
      trajectory = traj,
      models     = models_s
    )

    # ── Optional add_cols ──────────────────────────────────────────────────
    if (!is.null(add_cols) && length(add_cols) > 0) {
      add_df <- data |>
        dplyr::filter(.data[[id_col]] == subj) |>
        dplyr::select(all_of(c(id_col, add_cols))) |>
        dplyr::distinct()
      if (id_col != "id") add_df <- dplyr::rename(add_df, id = all_of(id_col))
      out$add_cols <- add_df
    }

    class(out) <- "egfr_traj_collection"
    result[[as.character(subj)]] <- out
  }

  result
}
