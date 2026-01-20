#' Classify a single eGFR trajectory
#'
#' Expects data for **exactly one subject** (id, time, egfr columns).
#' Fits linear, quadratic, and segmented (max 1 knot) models, selects best
#' using BIC (with tie-breaker favoring segmented), and returns predictions
#' for all models + detailed coefficient summary.
#'
#' @param df Data frame with columns: id (or custom name), time, egfr (one subject only)
#' @param id_col Name of the ID column (default: "id")
#' @param time_col Name of time column (default: "time")
#' @param egfr_col Name of eGFR column (default: "egfr")
#' @param bic_tie_threshold BIC difference threshold for preferring segmented (default: 4)
#' @param ci_level Confidence level for prediction intervals (default: 0.95)
#'
#' @return Named list with class "egfrtraj_single":
#'   - pattern: "insufficient_data", "linear", "quadratic", "with_breakpoint"
#'   - selected_model_name: "linear", "quadratic", "segmented"
#'   - selected_model: winning model object
#'   - aic_values, bic_values: named vectors for all models
#'   - breakpoint_time: estimated breakpoint (NA otherwise)
#'   - nobs: number of valid observations used
#'   - nmiss: number of missing observations
#'   - data: input data + columns for predictions/CIs from all models
#'     (lin_fitted, lin_conf.low, lin_conf.high, etc.)
#'   - model_summary: one-row tibble with id, nobs, nmiss, breakpoint_time,
#'     all coefficients/SEs/AIC/BIC for the three models
#'   - all_models: list of fitted objects (or NULL)
#' @export
classify_single_trajectory <- function(
    df,
    id_col             = "id",
    time_col           = "time",
    egfr_col           = "egfr",
    bic_tie_threshold  = 4,
    ci_level           = 0.95
) {
  # ── Input validation ───────────────────────────────────────────────────────
  required_cols <- c(id_col, time_col, egfr_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Input df must contain columns: ", id_col, ", ", time_col, ", ", egfr_col,
         ". Missing: ", paste(missing_cols, collapse = ", "))
  }

  unique_ids <- unique(df[[id_col]])
  if (length(unique_ids) != 1L) {
    stop("Input df must contain data for exactly one subject. Found ", length(unique_ids),
         " unique IDs: ", paste(unique_ids, collapse = ", "))
  }

  subject_id <- unique_ids[1]

  # Minimal clean data
  minimal_df <- df %>%
    dplyr::select(
      id   = !!id_col,
      time = !!time_col,
      egfr = !!egfr_col
    ) %>%
    dplyr::filter(!is.na(time), !is.na(egfr))

  n_valid <- nrow(minimal_df)
  n_miss  <- nrow(df) - n_valid

  # Early exit: insufficient data
  if (n_valid < 5L) {
    model_sum <- tibble::tibble(
      id                = subject_id,
      nobs              = n_valid,
      nmiss             = n_miss,
      breakpoint_time   = NA_real_,
      bp_se             = NA_real_
    )
    out <- list(
      pattern             = "insufficient_data",
      selected_model_name = NA_character_,
      selected_model      = NULL,
      aic_values          = NA_real_,
      bic_values          = NA_real_,
      breakpoint_time     = NA_real_,
      nobs                = n_valid,
      data                = minimal_df,
      model_summary       = model_sum,
      all_models          = list(linear = NULL, quadratic = NULL, segmented = NULL)
    )
    class(out) <- c("egfrtraj_single", "list")
    return(out)
  }

  # ── Fit all three models ───────────────────────────────────────────────────
  lin_formula <- reformulate("time", "egfr")
  lm_lin <- lm(lin_formula, data = minimal_df)

  quad_formula <- reformulate(c("time", "I(time^2)"), "egfr")
  lm_quad <- try(lm(quad_formula, data = minimal_df), silent = TRUE)

  seg_fit <- try(
    segmented::segmented(
      lm_lin,
      seg.Z = ~ time,
      psi   = median(minimal_df$time),
      control = segmented::seg.control(K = 1L, n.boot = 10L, it.max = 30L)
    ),
    silent = TRUE
  )

  # ── AIC & BIC for all models ───────────────────────────────────────────────
  models_list <- list(
    linear    = lm_lin,
    quadratic = if (!inherits(lm_quad, "try-error")) lm_quad else NULL,
    segmented = if (!inherits(seg_fit, "try-error") && !is.null(seg_fit$psi)) seg_fit else NULL
  )

  aic_vals <- sapply(models_list, function(m) if (!is.null(m)) AIC(m) else NA_real_)
  bic_vals <- sapply(models_list, function(m) if (!is.null(m)) BIC(m) else NA_real_)

  # ── Select best model by BIC ──────────────────────────────────────────────
  valid_models <- names(bic_vals)[!is.na(bic_vals)]
  if (length(valid_models) == 0) {
    best_model_name <- "linear"
    best_model <- lm_lin
  } else {
    min_bic <- min(bic_vals[valid_models], na.rm = TRUE)
    best_model_name <- valid_models[which.min(bic_vals[valid_models])]

    # Tie-breaker: prefer segmented if within threshold
    if ("segmented" %in% valid_models &&
        bic_vals["segmented"] <= min_bic + bic_tie_threshold &&
        !is.null(models_list$segmented)) {
      best_model_name <- "segmented"
    }

    best_model <- models_list[[best_model_name]]
  }

  # Breakpoint if segmented selected
  bp_time <- NA_real_
  bp_se   <- NA_real_
  if (best_model_name == "segmented" && !is.null(best_model$psi)) {
    bp_time <- round(best_model$psi[1L, "Est."], 1)
    ci_obj <- tryCatch(
      confint(best_model, level = 0.95, method = "delta"),
      error = function(e) NULL
    )
    if (!is.null(ci_obj) && nrow(ci_obj) > 0) {
      bp_se <- ci_obj[1, "St. Err."]
    }
  }

  # ── In-sample predictions + CI for ALL models ─────────────────────────────
  pred_linear <- suppressWarnings(predict(lm_lin, interval = "prediction", level = ci_level))
  pred_quad   <- if (!is.null(models_list$quadratic)) suppressWarnings(predict(models_list$quadratic, interval = "prediction", level = ci_level)) else NULL
  pred_seg    <- if (!is.null(models_list$segmented)) suppressWarnings(predict(models_list$segmented, interval = "prediction", level = ci_level)) else NULL

  pred_df <- minimal_df %>%
    mutate(
      lin_fitted     = pred_linear[, "fit"],
      lin_conf.low   = pred_linear[, "lwr"],
      lin_conf.high  = pred_linear[, "upr"],
      quad_fitted    = if (!is.null(pred_quad)) pred_quad[, "fit"] else NA_real_,
      quad_conf.low  = if (!is.null(pred_quad)) pred_quad[, "lwr"] else NA_real_,
      quad_conf.high = if (!is.null(pred_quad)) pred_quad[, "upr"] else NA_real_,
      seg_fitted     = if (!is.null(pred_seg)) pred_seg[, "fit"] else NA_real_,
      seg_conf.low   = if (!is.null(pred_seg)) pred_seg[, "lwr"] else NA_real_,
      seg_conf.high  = if (!is.null(pred_seg)) pred_seg[, "upr"] else NA_real_
    )

  # ── Model summary tibble (one row with coefficients, SEs, AIC, BIC) ────────
  model_sum <- tibble::tibble(
    id                  = subject_id,
    nobs                = n_valid,
    nmiss               = n_miss,
    breakpoint_time     = bp_time,
    bp_se               = bp_se,
    lin_intercept       = coef(lm_lin)[1],
    lin_slope           = coef(lm_lin)[2],
    lin_intercept_se    = summary(lm_lin)$coefficients[1, "Std. Error"],
    lin_slope_se        = summary(lm_lin)$coefficients[2, "Std. Error"],
    lin_AIC             = aic_vals["linear"],
    lin_BIC             = bic_vals["linear"],
    quad_intercept      = if (!is.null(lm_quad)) coef(lm_quad)[1] else NA_real_,
    quad_linear         = if (!is.null(lm_quad)) coef(lm_quad)[2] else NA_real_,
    quad_quadratic      = if (!is.null(lm_quad)) coef(lm_quad)[3] else NA_real_,
    quad_intercept_se   = if (!is.null(lm_quad)) summary(lm_quad)$coefficients[1, "Std. Error"] else NA_real_,
    quad_linear_se      = if (!is.null(lm_quad)) summary(lm_quad)$coefficients[2, "Std. Error"] else NA_real_,
    quad_quadratic_se   = if (!is.null(lm_quad)) summary(lm_quad)$coefficients[3, "Std. Error"] else NA_real_,
    quad_AIC            = aic_vals["quadratic"],
    quad_BIC            = bic_vals["quadratic"],
    seg_pre_intercept   = if (!is.null(seg_fit) && nrow(seg_fit$psi) > 0) coef(seg_fit)[1] else NA_real_,
    seg_pre_slope       = if (!is.null(seg_fit) && nrow(seg_fit$psi) > 0) coef(seg_fit)[2] else NA_real_,
    seg_post_slope      = if (!is.null(seg_fit) && nrow(seg_fit$psi) > 0) coef(seg_fit)[3] else NA_real_,
    seg_pre_intercept_se = if (!is.null(seg_fit) && nrow(seg_fit$psi) > 0) summary(seg_fit)$coefficients[1, "Std. Error"] else NA_real_,
    seg_pre_slope_se    = if (!is.null(seg_fit) && nrow(seg_fit$psi) > 0) summary(seg_fit)$coefficients[2, "Std. Error"] else NA_real_,
    seg_post_slope_se   = if (!is.null(seg_fit) && nrow(seg_fit$psi) > 0) summary(seg_fit)$coefficients[3, "Std. Error"] else NA_real_,
    seg_AIC             = aic_vals["segmented"],
    seg_BIC             = bic_vals["segmented"]
  )

  # ── Final return ──────────────────────────────────────────────────────────
  out <- list(
    pattern             = best_model_name,
    selected_model      = models_list[[best_model_name]],
    breakpoint_time     = bp_time,
    nobs                = n_valid,
    data                = pred_df,
    model_summary       = model_sum,
    all_models          = models_list
  )

  class(out) <- c("egfrtraj_single", "list")
  out
}
