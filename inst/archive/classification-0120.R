#' Create minimal dataset for one subject
#' Expects data for **exactly one subject** (id, time, egfr, egfr_add columns).
#' @param df Data frame with columns: id (or custom name), time, egfr (one subject only)
#' @param id_col Name of the ID column (default: "id")
#' @param time_col Name of time column (default: "time")
#' @param egfr_col Name of eGFR column (default: "egfr")
#' @param add_cols Names of additional columns (default: character(0))
minimal_1traj <- function(
    df,
    id_col   = "id",
    time_col = "time",
    egfr_col = "egfr",
    add_cols = character(0)
) {
  # ── Input validation ───────────────────────────────────────────────────────
  required_cols <- c(id_col, time_col, egfr_col)
  if (length(add_cols) > 0) {
    required_cols <- c(required_cols, add_cols)
  }
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Input df must contain columns: ", paste(required_cols, collapse = ", "),
      "\nMissing: ", paste(missing_cols, collapse = ", ")
    )
  }

  unique_ids <- unique(df[[id_col]])
  if (length(unique_ids) != 1L) {
    stop(
      "Input df must contain data for exactly one subject. ",
      "Found ", length(unique_ids), " unique IDs: ",
      paste(unique_ids, collapse = ", ")
    )
  }

  # ── Build minimal dataset ──────────────────────────────────────────────────
  select_cols <- c(
    id   = id_col,           # rename to standard "id"
    time = time_col,         # rename to standard "time"
    egfr = egfr_col          # rename to standard "egfr"
  )
  
  if (length(add_cols) > 0) {
    select_cols <- c(select_cols, add_cols)
  }

  minimal_df <- df %>%
    dplyr::select(all_of(select_cols)) %>%
    dplyr::filter(!is.na(time), !is.na(egfr))
  
  minimal_df
}

#' Classify a single eGFR trajectory
#'
#' Fits linear, quadratic, and segmented (max 1 knot) models to one subject's eGFR trajectory,
#' selects the best model using BIC (with tie-breaker favoring segmented),
#' and returns predictions + nested model summary.
#'
#' @param df Data frame with columns: id (or custom name), time, egfr (one subject only)
#' @param id_col Name of the ID column (default: "id")
#' @param time_col Name of time column (default: "time")
#' @param egfr_col Name of eGFR column (default: "egfr")
#' @param add_cols Names of additional columns to keep (default: character(0))
#' @param bic_tie_threshold BIC difference threshold for preferring segmented (default: 4)
#' @param ci_level Confidence level for prediction intervals (default: 0.95)
#'
#' @return Named list with:
#'   - pattern
#'   - breakpoint_time
#'   - nobs
#'   - data (augmented with predictions/CIs)
#'   - model_summary (one-row tibble with nested lists)
#' @export
classify_single_trajectory <- function(
    df,
    id_col             = "id",
    time_col           = "time",
    egfr_col           = "egfr",
    add_cols           = character(0),
    bic_tie_threshold  = 4,
    ci_level           = 0.95
) {
  minimal_df <- minimal_1traj(
    df         = df,
    id_col     = id_col,
    time_col   = time_col,
    egfr_col   = egfr_col,
    add_cols   = add_cols
  )
  
  n_valid <- nrow(minimal_df)
  n_miss  <- nrow(df) - n_valid
  subject_id <- unique(minimal_df[[id_col]])[1]

  # Early exit: insufficient data
  if (n_valid < 5L) {
    model_sum <- tibble::tibble(
      id              = subject_id,
      pattern         = "insufficient_data",
      info            = list(nobs = n_valid, nmiss = n_miss, breakpoint_time = NA_real_, bp_se = NA_real_),
      linear          = list(),
      quadratic       = list(),
      segmented       = list()
    )
    
    out <- list(
      pattern         = "insufficient_data",
      breakpoint_time = NA_real_,
      nobs            = n_valid,
      data            = minimal_df,
      model_summary   = model_sum
    )
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

  # ── Model list & ICs ───────────────────────────────────────────────────────
  models_list <- list(
    linear    = lm_lin,
    quadratic = if (!inherits(lm_quad, "try-error")) lm_quad else NULL,
    segmented = if (!inherits(seg_fit, "try-error") && !is.null(seg_fit$psi)) seg_fit else NULL
  )

  aic_vals <- purrr::map_dbl(models_list, ~ if (!is.null(.x)) AIC(.x) else NA_real_)
  bic_vals <- purrr::map_dbl(models_list, ~ if (!is.null(.x)) BIC(.x) else NA_real_)

  # ── Select best model ──────────────────────────────────────────────────────
  valid_models <- names(bic_vals)[!is.na(bic_vals)]
  if (length(valid_models) == 0) {
    best_model_name <- "linear"
  } else {
    min_bic <- min(bic_vals[valid_models], na.rm = TRUE)
    best_model_name <- valid_models[which.min(bic_vals[valid_models])]

    # Tie-breaker: prefer segmented if close enough
    if ("segmented" %in% valid_models &&
        bic_vals["segmented"] <= min_bic + bic_tie_threshold &&
        !is.null(models_list$segmented)) {
      best_model_name <- "segmented"
    }
  }

  # ── Breakpoint ─────────────────────────────────────────────────────────────
  bp_time <- NA_real_
  bp_se   <- NA_real_
  if (best_model_name == "segmented" && !is.null(models_list$segmented) &&
      !is.null(models_list$segmented$psi) && nrow(models_list$segmented$psi) > 0) {
    bp_time <- round(models_list$segmented$psi[1, "Est."], 1)
    bp_se   <- models_list$segmented$psi[1, "St.Err"]
  }

  # ── Predictions + intervals for all models ─────────────────────────────────
  pred_linear <- suppressWarnings(predict(lm_lin, interval = "prediction", level = ci_level))
  
  pred_quad <- if (!is.null(models_list$quadratic)) {
    suppressWarnings(predict(models_list$quadratic, interval = "prediction", level = ci_level))
  } else NULL
  
  pred_seg <- if (!is.null(models_list$segmented)) {
    suppressWarnings(predict(models_list$segmented, interval = "prediction", level = ci_level))
  } else NULL

  pred_df <- minimal_df %>%
    dplyr::mutate(
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

  # ── Nested model_summary using purrr::map ──────────────────────────────────
  extract_coefs <- function(model, model_name) {
    if (is.null(model)) return(list())
    
    s <- summary(model)
    cf <- coef(model)
    
    if (model_name == "linear") {
      list(
        intercept     = cf[1],
        slope         = cf[2],
        intercept_se  = s$coefficients[1, "Std. Error"],
        slope_se      = s$coefficients[2, "Std. Error"],
        AIC           = aic_vals[[model_name]],
        BIC           = bic_vals[[model_name]]
      )
    } else if (model_name == "quadratic") {
      list(
        intercept     = cf[1],
        linear        = cf[2],
        quadratic     = cf[3],
        intercept_se  = s$coefficients[1, "Std. Error"],
        linear_se     = s$coefficients[2, "Std. Error"],
        quadratic_se  = s$coefficients[3, "Std. Error"],
        AIC           = aic_vals[[model_name]],
        BIC           = bic_vals[[model_name]]
      )
    } else if (model_name == "segmented") {
      if (is.null(model$psi) || nrow(model$psi) == 0) return(list())
      list(
        intercept_pre  = cf[1],
        slope_pre      = cf[2],
        slope_post     = cf[3],
        intercept_se   = s$coefficients[1, "Std. Error"],
        slope_pre_se   = s$coefficients[2, "Std. Error"],
        slope_post_se  = s$coefficients[3, "Std. Error"],
        AIC            = aic_vals[[model_name]],
        BIC            = bic_vals[[model_name]]
      )
    } else {
      list()
    }
  }

  model_summary <- tibble::tibble(
  id        = subject_id,
  pattern   = best_model_name,
  
  info      = list(list(
    nobs            = n_valid,
    nmiss           = n_miss,
    breakpoint_time = bp_time,
    bp_se           = bp_se
  )),
  
  linear    = list(extract_coefs(models_list$linear,    "linear")),
  quadratic = list(extract_coefs(models_list$quadratic, "quadratic")),
  segmented = list(extract_coefs(models_list$segmented, "segmented"))
)
  # ── Final return ───────────────────────────────────────────────────────────
  list(
    pattern         = best_model_name,
    breakpoint_time = bp_time,
    nobs            = n_valid,
    data            = pred_df,
    model_summary   = model_summary
  )
}