# classify_single_trajectory.R

#' Classify a single eGFR trajectory
#'
#' Fits linear, quadratic, and segmented regression models to eGFR data from one subject.
#' Selects the best model using BIC with a tie-breaker favoring segmented if the difference
#' is within `bic_tie_threshold`. If a segmented model is selected but no reliable breakpoint
#' is estimated (`breakpoint_time` is NA), it is downgraded to linear.
#'
#' @param df A data frame containing data for **exactly one subject** with columns for
#'   time and eGFR (and optionally others).
#' @param id_col Name of the subject ID column (default: `"id"`).
#' @param time_col Name of the time column (default: `"time"`).
#' @param egfr_col Name of the eGFR column (default: `"egfr"`).
#' @param bic_tie_threshold BIC difference threshold for preferring segmented over simpler models
#'   (default: `4`).
#' @param ci_level Confidence level for prediction intervals (default: `0.95`).
#'
#' @return A list object of class `"egfr_traj"` with components:
#'   \describe{
#'     \item{id}{Subject identifier (factor or character).}
#'     \item{trajectory}{Tibble with observed data (`id`, `time`, `egfr`) and fitted values /
#'       prediction intervals for all models (`lin_fitted`, `lin_conf.low`, etc.).}
#'     \item{models}{One-row tibble with list-columns:
#'       \itemize{
#'         \item `info`: list with `pattern`, `nobs`, `nmiss`, `breakpoint_time`, `bp_se`
#'         \item `linear`, `quadratic`, `segmented`: lists with model coefficients, SEs, AIC/BIC
#'       }
#'     }
#'   }
#'
#' @details
#' - The function filters out rows with missing `time` or `egfr`.
#' - If fewer than 5 valid observations remain, returns a special `"insufficient_data"` result.
#' - Model selection prefers segmented when BIC improvement is small (tie-breaker).
#' - Segmented models with unreliable breakpoint (`NA`) are downgraded to linear.
#'
#' @examples
#' library(egfrtraj)
#'
#' # Select one subject (true pattern: linear)
#' single_df <- example_egfr_data %>% dplyr::filter(id == "1")
#'
#' # Classify
#' traj <- classify_single_trajectory(single_df)
#'
#' # Inspect results
#' class(traj)
#' methods(class = "egfr_traj")
#' summary(traj)
#' plot(traj)
#'
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

# Select best model by BIC (with tie-breaker favoring segmented)
valid_models <- names(bic_vals)[!is.na(bic_vals)]
best_pattern <- if (length(valid_models) == 0) "linear" else {
  min_bic <- min(bic_vals[valid_models])
  winner <- valid_models[which.min(bic_vals[valid_models])]

  if ("segmented" %in% valid_models && 
      bic_vals["segmented"] <= min_bic + bic_tie_threshold) {
    "segmented"
  } else winner
}

# Extract breakpoint (only for segmented)
bp_time <- NA_real_
bp_se   <- NA_real_
if (best_pattern == "segmented" && !is.null(models_list$segmented) && 
    !is.null(models_list$segmented$psi) && nrow(models_list$segmented$psi) > 0) {
  bp_time <- round(models_list$segmented$psi[1, "Est."], 1)
  bp_se   <- models_list$segmented$psi[1, "St.Err"]
}

# Downgrade segmented to linear if breakpoint is not reliable (NA)
if (best_pattern == "segmented" && is.na(bp_time)) {
  best_pattern <- "linear"
 # Flag it for downstream use (summary/plot)
  nested_info[[1]]$downgraded <- TRUE
  nested_info[[1]]$downgrade_reason <- "No reliable breakpoint estimated"
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


