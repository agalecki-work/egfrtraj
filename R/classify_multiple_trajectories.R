# classify_multiple_trajectories.R

# =============================================================================
# classify_multiple_trajectories
# =============================================================================

#' Classify multiple eGFR trajectories
#'
#' Applies `classify_single_trajectory()` to each subject in the dataset and returns
#' a simplified, tibble-based structure with two components: one row per subject for model
#' summaries, and one row per observation for fitted values.
#'
#' @param data A data frame containing longitudinal eGFR data for one or more subjects.
#' @param id_col Name of the subject ID column (default: `"id"`).
#' @param time_col Name of the time column (default: `"time"`).
#' @param egfr_col Name of the eGFR column (default: `"egfr"`).
#' @param bic_tie_threshold BIC difference threshold for preferring segmented (default: `4`).
#' @param ci_level Confidence level for prediction intervals (default: `0.95`).
#' @param best_model_only Logical. If `TRUE` (default), keep only the selected model's
#'   fitted values and coefficients.
#' @param estimates_only Logical. If `TRUE`, drop standard error columns from `models`.
#'   Default `FALSE`.
#'
#' @return A list of class `"egfr_traj_tibbles"` with two tibbles:
#'   \describe{
#'     \item{models}{One row per subject with columns:
#'       \itemize{
#'         \item `id`: subject identifier
#'         \item `info`: nested list with `pattern`, `nobs`, `nmiss`, `breakpoint_time`, `bp_se`
#'         \item `linear`, `quadratic`, `segmented`: nested lists with model coefficients,
#'           SEs, AIC/BIC (simplified if `best_model_only = TRUE`)
#'       }
#'     }
#'     \item{fitted}{One row per observation with columns:
#'       \itemize{
#'         \item `id`: subject identifier
#'         \item `time`: time point
#'         \item `egfr`: observed eGFR value
#'         \item `fitted_value`: fitted value from the selected model
#'         \item `conf_low`, `conf_high`: prediction interval bounds
#'       }
#'     }
#'   }
#'
#' @details
#' - Subjects with fewer than 5 valid observations are handled as `"insufficient_data"`.
#' - Segmented models with unreliable breakpoint are downgraded to linear in the single
#'   classification step.
#' - Use `summary()` or `plot()` on individual elements when `best_model_only = FALSE`.
#'
#' @examples
#' library(egfrtraj)
#'
#' # Classify all subjects
#' traj_all <- classify_multiple_trajectories(example_egfr_data)
#'
#' # Inspect structure
#' names(traj_all)               # "models" "fitted"
#' class(traj_all)               # "egfr_traj_tibbles"
#'
#' # Model summary (one row per subject)
#' traj_all$models
#' 
#' # Fitted values (long format, one row per observation)
#' head(traj_all$fitted)
#'
#' # Example: plot one subject from fitted data
#' single_id <- "11"
#' single_traj <- traj_all$fitted %>% dplyr::filter(id == single_id)
#' # (use plot.egfr_traj on single trajectory if needed)
#'
#' @export
classify_multiple_trajectories <- function(
    data,
    id_col             = "id",
    time_col           = "time",
    egfr_col           = "egfr",
    bic_tie_threshold  = 4,
    ci_level           = 0.95,
    best_model_only    = TRUE,
    estimates_only     = FALSE
) {
  subjects <- unique(data[[id_col]])
  if (length(subjects) == 0) stop("No subjects found")

  # Temporary storage of single results
  single_results <- setNames(vector("list", length(subjects)), as.character(subjects))

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

    single_results[[as.character(subj)]] <- single
  }

  # ── Extract models tibble (one row per subject) ─────────────────────────────
  models_tibble <- purrr::map_dfr(single_results, function(single) {
    tibble::tibble(
      id     = single$id,
      info   = single$models$info,
      linear = single$models$linear,
      quadratic = single$models$quadratic,
      segmented = single$models$segmented
    )
  })

  # Simplify if requested
  if (best_model_only) {
    # For each row, keep only info + selected model columns
    models_tibble <- models_tibble |>
      dplyr::rowwise() |>
      dplyr::mutate(
        selected_prefix = switch(
          info[[1]]$pattern,
          linear    = "lin",
          quadratic = "quad",
          segmented = "seg",
          "lin"
        ),
        across(
          c(linear, quadratic, segmented),
          ~ if (cur_column_name() == paste0(selected_prefix, "_") || cur_column_name() == "info") . else list()
        )
      ) |>
      dplyr::ungroup()
  }

  if (estimates_only) {
    models_tibble <- models_tibble |>
      dplyr::select(-dplyr::ends_with("_se"))
  }

  # ── Extract fitted tibble (long format) ─────────────────────────────────────
  fitted_tibble <- purrr::map_dfr(single_results, function(single) {
    if (best_model_only) {
      # Simplified: only selected model
      single$trajectory |>
        dplyr::select(id, time, egfr, fitted_value = fitted, conf_low, conf_high)
    } else {
      # Full version: keep all fitted columns
      single$trajectory
    }
  })

  # ── Return list with class ─────────────────────────────────────────────────
  result <- list(
    models = models_tibble,
    fitted = fitted_tibble
  )

  class(result) <- "egfr_traj_tibbles"   # or "egfr_traj_collection" if preferred

  result
}

