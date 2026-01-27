#' Create standardized eGFR trajectory object
#'
#' @description
#' Prepares a clean tibble with standardized column names (`id`, `time`, `egfr`)
#' from a longitudinal dataset. First removes rows with missing values in the key
#' columns, then applies user-specified or default filtering.
#'
#' The returned object has class `"egfrtraj"` plus one of:
#' - `"egfrtraj_1"`   (single subject / trajectory)
#' - `"egfrtraj_m"`   (multiple subjects / trajectories)
#'
#' Metadata is attached as an attribute (`origdata_info`).
#'
#' @param data A data frame or tibble containing longitudinal eGFR data
#' @param key_vars A named list specifying the column names for the key variables.
#'   Default: `list(id_col = "id", time_col = "time", egfr_col = "egfr")`
#' @param filter A function that takes a data frame as first argument and returns
#'   a logical vector (rows to keep).  
#'   If `NULL` (default), keeps subjects with at least 5 non-missing `egfr` values.
#' @param ... Additional arguments passed to the `filter` function
#'
#' @return An object of class `"egfrtraj"` + `"egfrtraj_1"` or `"egfrtraj_m"`,
#'   inheriting from tibble. Columns are `id`, `time`, `egfr`.
#'   Metadata is available via `attr(x, "origdata_info")`.
#'
#' @examples
#' \dontrun{
#' # ── Basic usage (multiple subjects) ──────────────────────────────────────────
#' res_m <- create_egfrtraj(example_egfr_data)
#' 
#' class(res_m)                               # "egfrtraj_m" "egfrtraj" ...
#' attr(res_m, "origdata_info")$final_subjects   # "All"
#' 
#' 
#' # ── Single-subject data → different class ────────────────────────────────
#' one_id <- unique(example_egfr_data$id)[1]
#' one_subject <- example_egfr_data |>
#'   dplyr::filter(id == one_id)
#' 
#' res_1 <- create_egfrtraj(one_subject)
#' class(res_1)                               # "egfrtraj_1" "egfrtraj" ...
#' attr(res_1, "origdata_info")$final_subjects   # "one"
#' 
#' 
#' # ── Explicit key_vars (when names differ from defaults) ──────────────────────
#' res_explicit <- create_egfrtraj(
#'   data = example_egfr_data,
#'   key_vars = list(
#'     id_col   = "id",
#'     time_col = "time",
#'     egfr_col = "egfr"
#'   )
#' )
#' identical(res_explicit, res_m)             # TRUE when names match
#' 
#' 
#' # ── Custom filter example ────────────────────────────────────────────────────
#' res_custom <- create_egfrtraj(
#'   example_egfr_data,
#'   filter = function(d) {
#'     d |>
#'       dplyr::group_by(id) |>
#'       dplyr::summarise(max_e = max(egfr, na.rm = TRUE)) |>
#'       dplyr::mutate(keep = max_e >= 100) |>
#'       dplyr::right_join(d, by = "id") |>
#'       dplyr::pull(keep)
#'   }
#' )
#' attr(res_custom, "origdata_info")$final_subjects   # "All" or "one"
#' 
#' 
#' # ── Invalid input (caught early) ─────────────────────────────────────────────
#' \dontrun{
#'   # Duplicate source column → error
#'   create_egfrtraj(
#'     example_egfr_data,
#'     key_vars = list(id_col = "id", time_col = "id", egfr_col = "egfr")
#'   )
#'   # → Error: Duplicate column names detected in key_vars: id
#' }
#' }
#' @export
create_egfrtraj <- function(
    data,
    key_vars = list(id_col = "id", time_col = "time", egfr_col = "egfr"),
    filter = NULL,
    ...
) {
  # ── 1. Input validation ───────────────────────────────────────────────────────
  stopifnot(
    is.data.frame(data),
    is.list(key_vars),
    all(c("id_col", "time_col", "egfr_col") %in% names(key_vars)),
    all(c(key_vars$id_col, key_vars$time_col, key_vars$egfr_col) %in% colnames(data))
  )

  id_col   <- key_vars$id_col
  time_col <- key_vars$time_col
  egfr_col <- key_vars$egfr_col

  # Check for duplicate source columns
  src_cols <- c(id_col, time_col, egfr_col)
  if (anyDuplicated(src_cols)) {
    dup_cols <- src_cols[duplicated(src_cols)]
    stop(
      "Duplicate column names detected in key_vars: ",
      paste(dup_cols, collapse = ", "),
      ". The id_col, time_col, and egfr_col must refer to distinct columns."
    )
  }

  # ── 2. Select key columns + drop rows with any missing key variable ───────────
  valid_data <- data |>
    dplyr::select(dplyr::all_of(src_cols)) |>
    dplyr::filter(
      !is.na(.data[[id_col]]),
      !is.na(.data[[time_col]]),
      !is.na(.data[[egfr_col]])
    )

  # ── 3. Default filter: at least 5 valid (non-missing) egfr values per subject ─
  default_filter <- function(d) {
    d |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::summarise(
        n_valid = sum(!is.na(.data[[egfr_col]])),
        .groups = "drop"
      ) |>
      dplyr::mutate(keep = n_valid >= 5L) |>
      dplyr::right_join(d, by = id_col) |>
      dplyr::pull(keep)
  }

  # Apply filter (user or default)
  if (is.null(filter)) {
    keep_rows <- default_filter(valid_data)
    filter_desc <- "≥ 5 non-missing egfr per subject"
  } else {
    keep_rows <- filter(valid_data, ...)
    filter_desc <- paste0("user-defined: ", deparse(substitute(filter)))
  }

  filtered_data <- valid_data[keep_rows, , drop = FALSE]

  # ── 4. Standardize column names – safe and modern ─────────────────────────────
  out <- filtered_data |>
    dplyr::rename_with(
      ~ c("id", "time", "egfr"),
      .cols = dplyr::all_of(c(id_col, time_col, egfr_col))
    )

    # ── 5. Determine final_subjects flag & class suffix ───────────────────────────
  n_final <- dplyr::n_distinct(out$id)
  final_flag <- if (n_final == 0L) "empty" else if (n_final == 1L) "one" else "All"

  # Class prefix
  if (final_flag == "empty") {
    class_suffix <- "egfrtraj_empty"
  } else if (final_flag == "one") {
    class_suffix <- "egfrtraj_1"
  } else {
    class_suffix <- "egfrtraj_m"
  }

  # ── 6. Collect original data info ────────────────────────────────────────────
  origdata_info <- list(
    key_vars = list(
      id_col   = id_col,
      time_col = time_col,
      egfr_col = egfr_col
    ),
    n_rows_input      = nrow(data),
    n_subjects_input  = dplyr::n_distinct(data[[id_col]], na.rm = TRUE),
    n_rows_after_key_na_remove = nrow(valid_data),
    n_subjects_after_key_na_remove = dplyr::n_distinct(valid_data[[id_col]]),
    filter_description = filter_desc,
    n_rows_after_filter = nrow(out),
    n_subjects_final  = n_final,
    final_subjects    = final_flag,
    extra_colnames_input = setdiff(colnames(data), src_cols)
  )

  # ── 7. Set class and attribute ───────────────────────────────────────────────
  class(out) <- c(class_suffix, "egfrtraj", class(out))
  attr(out, "origdata_info") <- origdata_info

  out
}
