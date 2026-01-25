#' Check structure of longitudinal data in long format
#'
#' @description
#' Verifies whether a longitudinal dataset in long format has:
#' - one row per time point per subject
#' - rows properly grouped and sorted by subject ID and time
#' - no duplicate (id, time) combinations
#' - classifies remaining variables as time-varying or time-invariant
#'
#' @details
#' The function:
#' 1. Converts the ID column to character (if it is a factor). Message is issued
#' 2. Sorts the data by ID and time
#' 3. Checks whether rows with the same ID are consecutive and times are non-decreasing within each ID
#' 4. Detects duplicate (ID, time) pairs
#' 5. Determines which variables vary within subjects (time-varying) vs. are constant per subject (time-invariant)
#'
#' @note
#' - The function assumes `time_col` is numeric (or something `diff()` can handle). If `time` is a Date/POSIXct,
#'   consider converting it to numeric days beforehand.
#' - Because the package depends on `segmented` (which loads `MASS`), we deliberately avoid broad preferences
#'   for `dplyr::select`. Use explicit `dplyr::select()` when needed in your code.
#'
#' @param data A data frame in long format (one row per observation per subject)
#' @param id_col String. Name of the subject identifier column (default: `"id"`)
#' @param time_col String. Name of the time/visit/measurement occasion column (default: `"time"`)
#'
#' @return A named list with three components:
#'   \describe{
#'     \item{`is_properly_sorted`}{Logical. Are rows properly grouped by ID and sorted by time within groups?}
#'     \item{`has_no_duplicate_id_time`}{Logical. Are there no duplicate (id, time) combinations?}
#'     \item{`classification`}{A tibble with columns `variable`, `is_time_varying`, `status` showing
#'       whether each variable (excluding `id_col` and `time_col`) is time-varying or time-invariant.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' data(example_egfr_data)
#' result <- check_long_structure(example_egfr_data)
#' str(result)
#' result$classification
#'
#' }
#'
#' @importFrom dplyr arrange group_by summarise pull n_distinct any_of across
#'   mutate if_else count filter
#" @importFrom tidyr pivot_longer 
#' @importFrom rlang .data
#' @importFrom stats na.omit
#'
#' @export
check_long_structure <- function(data, id_col = "id", time_col = "time") {

  # ── Input validation ──────────────────────────────────────────────────────────
  required_cols <- c(id_col, time_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Required column(s) not found: ", paste(missing_cols, collapse = ", "))
  }

  # ── Convert id to character (handles factors safely) ──────────────────────────
  if (is.factor(data[[id_col]])) {
    message("Converting '", id_col, "' from factor to character")
    data[[id_col]] <- as.character(data[[id_col]])
  }

  # Sort once (safer and needed for consecutive check)
  data <- data |>
    dplyr::arrange(.data[[id_col]], .data[[time_col]])

  # ── 1. Check proper grouping & sorting ───────────────────────────────────────

  # Are same-ID rows consecutive? (no interleaving of subjects)
  consecutive_groups <- {
    id_vec <- data[[id_col]]
    same_as_prev <- id_vec == dplyr::lag(id_vec, default = id_vec[1])
    same_as_next <- id_vec == dplyr::lead(id_vec, default = id_vec[nrow(data)])
    all(same_as_prev | same_as_next)
  }

  # Within each id: time non-decreasing
  time_sorted_within <- data |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(
      time_diff_ok = all(diff(.data[[time_col]]) >= 0 | is.na(diff(.data[[time_col]]))),
      .groups = "drop"
    ) |>
    dplyr::pull(time_diff_ok) |>
    all(na.rm = TRUE)

  is_properly_sorted <- consecutive_groups && time_sorted_within

  # ── 2. Check for duplicate (id, time) pairs ──────────────────────────────────
  dups <- data |>
    dplyr::count(.data[[id_col]], .data[[time_col]], name = "cnt") |>
    dplyr::filter(cnt > 1L)

  has_duplicates <- nrow(dups) > 0

  # ── 3. Classify time-varying vs time-invariant variables ─────────────────────
  exclude_cols <- c(id_col, time_col)

  varying <- data |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(
      dplyr::across(
        .cols = -dplyr::any_of(exclude_cols),
        .fns  = ~ n_distinct(.) > 1L | any(is.na(.)) != all(is.na(.))
      ),
      .groups = "drop"
    ) |>
    dplyr::summarise(dplyr::across(-dplyr::any_of(id_col), any)) |>
    tidyr::pivot_longer(
      cols      = dplyr::everything(),
      names_to  = "variable",
      values_to = "is_time_varying"
    ) |>
    dplyr::mutate(status = dplyr::if_else(
      is_time_varying,
      "time-varying",
      "time-invariant"
    ))

  # ── Return ───────────────────────────────────────────────────────────────────
  list(
    is_properly_sorted       = is_properly_sorted,
    has_no_duplicate_id_time = !has_duplicates,
    classification           = varying
  )
}