#' Subset an egfrtibble object by subject IDs
#'
#' @description
#' Creates a subset of an `"egfrtibble"` object containing only the specified
#' subject IDs. If `ids = NULL`, selects only the first subject present in the
#' object. Invalid/missing IDs are silently skipped. The `final_subjects`
#' attribute is updated accordingly.
#'
#' If the resulting subset is empty, returns an empty tibble of class
#' `"egfrtibble"` with updated metadata (and a warning).
#'
#' @param x An object of class `"egfrtibble"`
#' @param ids Character vector of subject IDs to keep.  
#'   If `NULL` (default), only the first subject in the object is selected.
#'
#' @return An `"egfrtibble"` object (possibly empty) with updated
#'   `origdata_info$final_subjects` attribute:
#'   - `"one"`    if exactly 1 subject remains
#'   - `"subset"` if 1 < n_subjects < original n_subjects_final
#'   - `"All"`    if all original subjects are kept (rare)
#'   - `"empty"`  if no subjects remain
#'
#' @examples
#' \dontrun{
#' # First create a full egfrtibble object
#' res <- extract_egfrtibble(example_egfr_data)
#' 
#' # Get all unique IDs from the original data
#' all_ids <- unique(example_egfr_data$id)
#' 
#' 
#' # ── 1. Subset to subjects with specific true_pattern (using original data) ───
#' # Note: true_pattern is not kept in res → we use the original data frame
#' selected_ids <- example_egfr_data |>
#'   dplyr::filter(true_pattern %in% c("linear", "with_breakpoint")) |>
#'   dplyr::distinct(id) |>
#'   dplyr::pull(id)
#' 
#' print(selected_ids)   # show which IDs were selected
#' 
#' sub <- subset_egfrtibble(res, ids = selected_ids)
#' 
#' nrow(sub)                              # number of rows in subset
#' n_distinct(sub$id)                     # number of subjects kept
#' attr(sub, "origdata_info")$final_subjects   # usually "subset"
#' 
#' 
#' # ── 2. Select only the first subject (ids = NULL) ────────────────────────────
#' sub_first <- subset_egfrtibble(res)
#' 
#' n_distinct(sub_first$id)               # should be 1
#' attr(sub_first, "origdata_info")$final_subjects   # "one"
#' sub_first$id[1]                        # which subject was chosen
#' 
#' 
#' # ── 3. Request IDs that don't exist (silently skipped) ───────────────────────
#' bad_ids <- c("nonexistent1", "nonexistent2", selected_ids[3])
#' sub_bad <- subset_egfrtibble(res, ids = bad_ids)
#' 
#' n_distinct(sub_bad$id)                 # only 1 subject kept
#' attr(sub_bad, "origdata_info")$final_subjects   # "one"
#' 
#' 
#' # ── 4. Empty subset (no matching IDs) ────────────────────────────────────────
#' sub_empty <- subset_egfrtibble(res, ids = c("9999", "ZZZ"))
#' 
#' nrow(sub_empty)                        # 0 rows
#' n_distinct(sub_empty$id)               # 0 subjects
#' attr(sub_empty, "origdata_info")$final_subjects   # "empty"
#' # Note: warning is shown when result is empty
#' 
#' 
#' # ── 5. Chain with other operations ───────────────────────────────────────────
#' # Keep only one subject and prepare for single-case analysis
#' one_sub <- subset_egfrtibble(res, ids = selected_ids[1])
#' 
#' 
#' # ── Invalid usage examples (will stop with error) ────────────────────────────
#' # \dontrun{
#' #   subset_egfrtibble(example_egfr_data)          # not an egfrtibble → error
#' #   subset_egfrtibble(res, ids = 123)             # non-character ids → error
#' # }
#' }
#' @export


subset_egfrtibble <- function(x, ids = NULL) {

  # ── Input validation ──────────────────────────────────────────────────────────
  stopifnot(
    inherits(x, "egfrtibble"),
    is.null(ids) || is.character(ids)
  )

  info <- attr(x, "origdata_info")

  if (is.null(info) || !is.list(info)) {
    stop("Input object is missing 'origdata_info' attribute")
  }

  original_n_subjects <- info$n_subjects_final
  original_flag       <- info$final_subjects

  # ── Handle ids = NULL → select first subject ────────────────────────────────
  if (is.null(ids)) {
    if (nrow(x) == 0) {
      warning("Input egfrtibble is already empty → returning empty subset")
      new_flag <- "empty"
    } else {
      first_id <- x$id[1]
      ids <- first_id
      message("No ids provided → selecting first subject: ", first_id)
    }
  }

  # ── Subset the data ──────────────────────────────────────────────────────────
  if (length(ids) == 0) {
    # empty request → return empty
    subset_data <- x[0, , drop = FALSE]
  } else {
    # keep only rows matching requested ids (skip non-existing ones)
    subset_data <- x |>
      dplyr::filter(.data$id %in% ids)
  }

  # ── Determine new flag ───────────────────────────────────────────────────────
  n_after <- dplyr::n_distinct(subset_data$id)

  if (n_after == 0) {
    new_flag <- "empty"
    warning("No matching subjects found → returning empty egfrtibble")
  } else if (n_after == 1) {
    new_flag <- "one"
  } else if (n_after < original_n_subjects) {
    new_flag <- "subset"
  } else {
    new_flag <- "All"   # rare: all original subjects were requested and kept
  }

  # ── Update metadata ──────────────────────────────────────────────────────────
  new_info <- info
  new_info$n_rows_after_filter   <- nrow(subset_data)
  new_info$n_subjects_final      <- n_after
  new_info$final_subjects        <- new_flag
  new_info$subset_ids_requested  <- if (is.null(ids)) "first subject" else ids
  new_info$subset_ids_kept       <- unique(subset_data$id)

  # ── Return updated object ────────────────────────────────────────────────────
  class(subset_data) <- c("egfrtibble", class(subset_data))
  attr(subset_data, "origdata_info") <- new_info

  subset_data
}

#   # Example future method dispatch based on final_subjects
# if (attr(one_sub, "origdata_info")$final_subjects == "one") {
#   message("Single-subject trajectory selected – special single-case handling")
#}
 
