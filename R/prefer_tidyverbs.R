#' Prefer common tidyverse verbs (dplyr + tidyr)
#'
#' Gives priority to dplyr and tidyr verbs, including select/rename/distinct.
#'
#' @export
#' @examples
#' \dontrun{
#' library(conflicted)
#' library(dplyr)
#' library(tidyr)
#' library(segmented)
#' library(egfrtraj)
#' data(example_egfr_data)
#' egfrtraj::prefer_tidyverbs()
#'
#' # Now filter(), lag(), select() etc. will use dplyr versions
#' example_egfr_data |> filter(time > 0)
#' }
prefer_tidyverbs <- function() {
  if (!requireNamespace("conflicted", quietly = TRUE)) {
    stop("Package 'conflicted' is required.", call. = FALSE)
  }

  conflicted::conflicts_prefer(
    dplyr::filter, dplyr::mutate, dplyr::summarise, dplyr::summarize,
    dplyr::group_by, dplyr::ungroup, dplyr::arrange, dplyr::count,
    dplyr::add_count, dplyr::slice, dplyr::slice_head, dplyr::slice_tail,
    dplyr::slice_min, dplyr::slice_max, dplyr::slice_sample,
    dplyr::lag, dplyr::lead, dplyr::across, dplyr::if_else,
    dplyr::case_when, dplyr::between, dplyr::first, dplyr::last,
    dplyr::n, dplyr::n_distinct,

    dplyr::select, dplyr::rename, dplyr::distinct, dplyr::transmute,

    tidyr::pivot_longer, tidyr::pivot_wider, tidyr::drop_na,
    tidyr::replace_na, tidyr::fill, tidyr::complete, tidyr::expand,
    tidyr::nest, tidyr::unnest, tidyr::separate, tidyr::unite,
    tidyr::hoist,

    dplyr::everything(),

    .quiet = TRUE
  )
  invisible(NULL)
}


