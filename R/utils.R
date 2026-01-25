
#' Inspect egfr_traj object in detail (internal helper)
#'
#' Prints object structure, ID, and unnested contents of all list-columns
#' in the models tibble (info, linear, quadratic, segmented).
#'
#' @param o An object of class `"egfr_traj"`
#' @param verbose Logical. If `TRUE`, also prints the full trajectory tibble head/tail.
#'   Default: `FALSE`.
#'
#' @return Invisibly returns `o` (for piping)
#' @keywords internal
#' @noRd
inspect_egfr_traj <- function(o, verbose = FALSE) {
  if (!inherits(o, "egfr_traj")) {
    stop("Object must be of class 'egfr_traj'")
  }

  cat("\n", crayon::bold("Object structure:"), "\n")
  str(o, max.level = 1)

  cat("\n", crayon::bold("ID:"), "\n")
  print(o$id)

  cat("\n", crayon::bold("Trajectory columns:"), "\n")
  print(names(o$trajectory))

  if (verbose) {
    cat("\n", crayon::bold("Trajectory head (first 6 rows):"), "\n")
    print(head(o$trajectory))
    cat("\n", crayon::bold("Trajectory tail (last 6 rows):"), "\n")
    print(tail(o$trajectory))
  }

  summ <- o$models

  cat("\n", crayon::bold("Models tibble columns:"), "\n")
  print(names(summ))

  cat("\n", crayon::bold("Unnested $info:"), "\n")
  print(as_tibble(summ$info[[1]]))

  cat("\n", crayon::bold("Unnested $linear:"), "\n")
  print(as_tibble(summ$linear[[1]]))

  cat("\n", crayon::bold("Unnested $quadratic:"), "\n")
  print(as_tibble(summ$quadratic[[1]]))

  cat("\n", crayon::bold("Unnested $segmented:"), "\n")
  print(as_tibble(summ$segmented[[1]]))

  cat("\n", crayon::silver(strrep("─", 60)), "\n", sep = "")

  invisible(o)
}



#' Check if segmented package is available and loaded
#'
#' @return Logical
#' @keywords internal
has_segmented <- function() {
  requireNamespace("segmented", quietly = TRUE)
}

