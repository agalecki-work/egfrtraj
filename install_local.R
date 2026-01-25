
# Make sure you have devtools (once)
# install.packages("devtools")

# Then, from any working directory:
# devtools::install("C:/ATG/github/egfrtraj")
# or using Windows-style backslashes (both usually work)
devtools::install("C:\\ATG\\github\\egfrtraj") # skips tests by default

library(conflicted)
library(dplyr)
library(tidyr)     # ← this brings pivot_longer() and pivot_wider()
library(egfrtraj)
library(segmented)
conflicts_prefer(
  dplyr::filter,
  dplyr::lag,
  dplyr::select,
  dplyr::first, dplyr::last, dplyr::between,   # common masked ones
  dplyr::everything()                          # catch the rest
)
conflicts_prefer(dplyr::everything())

getNamespaceExports("egfrtraj")
data(example_egfr_data)
colnames(example_egfr_data)

check_long_structure <- function(data, id_col = "id", time_col = "time") {
  
  required_cols <- c(id_col, time_col)
  if (!all(required_cols %in% names(data))) {
    stop("Required columns not found: ", paste(setdiff(required_cols, names(data)), collapse = ", "))
  }
  
  if (is.factor(data[[id_col]])) {
    message("Converting '", id_col, "' from factor to character")
  }
  
  data <- data |>
    mutate("{id_col}" := as.character(.data[[id_col]])) |>
    arrange(.data[[id_col]], .data[[time_col]])
  
  # 1. Check proper grouping & sorting
  consecutive_groups <- all(
    data[[id_col]] == dplyr::lag(data[[id_col]], default = data[[id_col]][1]) |
    data[[id_col]] == dplyr::lead(data[[id_col]], default = data[[id_col]][nrow(data)])
  )
  
  time_sorted_within <- data |>
    group_by(.data[[id_col]]) |>
    summarise(
      time_diff_ok = all(diff(.data[[time_col]]) >= 0 | is.na(diff(.data[[time_col]]))),
      .groups = "drop"
    ) |>
    pull(time_diff_ok) |>
    all(na.rm = TRUE)
  
  is_properly_sorted <- consecutive_groups && time_sorted_within
  
  # 2. Duplicates – use explicit dplyr:: to be safe against conflicts
  dups <- data |>
    dplyr::count(.data[[id_col]], .data[[time_col]], name = "cnt") |>
    dplyr::filter(cnt > 1)
  
  has_duplicates <- nrow(dups) > 0
  
  # 3. Classification of variables
  exclude_cols <- c(id_col, time_col)
  
  varying <- data |>
    group_by(.data[[id_col]]) |>
    summarise(
      across(
        .cols = -any_of(exclude_cols),
        .fns  = ~ n_distinct(.) > 1L | any(is.na(.)) != all(is.na(.))
      ),
      .groups = "drop"
    ) |>
    summarise(across(-any_of(id_col), any)) |>
    pivot_longer(
      cols      = everything(),
      names_to  = "variable",
      values_to = "is_time_varying"
    ) |>
    mutate(status = if_else(is_time_varying, "time-varying", "time-invariant"))
  
  list(
    is_properly_sorted       = is_properly_sorted,
    has_no_duplicate_id_time = !has_duplicates,
    classification           = varying
  )
}


# Usage
result <- check_long_structure(example_egfr_data, id_col = "id", time_col = "time")
result

example_subjects <- example_egfr_data |>
  group_by(id) |>
  slice_min(time, with_ties = FALSE) |>   # gets earliest time explicitly
  ungroup() |>
  select(id, true_pattern)

table(example_subjects$true_pattern)

df_single <- example_egfr_data %>%
  filter(id == "11")

# Classify all subjects
traj_all <- classify_multiple_trajectories(example_egfr_data)

# Inspect structure
names(traj_all)               # "models" "fitted"
class(traj_all)               # "egfr_traj_tibbles"

# Model summary (one row per subject)
traj_all$models

# Fitted values (long format, one row per observation)
head(traj_all$fitted)


