test_that("prepare_egfr_data standardizes columns and filters correctly", {
  # Mock data with different column names
  df <- tibble::tribble(
    ~participant, ~visit_year, ~egfr_value, ~extra,
    "A", 0, 125, 1,
    "A", 1, 124, 2,
    "A", 2, NA,   3,
    "B", 0, 120, 4
  )

  prepared <- prepare_egfr_data(
    df,
    id_col   = "participant",
    time_col = "visit_year",
    egfr_col = "egfr_value",
    min_obs  = 2
  )

  expect_equal(names(prepared), c("id", "time", "egfr"))
  expect_equal(nrow(prepared), 2)           # A has 2 valid, B has 1 -> only A kept
  expect_true(all(prepared$id == "A"))
  expect_true(is.factor(prepared$id))
})

test_that("prepare_egfr_data respects min_obs", {
  df_small <- tibble::tibble(
    id = rep("X", 4),
    time = 0:3,
    egfr = c(125, 124, NA, 123)
  )

  prepared <- prepare_egfr_data(
    df_small,
    id_col   = "id",
    time_col = "time",
    egfr_col = "egfr",
    min_obs  = 4
  )
  expect_equal(nrow(prepared), 0)

  prepared <- prepare_egfr_data(
    df_small,
    id_col   = "id",
    time_col = "time",
    egfr_col = "egfr",
    min_obs  = 3
  )
  expect_equal(nrow(prepared), 3)  # NA removed
})
