# Loosen for synthetic noise
slope_p  = 0.05
quad_p   = 0.25
delta_r2 = 0.01

# Uses generate_synthetic_egfr from the package

test_that("stable trajectory is correctly classified", {
  stable <- generate_one_trajectory("S1", "stable", years = 0:10)
  res <- classify_single_trajectory(
    stable,   
    time_col = "time",
    egfr_col = "egfr")

  expect_equal(res$pattern, "stable")
  expect_true(is.na(res$breakpoint_time))
  expect_equal(res$n_knots, 0L)
})

test_that("slow_linear is detected", {
  slow <- generate_one_trajectory("SL1", "slow_linear", years = 0:12)
  res <- classify_single_trajectory(
    slow,
    time_col = "time",
    egfr_col = "egfr")

  expect_equal(res$pattern, "slow_linear")
  expect_true(is.na(res$breakpoint_time))
})

test_that("rapid_with_breakpoint is detected", {
  set.seed(42)
  rapid <- generate_one_trajectory("R1", "rapid_breakpoint", years = 0:20)
  res <- classify_single_trajectory(
    rapid, 
    time_col = "time",
    egfr_col = "egfr",
    slope_p = 0.1,
    quad_p = 0.3,
    delta_r2 = 0.02)

  expect_equal(res$pattern, "rapid_with_breakpoint")
  expect_true(!is.na(res$breakpoint_time))
  expect_true(res$n_knots %in% 1:2)
  expect_true(res$breakpoint_time >= 10 && res$breakpoint_time <= 18)
})

test_that("insufficient data returns correct category", {
  insuff <- generate_one_trajectory("I1", "insufficient", years = 0:26)
  res <- classify_single_trajectory(
    insuff,  
    time_col = "time",
    egfr_col = "egfr")
  expect_equal(res$pattern, "insufficient_data")
})
