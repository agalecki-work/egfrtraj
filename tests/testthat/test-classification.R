test_that("classify_egfr_trajectories works end-to-end with synthetic data", {
  sim_data <- generate_synthetic_egfr(n_per_pattern = 2, seed = 42)

  result <- classify_egfr_trajectories(
    sim_data,
    id_col   = "id",
    time_col = "time",
    egfr_col = "egfr",
    min_obs  = 5,
    slope_p  = 0.10,
    quad_p   = 0.30,
    delta_r2 = 0.02,
    verbose  = FALSE
  )

  expect_s3_class(result, "egfr_trajectories")
  expect_true("results" %in% names(result))
  expect_true(nrow(result$results) > 0)

  patterns <- result$results$pattern
  expect_true("stable" %in% patterns)
  expect_true("slow_linear" %in% patterns)
  expect_true(any(grepl("rapid_with_breakpoint", patterns)))
  expect_true(any(grepl("insufficient_data", patterns)))
})

test_that("classification respects parameter changes", {
  sim_data <- generate_synthetic_egfr(n_per_pattern = 3, seed = 123)

  # Very strict → should get more "slow_linear"
  strict <- classify_egfr_trajectories(
    sim_data,
    id_col   = "id",
    time_col = "time",
    egfr_col = "egfr",
    slope_p  = 0.10,    
    quad_p   = 0.30,
    delta_r2 = 0.02, 
    verbose  = FALSE
)
  expect_true(sum(strict$results$pattern == "slow_linear") >= 8)

 # Lenient → more rapid detected
 lenient <- classify_egfr_trajectories(
   sim_data,
   id_col   = "id",
   time_col = "time",
   egfr_col = "egfr",
   slope_p  = 0.10,   
   quad_p   = 0.30,
   delta_r2 = 0.02,  
   verbose  = FALSE
)
  expect_true(sum(lenient$results$pattern == "rapid_with_breakpoint") >= 3)
})
