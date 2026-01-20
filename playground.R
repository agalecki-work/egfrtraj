# ── 1. Restart R session first (strongly recommended to clear everything)
# In RStudio: Session → Restart R (or close and reopen R)

# ── 2. Set working directory to your package folder
setwd("C:/ATG/github/egfrtraj")

# ── 3. Load devtools (install first if missing: install.packages("devtools"))
library(devtools)

# ── 4. Re-document the package (updates docs, NAMESPACE, Rd files)
devtools::document()

# ── 5. Reload the package (makes all functions and current saved data available)
devtools::load_all()

# ── 6. Load required dependencies
library(dplyr)
library(segmented)

# ── 7. RECREATE the example data from scratch (just in case)
set.seed(20250116)  # for reproducibility

example_egfr_trajectories <- generate_synthetic_egfr(
  n_per_pattern = 1,       # one trajectory per pattern → 108 rows total
  years         = 0:26,
  baseline_egfr = 125,
  noise_sd      = 2
)

# Quick verification (should show 108 rows, 27 per pattern)
nrow(example_egfr_trajectories)
table(example_egfr_trajectories$true_pattern)

# Save the freshly generated data back to data/ folder (overwrites old version)
usethis::use_data(example_egfr_trajectories, overwrite = TRUE)

# ── 8. Re-document after data recreation
devtools::document()

# ── 9. Reload package one more time (now with the fresh data)
devtools::load_all()

# ── 10. Optional: Quick test on one subject (ID "3" for example)
df_single <- example_egfr_trajectories %>%
  filter(id == "3")

fit <- classify_single_trajectory(
  df_single,
  id_col        = "id",
  time_col      = "time",
  egfr_col      = "egfr",
  estimate_se   = TRUE
)

# Quick results check
fit$pattern
fit$breakpoint_time
fit$data %>% head()


