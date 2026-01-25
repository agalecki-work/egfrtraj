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
library(tidyr)
library(ggplot2)
library(dplyr)
library(segmented)


# ── 7. RECREATE the example data from scratch (just in case)
set.seed(20250116)  # for reproducibility

example_egfr_trajectories <- generate_synthetic_egfr(
  n_per_pattern = 5,       # 5 subjects per pattern
  years         = 0:26,    # 27 rows per person
  baseline_egfr = 125,
  noise_sd      = 2
)

# Quick verification (should show  4 patterns * 5 subjects per pattern -> 20 subjects, 27 per subject)
nrow(example_egfr_trajectories) #  27*20 = 540 rows in total
table(example_egfr_trajectories$true_pattern) # 5 *27 = 135 rows per pattern 

# Save the freshly generated data back to data/ folder (overwrites old version)
usethis::use_data(example_egfr_trajectories, overwrite = TRUE)

# ── 8. Re-document after data recreation
devtools::document()

# ── 9. Reload package one more time (now with the fresh data)
devtools::load_all()

# ── 10. Optional: Quick test on one subject (ID "?" for example)

traj1 <- extract_and_classify(example_egfr_trajectories, "1") # linear trajectory
traj6 <- extract_and_classify(example_egfr_trajectories, "6") # quadratic
traj11 <- extract_and_classify(example_egfr_trajectories, "11")
traj16 <- extract_and_classify(example_egfr_trajectories, "16")

traj <- traj1



plot(traj) # Observed values in black

plot(traj, fitted = c("best")) 
plot(traj, fitted = c("best"), show_ci = FALSE) 

plot(traj, fitted = c("linear"))
plot(traj, fitted = c("linear"), show_ci = FALSE)

plot(traj, fitted = c("quad"))
plot(traj, fitted = c("seg"))



plot(traj, fitted = c("best", "quad")) # fitted line (use different color, CI included. In the title include "Best model: linear. slope = xx(xx, xx)"
plot(traj, fitted = c("best", "quad", "seg")) # fitted line in black (use different color, CI included. In the title include "Best model: linear. slope = xx(xx, xx)"

plot(traj, fitted = c("quad")) # fitted quadratic function (different color, dashed line, good idea) , CI included
plot(traj, fitted = c("seg")) # fitted segmented (dashed line red), CI included add vertical dotted (descrete) vertical line for breakpoint and CI on x-axis)















# Extract multiple subjects from a dataset

trajx <- extract_and_classify(example_egfr_trajectories, c("1","6","11","16"))


# inspect_egfr_traj(traj)

