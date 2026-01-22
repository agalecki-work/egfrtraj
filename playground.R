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
df_single <- example_egfr_trajectories %>%
  filter(id == "11") # 1, 6, 11, 16

head(df_single)
tail(df_single)

res <- classify_single_trajectory(
    df_single,
    id_col             = "id",
    time_col           = "time",
    egfr_col           = "egfr",
    bic_tie_threshold  = 4,
    ci_level           = 0.95
)
summary(res)

head(df_single)
tail(df_single)

res <- classify_single_trajectory(
    df_single,
    id_col             = "id",
    time_col           = "time",
    egfr_col           = "egfr",
    bic_tie_threshold  = 4,
    ci_level           = 0.95
)
summary(res)
plot(res)
plot(res, plot_model= "")


names(res)
res$id
res$models
names(res$trajectory)

model_summ <- res$models
names(model_summ)
model_summ
as_tibble(model_summ$info[[1]]) # unnested $info
as_tibble(model_summ$linear[[1]])# unnested $linear
as_tibble(model_summ$quadratic[[1]])# unnested $quadratic
as_tibble(model_summ$segmented[[1]])# unnested $segmented

traj_all <- classify_multiple_trajectories(
  example_egfr_trajectories,
  id_col             = "id",
  time_col           = "time",
  egfr_col           = "egfr",
  bic_tie_threshold  = 4,
  ci_level           = 0.95
)



traj_all <- classify_multiple_trajectories(example_egfr_trajectories,
    id_col             = "id",
    time_col           = "time",
    egfr_col           = "egfr",
    bic_tie_threshold  = 4,
    ci_level           = 0.95
)



