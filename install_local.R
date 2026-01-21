
# Make sure you have devtools (once)
# install.packages("devtools")

# Then, from any working directory:
# devtools::install("C:/ATG/github/egfrtraj")
# or using Windows-style backslashes (both usually work)
devtools::install("C:\\ATG\\github\\egfrtraj") # skips tests by default)



library(egfrtraj)
library(dplyr)
library(segmented)
data(example_egfr_trajectories)
summary(example_egfr_trajectories)

df_single <- example_egfr_trajectories %>%
  filter(id == "11")
summary(df_single)
head(df_single)
tail(df_single)


min_df1 <- egfrtraj:::minimal_1traj(
  df_single,
  id_col        = "id",
  time_col      = "time",
  egfr_col      = "egfr"
)

head(min_df1)

min_df2 <- egfrtraj:::minimal_1traj(
  df_single,
  id_col        = "id",
  time_col      = "time",
  egfr_col      = "egfr",
  add_cols      = "true_pattern"
)

head(min_df2)




