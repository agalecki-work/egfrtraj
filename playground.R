# ── 1. Restart R session first (strongly recommended to clear everything)
# In RStudio: Session → Restart R (or close and reopen R)

# ── 2. Set working directory to your package folder
setwd("C:/ATG/github/egfrtraj")

# ── 3. Load devtools (install first if missing: install.packages("devtools"))
library(devtools)

# ── 4. Re-document the package (updates docs, NAMESPACE, Rd files)
devtools::document()
### DO NOT usethis::use_readme_rmd() #???
# devtools::build_readme()

# ── 5. Reload the package (makes all functions and current saved data available)
devtools::load_all()

# ── 6. Load required dependencies
library(tidyr)
library(ggplot2)
library(dplyr)
library(segmented)


# ── 8. Re-document after data recreation
devtools::document()

# ── 9. Reload package one more time (now with the fresh data)
devtools::load_all()

# ── 10. Optional: Quick test on one subject (ID "?" for example)

single_df1    <- example_egfr_data %>% filter(id == "1")     # True pattern: linear
single_df6    <- example_egfr_data %>% filter(id == "6")     # true pattern: quadratic
single_df11   <- example_egfr_data %>% filter(id == "11")    # true pattern: segmented
single_df16   <- example_egfr_data %>% filter(id == "16")    # true pattern: insufficient data


single_df <- single_df1  # Data with one subject
traj = classify_single_trajectory(single_df)




traj1 <- extract_and_classify(example_egfr_trajectories, "1") # linear trajectory
traj6 <- extract_and_classify(example_egfr_trajectories, "6") # quadratic
traj11 <- extract_and_classify(example_egfr_trajectories, "11")
traj16 <- extract_and_classify(example_egfr_trajectories, "16")

traj <- traj1

# requirements
$ xlabel should be "Timr (in years)"
# ylabel should be "eGFR (<provide units>)
# if fitted argument is not NULL: Best model name should be included in the subtitle, for example: Best model: linear. Slope = xx (SE = xx)
# Observed data points always black

# linear fit always black. 
# If linear model is the best use solid thick black line, 
# if not the best use dashed thin line. In the subtitle indicate that fitted line (graphic representation) is from <different model> specified in fitted 
@ argumwnt. This will serve as a legend.

# Similar approach for quad fit
# If quad model is the best use solid thick purple line, 
# if not the best use dashed thin line. In the subtitle indicate that fitted line (graphic representation) is from different model specified in fitted 
@ argumwnt. This will serve as a legend.

# Similar approach for seg fit
# If seg model is the best use solid thick red line, 
# if not the best use dashed thin line. In the subtitle indicate that fitted line (graphic representation) is from different model specified in fitted 
@ argumwnt. This will serve as a legend.




"plot(traj) # / Observed values in black, fitted line solid black(because linear is the best)

plot(traj, fitted = c("best")) $ same
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

