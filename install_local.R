
# Make sure you have devtools (once)
# install.packages("devtools")

# Then, from any working directory:
# devtools::install("C:/ATG/github/egfrtraj")
# or using Windows-style backslashes (both usually work)
devtools::install("C:\\ATG\\github\\egfrtraj") # skips tests by default

library(conflicted)
# library(dplyr)
# library(tidyr)     # ← this brings pivot_longer() and pivot_wider()
library(tidyverse)
library(egfrtraj)
library(segmented)
#  conflicts()
prefer_tidyverbs()

sort(getNamespaceExports("egfrtraj"))
ls("package:egfrtraj")

# ?prefer_tidyverbs
# ?check_long_structure
# ?generate_synthetic_egfr

data(example_egfr_data)
colnames(example_egfr_data)

# ?check_long_structure
# ?extract_egfrtibble
res <- extract_egfrtibble(example_egfr_data)

# Quick inspection
class(res)                                 # "egfrtibble" + tibble classes
dim(res)                                   # rows × 3 (id, time, egfr)
attr(res, "origdata_info")$final_subjects  # usually "All"

# ?subset_egfrtibble

res <- extract_egfrtibble(example_egfr_data)

# Get all unique IDs from the original data
all_ids <- unique(example_egfr_data$id)




