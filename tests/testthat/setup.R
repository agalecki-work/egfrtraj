# setup.R - run before tests
library(testthat)
library(egfrtraj)

# Ensure segmented is available
if (!requireNamespace("segmented", quietly = TRUE)) {
  skip("segmented package not available - skipping classification tests")
}
