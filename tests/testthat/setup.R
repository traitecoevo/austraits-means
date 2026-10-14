# Load the project library from the repository root (R/init.R uses paths
# relative to the root). Run the tests from the root with
#   Rscript -e 'testthat::test_dir("tests/testthat")'
project_root <- normalizePath(file.path("..", ".."))
withr::with_dir(project_root, source("R/init.R"))
