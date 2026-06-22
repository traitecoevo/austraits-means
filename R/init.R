# Project function library loader.
#
# Source this once at the top of a notebook (`source("R/init.R")`) to load every
# function in the R/ library. Run-scripts in R/scripts/ are intentionally NOT
# sourced here (they execute work on source); run those directly with Rscript.
#
# Sourcing only defines functions, so package attach order does not matter —
# load the packages each notebook needs separately.

local({
  lib_dirs <- c("utils", "data", "medoid", "model", "process", "plot", "means", "panels")
  files <- list.files(
    file.path("R", lib_dirs),
    pattern = "\\.R$", full.names = TRUE, recursive = TRUE
  )
  invisible(lapply(files, source))
})
