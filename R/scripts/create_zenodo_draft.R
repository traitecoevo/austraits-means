# Upload a prepared species means package to Zenodo as a new DRAFT deposition.
# It never publishes; review and publish by hand on Zenodo, with the
# maintainers' go-ahead (#30).
#
# Run from the project root, after tar_make() has built the package:
#   ZENODO_TOKEN=... Rscript R/scripts/create_zenodo_draft.R outputs/zenodo/7.0.0 sandbox
# The second argument is "sandbox" (sandbox.zenodo.org, the default) or
# "zenodo" (zenodo.org).

source("R/init.R")

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[[1]] else "outputs/zenodo/7.0.0"
target <- if (length(args) >= 2) args[[2]] else "sandbox"
stopifnot(target %in% c("sandbox", "zenodo"), dir.exists(dir))

draft <- create_zenodo_draft(dir, sandbox = target == "sandbox")
message("Draft created (not published): ", draft$url, " (id ", draft$id, ")")
