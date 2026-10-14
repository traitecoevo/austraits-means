# Upload a prepared species means package to Zenodo as a new DRAFT deposition.
# It never publishes; review and publish by hand on Zenodo, with the
# maintainers' go-ahead (#30).
#
# Run from the project root, after tar_make() has built the package:
#   Rscript R/scripts/create_zenodo_draft.R outputs/zenodo/7.0.0 sandbox [id]
# with ZENODO_TOKEN set (e.g. in ~/.Renviron, never in the repository). The
# second argument is "sandbox" (sandbox.zenodo.org, the default) or "zenodo"
# (zenodo.org); the optional third resumes the draft with that id.

source("R/init.R")

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[[1]] else "outputs/zenodo/7.0.0"
target <- if (length(args) >= 2) args[[2]] else "sandbox"
id <- if (length(args) >= 3) args[[3]] else NULL
stopifnot(target %in% c("sandbox", "zenodo"), dir.exists(dir))

draft <- create_zenodo_draft(dir, sandbox = target == "sandbox", id = id)
message("Draft created (not published): ", draft$url, " (id ", draft$id, ")")
