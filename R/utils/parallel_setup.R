# furrr/future parallel setup, macOS-safe.
#
# On macOS (especially inside RStudio) forked `multicore` is unstable with
# graphics/ObjC-loaded libraries and is silently disabled, so we always use
# `multisession` (separate R processes). Note `.Rprofile` sometimes sets
# `mc.cores` low, which throttles `availableCores()` and makes `parallelly`
# refuse to start more workers; we read the true core count with
# `methods = "system"` and raise `mc.cores` to match.

#' Set up a multisession future plan.
#'
#' @param workers number of workers; defaults to (system cores - reserve), capped.
#' @param reserve cores to leave free for the OS / RStudio.
#' @param max_workers upper cap on workers.
#' @return (invisibly) the number of workers configured.
setup_parallel <- function(workers = NULL, reserve = 2, max_workers = Inf) {
  avail <- parallelly::availableCores(methods = "system")
  if (is.null(workers)) {
    workers <- max(1, min(max_workers, avail - reserve))
  }
  options(mc.cores = workers)
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    options(future.fork.enable = FALSE)
  }
  future::plan(future::multisession, workers = workers)
  message("Parallel plan: multisession with ", workers,
          " workers (", avail, " system cores available)")
  invisible(workers)
}

#' Reset the future plan back to sequential.
teardown_parallel <- function() {
  future::plan(future::sequential)
  invisible(TRUE)
}
