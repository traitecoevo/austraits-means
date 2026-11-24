# Continuous spread (#41): add the RTMB models with the species SD a smooth
# function of the species' climate centre (R/model/spread_rtmb.R) to the
# location prototype result, outputs/local_range/prototype_<cell_km>km.rds:
# the free quadratic under `spread$continuous`, the shrunk one (from the
# free one, fit_spread_shrunk()) under `spread$shrunk`. Fits already in the
# file are kept unless --refit is given. Then draw the spread and cell figures. Run after
# R/scripts/run_location_prototype.R.
#
# Usage, from the project root:
#   Rscript R/scripts/run_spread_rtmb.R [--refit] [cell_km] [trait ...]

source("R/init.R")
args <- commandArgs(trailingOnly = TRUE)
refit <- "--refit" %in% args
args <- setdiff(args, "--refit")
cell_km <- if (length(args) >= 1) as.numeric(args[[1]]) else 50
res_file <- sprintf("outputs/local_range/prototype_%gkm.rds", cell_km)
res <- readRDS(res_file)
traits <- if (length(args) >= 2) args[-1] else names(res)
# height last: the slow one
traits <- c(setdiff(traits, "plant_height"), intersect(traits, "plant_height"))

unified <- add_location_key(targets::tar_read(unified), "cell", cell_km = cell_km)
woody_lookup <- targets::tar_read(woody_lookup)

summarise_spread_fit <- function(m, vc_new) list(
  aic = m$aic, converged = m$converged, pd_hessian = m$pd_hessian,
  seconds = m$seconds, message = m$message,
  terms = m$spread$terms, theta = m$spread$theta, cov = m$spread$cov,
  shrinkage = m$spread$shrinkage, varcomp = m$vc,
  curves = dplyr::bind_rows(spread_curve(m, "precipitation"), spread_curve(m, "temperature")),
  benchmark = summarise_community_benchmark(benchmark_community_variance(m, noise_vc = vc_new)),
  cells = community_diagnostics(m)$cells)

for (trait in traits) {
  if (is.null(res[[trait]])) next
  for (fit in c("continuous", "shrunk")) {
    if (!refit && !is.null(res[[trait]]$spread[[fit]])) next
    message(Sys.time(), "  ", trait, ", ", fit)
    m <- if (fit == "continuous") fit_spread_rtmb(unified, trait, woody_lookup) else
      fit_spread_shrunk(unified, trait, woody_lookup, free = res[[trait]]$spread$continuous)
    res[[trait]]$spread[[fit]] <- summarise_spread_fit(m, res[[trait]]$varcomp_new)
    res[[trait]]$aic[[paste0("spread_", fit)]] <- m$aic
    saveRDS(res, res_file)
    message(sprintf("  %s, %s: %.0f s, converged %s, AIC location %.0f -> %.0f%s",
                    trait, fit, m$seconds, m$converged, res[[trait]]$aic[["location"]], m$aic,
                    if (fit == "shrunk") paste0("; tau ", paste(signif(m$spread$shrinkage$tau, 2),
                                                                collapse = ", ")) else ""))
  }
}

for (r in res) {
  if (is.null(r$spread$continuous)) next
  save_figure(plot_spread_by_climate(r),
              sprintf("outputs/local_range/figures/%s_spread_%gkm", r$trait, cell_km),
              width = 11, height = 6.5)
  save_figure(plot_community_range_cells(r),
              sprintf("outputs/local_range/figures/%s_cells_%gkm", r$trait, cell_km),
              width = 11, height = 6.5)
  save_figure(plot_spread_curves(r),
              sprintf("outputs/local_range/figures/%s_spread_curve_%gkm", r$trait, cell_km),
              width = 9, height = 6)
}
