# Continuous spread (#41): add the RTMB model with the species SD a smooth
# function of the species' climate centre (R/model/spread_rtmb.R) to the
# location prototype result, outputs/local_range/prototype_<cell_km>km.rds,
# under `spread$continuous`, and draw the spread figures. Run after
# R/scripts/run_location_prototype.R.
#
# Usage, from the project root:
#   Rscript R/scripts/run_spread_rtmb.R [cell_km] [trait ...]

source("R/init.R")
args <- commandArgs(trailingOnly = TRUE)
cell_km <- if (length(args) >= 1) as.numeric(args[[1]]) else 50
res_file <- sprintf("outputs/local_range/prototype_%gkm.rds", cell_km)
res <- readRDS(res_file)
traits <- if (length(args) >= 2) args[-1] else names(res)
# height last: the slow one
traits <- c(setdiff(traits, "plant_height"), intersect(traits, "plant_height"))

unified <- add_location_key(targets::tar_read(unified), "cell", cell_km = cell_km)
woody_lookup <- targets::tar_read(woody_lookup)

for (trait in traits) {
  if (is.null(res[[trait]])) next
  message(Sys.time(), "  ", trait)
  m <- fit_spread_rtmb(unified, trait, woody_lookup)
  vc_new <- res[[trait]]$varcomp_new
  res[[trait]]$spread$continuous <- list(
    aic = m$aic, converged = m$converged, seconds = m$seconds, message = m$message,
    theta = m$spread$theta, cov = m$spread$cov, varcomp = m$vc,
    curves = dplyr::bind_rows(spread_curve(m, "precipitation"), spread_curve(m, "temperature")),
    benchmark = summarise_community_benchmark(benchmark_community_variance(m, noise_vc = vc_new)),
    cells = community_diagnostics(m)$cells)
  res[[trait]]$aic[["spread_continuous"]] <- m$aic
  saveRDS(res, res_file)
  message(sprintf("  %s: %.0f s, converged %s, AIC location %.0f -> continuous spread %.0f",
                  trait, m$seconds, m$converged, res[[trait]]$aic[["location"]], m$aic))
}

for (r in res) {
  if (is.null(r$spread$continuous)) next
  save_figure(plot_spread_by_climate(r),
              sprintf("outputs/local_range/figures/%s_spread_%gkm", r$trait, cell_km),
              width = 11, height = 6.5)
  save_figure(plot_spread_curves(r),
              sprintf("outputs/local_range/figures/%s_spread_curve_%gkm", r$trait, cell_km),
              width = 9, height = 6)
}
