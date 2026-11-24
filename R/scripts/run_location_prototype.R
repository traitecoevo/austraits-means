# Location model prototype (#41): fit the pooled woodiness model (species effects only) and the
# location model for each trait, and write the comparison to
# outputs/local_range/prototype_<cell_km>km.rds, read by
# docs/local-range-model.qmd. Not part of the product targets.
#
# Usage, from the project root:
#   Rscript R/scripts/run_location_prototype.R [cell_km] [trait ...]
# Defaults: 50 km cells, every modelled trait, height last (the slow one).

source("R/init.R")
args <- commandArgs(trailingOnly = TRUE)
cell_km <- if (length(args) >= 1) as.numeric(args[[1]]) else 50
traits <- if (length(args) >= 2) args[-1] else
  c("leaf_mass_per_area", "leaf_area", "leaf_N_per_dry_mass", "leaf_N_per_area",
    "leaf_P_per_area", "leaf_P_per_dry_mass", "seed_dry_mass", "wood_density",
    "plant_height")

unified <- targets::tar_read(unified)
woody_lookup <- targets::tar_read(woody_lookup)
dir.create("outputs/local_range", showWarnings = FALSE, recursive = TRUE)
out_file <- sprintf("outputs/local_range/prototype_%gkm.rds", cell_km)
res <- if (file.exists(out_file)) readRDS(out_file) else list()

for (trait in traits) {
  message(Sys.time(), "  ", trait)
  res[[trait]] <- run_location_prototype(unified, trait, woody_lookup, cell_km = cell_km)
  saveRDS(res, out_file)  # after each trait, so a long height fit loses nothing
  r <- res[[trait]]
  if (!is.null(r)) message(sprintf("  %s: %.0f s, converged %s, AIC %.0f -> %.0f",
                                   trait, r$seconds, r$converged, r$aic[[1]], r$aic[[2]]))
}

# Diagnostic figures: records, species means and modelled range in 12 cells
for (r in res) {
  if (is.null(r$diagnostics)) next
  save_figure(plot_community_range_cells(r),
              sprintf("outputs/local_range/figures/%s_cells_%gkm", r$trait, cell_km),
              width = 11, height = 6.5)
  if (!is.null(r$spread))
    save_figure(plot_spread_by_climate(r),
                sprintf("outputs/local_range/figures/%s_spread_%gkm", r$trait, cell_km),
                width = 11, height = 6.5)
}
