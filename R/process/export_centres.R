# Export a per-trait species-centres table + a data-dictionary, mirroring the
# height deliverable (one row per species, gzip-compressed CSV).

#' Column definitions for the species-centres export.
species_centres_meta <- function(centres) {
  defs <- c(
    taxon_name = "species name",
    family = "plant family",
    genus = "genus",
    n = "number of trait records used for the species",
    log_annual_precipitation_mn = "log mean annual precipitation (mm) at the species climate centre (medoid)",
    average_temperature_mn = "mean annual temperature (degC) at the species climate centre (medoid)",
    annual_precipitation_centre = "mean annual precipitation (mm) at the medoid",
    average_temperature_centre = "mean annual temperature (degC) at the medoid",
    log_annual_precipitation_centre = "log mean annual precipitation (mm) at the medoid",
    log_estimate = "predicted log trait value at the climate centre",
    log_conf.low = "lower bound on predicted log trait value",
    log_conf.high = "upper bound on predicted log trait value",
    estimate = "predicted trait value at the climate centre (original units)",
    conf.low = "lower bound on predicted trait value (original units)",
    conf.high = "upper bound on predicted trait value (original units)",
    g_0 = "species random intercept",
    g_log_annual_precipitation_c = "species random slope for within-species precipitation",
    g_average_temperature_c = "species random slope for within-species temperature"
  )
  tibble::tibble(
    name = names(centres),
    class = vapply(centres, \(x) class(x)[1], character(1)),
    definition = unname(defs[names(centres)])
  )
}

#' Write `<trait>_climate_response.csv.gz` + `_meta.csv` to `out_dir`.
#'
#' @return (invisibly) the data CSV path.
export_species_centres <- function(centres, trait, out_dir = "outputs/traits") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(out_dir, paste0(clean_species_filename(trait), "_climate_response"))
  readr::write_csv(centres, paste0(stem, ".csv.gz"))
  readr::write_csv(species_centres_meta(centres), paste0(stem, "_meta.csv"))
  invisible(paste0(stem, ".csv.gz"))
}
