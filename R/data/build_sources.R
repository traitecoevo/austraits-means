# Builders for the heavy upstream source artefacts. These are expensive and
# rarely rebuilt (GBIF cleaning ~1.2GB; AusTraits flatten is the slow step), so
# in the targets pipeline (`_targets.R`) the resulting parquets are tracked as
# *file inputs* rather than recomputed on every `tar_make()`. Use these builders
# to (re)generate the cached inputs, then commit the refreshed parquet.
#
# Regeneration (from the project root):
#   source("R/init.R")
#   arrow::write_parquet(build_gbif_climate(),
#                        "data-processed/data_gbif_cleaned.parquet")
#   arrow::write_parquet(build_austraits_climate(),
#                        "data-processed/data_traits_climate.parquet")
# Height is built by `load_data()` (R/data/load_data_height_avh.R); medoids by
# `climate_centres.qmd` (process_taxon_medoid over all taxa); the
# AVH<->GBIF crosswalk by `Data_height_AVH.qmd` (see AGENTS.md).

#' Clean the raw GBIF/ALA occurrences, align names, attach climate, centre.
#'
#' Lifted from the GBIF chunk of the original `dataset_construction.qmd`. Reads
#' the raw GBIF parquet + the aligned-names lookup, extracts climate at each
#' occurrence, drops rogue points, and standardises the climate predictors.
#'
#' @return tibble of cleaned, climate-tagged occurrences (one row per record).
build_gbif_climate <- function(
    gbif_raw = "data/GBIF/2025/GBIF_aus_clean_0001494-250914085247600.parquet",
    gbif_names = "data/GBIF/2025/gbif_names_aligned_2025.parquet",
    rasters = load_climate_rasters()) {

  data_gbif <- arrow::read_parquet(gbif_raw) |>
    janitor::clean_names()

  data_gbif_vect <- data_gbif |>
    terra::vect(
      geom = c("decimal_longitude", "decimal_latitude"),
      crs = terra::crs(rasters)
    )

  data_gbif_climate <-
    terra::extract(x = rasters, y = data_gbif_vect, method = "simple",
                   xy = TRUE, ID = FALSE, dataframe = TRUE) |>
    tibble::as_tibble() |>
    dplyr::select(x, y, dplyr::everything())

  aligned_species_list <- arrow::read_parquet(gbif_names)

  data_gbif |>
    dplyr::rename(original_name = species) |>
    dplyr::select(-taxon_rank) |>
    dplyr::left_join(aligned_species_list, by = "original_name") |>
    dplyr::bind_cols(data_gbif_climate |> dplyr::select(-x, -y)) |>
    dplyr::select(
      original_name, taxon_name = accepted_name, genus, scientific_name,
      dplyr::starts_with("decimal"),
      annual_precipitation = wc2.1_2.5m_bio_12,
      average_temperature = wc2.1_2.5m_bio_1,
      dplyr::everything()
    ) |>
    dplyr::filter(
      !is.na(annual_precipitation), !is.na(average_temperature),
      annual_precipitation > 0,
      !(decimal_longitude == 134.1 & decimal_latitude == -25.7)
    ) |>
    dplyr::mutate(
      log_annual_precipitation = log(annual_precipitation),
      log_annual_precipitation_centered =
        (log_annual_precipitation - mean(log_annual_precipitation)) /
          stats::sd(log_annual_precipitation),
      average_temperature_centered =
        (average_temperature - mean(average_temperature)) /
          stats::sd(average_temperature)
    ) |>
    dplyr::relocate(dplyr::contains("_centered"), .after = "average_temperature")
}

#' Flatten AusTraits and attach climate for the analysis traits.
#'
#' Lifted from the AusTraits chunk of the original `dataset_construction.qmd`.
#' The `flatten_database()` step is the expensive one, hence the file cache.
#'
#' @return tibble of field, georeferenced trait records with climate attached.
build_austraits_climate <- function(
    austraits_rds = "data/data-ausplants/austraits/austraits-7.0.0.rds",
    traits = c("leaf_mass_per_area", "leaf_area", "leaf_N_per_dry_mass",
               "leaf_N_per_area", "leaf_P_per_area", "leaf_P_per_dry_mass",
               "seed_dry_mass", "wood_density"),
    rasters = load_climate_rasters()) {

  austraits <- readRDS(austraits_rds) |>
    austraits::flatten_database()

  austraits |>
    dplyr::filter(
      trait_name %in% traits,
      basis_of_record == "field",
      !is.na(`latitude (deg)`), !is.na(`longitude (deg)`)
    ) |>
    dplyr::mutate(dplyr::across(
      c("latitude (deg)", "longitude (deg)", value),
      ~ suppressWarnings(as.numeric(.x))
    )) |>
    dplyr::filter(!is.na(`latitude (deg)`), !is.na(`longitude (deg)`),
                  !is.na(value)) |>
    dplyr::rename(decimal_longitude = `longitude (deg)`,
                  decimal_latitude = `latitude (deg)`) |>
    extract_climate(rasters = rasters) |>
    dplyr::filter(!is.na(annual_precipitation), !is.na(average_temperature))
}
