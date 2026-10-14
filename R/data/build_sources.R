# Builders for the upstream source artefacts.
#
# The GBIF occurrences (~1.2GB) are expensive and rarely rebuilt, so in the
# targets pipeline (`_targets.R`) the resulting parquet is tracked as a *file
# input* rather than recomputed on every `tar_make()`. Regenerate it from the
# project root with
#   source("R/init.R")
#   arrow::write_parquet(build_gbif_climate(),
#                        "data-processed/data_gbif_cleaned.parquet")
# Medoids are built by `climate_centres.qmd` (process_taxon_medoid over all
# taxa); the name crosswalk by `Data_height_AVH.qmd` (see AGENTS.md).
#
# Trait records, height included, come from an austraits.build export
# (`data-raw/austraits/README.md`); `build_austraits_climate()` runs inside the
# pipeline, because flattening a build takes seconds.

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

#' Version string and commit of an austraits.build export.
#'
#' @param austraits_rds an rds written by austraits.build (`build.R`).
#' @return list(austraits_version, austraits_build_commit).
austraits_build_info <- function(austraits_rds) {
  info <- readRDS(austraits_rds)$build_info
  list(austraits_version = info$version,
       austraits_build_commit = info$git_SHA)
}

#' Flatten AusTraits and attach climate for the analysis traits.
#'
#' Reads an austraits.build export (a development build since #12; it carries
#' the AVH height dataset `AVH_2026`), keeps field records of the analysis
#' traits that have coordinates, and attaches WorldClim climate. Value-type
#' rules are applied later, in `load_austraits_records()`.
#'
#' @param austraits_rds the austraits.build export.
#' @param traits trait names to keep.
#' @param rasters climate rasters (`load_climate_rasters()`).
#' @return tibble of field, georeferenced trait records with climate attached.
build_austraits_climate <- function(
    austraits_rds,
    traits = c("plant_height", "leaf_mass_per_area", "leaf_area",
               "leaf_N_per_dry_mass", "leaf_N_per_area", "leaf_P_per_area",
               "leaf_P_per_dry_mass", "seed_dry_mass", "wood_density"),
    rasters = load_climate_rasters()) {

  austraits <- readRDS(austraits_rds) |>
    austraits::flatten_database()

  austraits |>
    dplyr::filter(
      trait_name %in% traits,
      basis_of_record == "field",
      !is.na(`latitude (deg)`), !is.na(`longitude (deg)`)
    ) |>
    dplyr::select(
      dataset_id, taxon_name, observation_id, trait_name, value, unit,
      value_type, basis_of_record, life_stage, location_id, collection_date,
      original_name, taxon_rank, genus, family,
      `latitude (deg)`, `longitude (deg)`
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
