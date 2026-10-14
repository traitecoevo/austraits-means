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
# taxa). The name crosswalk (`build_name_crosswalk()`) runs inside the
# pipeline, with the APC release pinned in `_targets.R`.
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

#' Name crosswalk from trait-record names to GBIF names, via one APC release.
#'
#' The trait records (austraits.build) and the GBIF occurrences, and so the
#' medoids, were name-aligned against different APC releases, so a species can
#' carry different names in the two (#1). Both name sets are aligned to the
#' same pinned APC release with APCalign and matched on the accepted name.
#' Names that do not align to an accepted name are left out; such names can
#' still match a GBIF name directly in `build_medoid_lookup()`.
#'
#' @param taxon_names trait-record taxon names (all taxa in the build).
#' @param gbif_names taxon names in the GBIF occurrences (the medoid keys).
#' @param apc_version APC/APNI release to align against (an APCalign release
#'   tag, e.g. "2026-07-30").
#' @return tibble(taxon_name, accepted_name, gbif_taxon_name, apc_version); a
#'   trait name can map to more than one GBIF name (synonyms in GBIF).
build_name_crosswalk <- function(taxon_names, gbif_names, apc_version) {
  resources <- APCalign::load_taxonomic_resources(version = apc_version,
                                                  quiet = TRUE)
  align_accepted <- function(names) {
    APCalign::create_taxonomic_update_lookup(
      sort(unique(stats::na.omit(names))), resources = resources,
      quiet = TRUE
    ) |>
      dplyr::select(original_name, accepted_name = suggested_name) |>
      dplyr::filter(!is.na(accepted_name))
  }

  trait <- align_accepted(taxon_names) |> dplyr::rename(taxon_name = original_name)
  gbif <- align_accepted(gbif_names) |> dplyr::rename(gbif_taxon_name = original_name)

  trait |>
    dplyr::inner_join(gbif, by = "accepted_name", relationship = "many-to-many") |>
    dplyr::mutate(apc_version = apc_version) |>
    dplyr::arrange(taxon_name, gbif_taxon_name)
}

#' Taxon names in the GBIF occurrence parquet.
gbif_taxon_names <- function(gbif_path) {
  arrow::open_dataset(gbif_path) |>
    dplyr::distinct(taxon_name) |>
    dplyr::collect() |>
    dplyr::pull(taxon_name) |>
    stats::na.omit() |>
    sort()
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
