# Build the analysis dataset: trait records + climate + species climate centre,
# with within-group centring on the medoid (Option B; see refactor_plan.md §3.2).
#
# Phase 1 covers height (AVH). The same machinery generalizes to the unified
# multi-trait table in Phase 2 (height becomes one trait among many).

#' Per-species climate centre (medoid), keyed by trait-data `taxon_name`.
#'
#' The medoids were computed on GBIF accepted names; trait data uses AVH/AusTraits
#' names aligned to a different APC release. We bridge with the AVH<->GBIF
#' crosswalk and fall back to a direct name match for names absent from it.
#' Result has one centre per `taxon_name` (crosswalk preferred over direct).
#'
#' @return tibble(taxon_name, annual_precipitation_centre,
#'   average_temperature_centre, log_annual_precipitation_centre).
build_medoid_lookup <- function(
    medoids = "data-processed/climate_centres.parquet",
    crosswalk = "data-processed/gbif_avh_name_crosswalk.parquet") {

  med <- arrow::read_parquet(medoids) |>
    dplyr::transmute(
      gbif_taxon_name = taxon_name,
      annual_precipitation_centre = annual_precipitation,
      average_temperature_centre = average_temperature,
      log_annual_precipitation_centre = log(annual_precipitation)
    ) |>
    dplyr::filter(!is.na(annual_precipitation_centre), annual_precipitation_centre > 0)

  centre_cols <- c("taxon_name", "annual_precipitation_centre",
                   "average_temperature_centre", "log_annual_precipitation_centre")

  # Preferred: AVH name -> GBIF name (crosswalk) -> medoid
  via_crosswalk <-
    arrow::read_parquet(crosswalk) |>
    dplyr::select(taxon_name, gbif_taxon_name) |>
    dplyr::filter(!is.na(gbif_taxon_name)) |>
    dplyr::distinct() |>
    dplyr::inner_join(med, by = "gbif_taxon_name") |>
    dplyr::select(dplyr::all_of(centre_cols))

  # Fallback: trait name matches a GBIF/medoid name directly
  direct <- med |>
    dplyr::mutate(taxon_name = gbif_taxon_name) |>
    dplyr::select(dplyr::all_of(centre_cols))

  # Crosswalk wins when both exist; one centre per taxon_name
  dplyr::bind_rows(via_crosswalk, direct) |>
    dplyr::group_by(taxon_name) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
}

#' Attach the medoid centre to records and centre the climate predictors on it.
#'
#' Option B centring: `<var>_mn` = medoid, `<var>_c` = obs - medoid. Records whose
#' species has no medoid are dropped (inner join) and reported.
#'
#' @param data record-level data with `log_annual_precipitation`,
#'   `average_temperature`, `taxon_name`.
#' @param medoid_lookup output of `build_medoid_lookup()` (built by default).
add_medoid_centring <- function(data, medoid_lookup = build_medoid_lookup()) {
  n_spp_in <- dplyr::n_distinct(data$taxon_name)

  out <- data |>
    # drop any stale centring columns from upstream caches
    dplyr::select(-dplyr::any_of(c(
      "log_annual_precipitation_mn", "average_temperature_mn",
      "log_annual_precipitation_c", "average_temperature_c",
      "decimal_latitude_mn", "decimal_longitude_mn", "log_value_c"
    ))) |>
    dplyr::inner_join(medoid_lookup, by = "taxon_name") |>
    dplyr::mutate(
      log_annual_precipitation_mn = log_annual_precipitation_centre,
      average_temperature_mn = average_temperature_centre,
      log_annual_precipitation_c = log_annual_precipitation - log_annual_precipitation_mn,
      average_temperature_c = average_temperature - average_temperature_mn
    )

  n_spp_out <- dplyr::n_distinct(out$taxon_name)
  message("Medoid centring: kept ", n_spp_out, " / ", n_spp_in,
          " species (dropped ", n_spp_in - n_spp_out, " with no medoid)")
  out
}

# Common record-level schema shared by all sources.
.unified_cols <- c(
  "source", "dataset_id", "taxon_name", "original_name", "family", "genus",
  "trait_name", "value", "log_value", "unit", "value_type",
  "decimal_longitude", "decimal_latitude", "state", "year",
  "annual_precipitation", "average_temperature", "log_annual_precipitation"
)

#' Load AVH height records (uncentred), in the unified schema.
#'
#' Height is AVH-only (`trait_name = "plant_height"`); `dataset_id` is one per
#' AVH collector (placeholder: institution code, refined in round 3).
load_height_records <- function(
    height_climate = "data-processed/data_height_climate.parquet") {
  arrow::read_parquet(height_climate) |>
    dplyr::filter(!is.na(value), value > 0) |>
    dplyr::mutate(
      source = "AVH",
      trait_name = "plant_height",
      unit = "m",
      log_value = log(value),
      dataset_id = paste0("AVH:", dplyr::coalesce(institutionCode, "unknown"))
    ) |>
    assign_state_fast() |>
    dplyr::select(dplyr::any_of(.unified_cols))
}

#' Load AusTraits multi-trait field records (uncentred), in the unified schema.
#'
#' @param value_types kept value types (field raw/mean/median; flora min/max go
#'   to the round-2 means path).
load_austraits_records <- function(
    austraits_climate = "data-processed/data_traits_climate.parquet",
    value_types = c("raw", "mean", "median")) {
  arrow::read_parquet(austraits_climate) |>
    dplyr::filter(
      basis_of_record == "field",
      value_type %in% value_types,
      !is.na(value), value > 0,
      !is.na(annual_precipitation), !is.na(average_temperature)
    ) |>
    dplyr::mutate(
      source = "AusTraits",
      log_value = log(value),
      log_annual_precipitation = log(annual_precipitation)
    ) |>
    assign_state_fast() |>
    dplyr::select(dplyr::any_of(.unified_cols))
}

#' Build the unified analysis dataset (AVH height + AusTraits traits).
#'
#' Binds the sources into one long table (one row per trait record), attaches
#' the medoid centre with Option-B centring, and drops species with no medoid.
#' Climate centring is per species (shared across traits).
#'
#' @param traits optional character vector to restrict to; NULL keeps all.
#' @param sources which sources to include.
#' @param medoid_lookup output of `build_medoid_lookup()`.
build_dataset <- function(traits = NULL,
                          sources = c("AVH", "AusTraits"),
                          medoid_lookup = build_medoid_lookup()) {
  sources <- match.arg(sources, several.ok = TRUE)
  parts <- list()
  if ("AVH" %in% sources) parts$avh <- load_height_records()
  if ("AusTraits" %in% sources) parts$austraits <- load_austraits_records()

  data <- dplyr::bind_rows(parts)
  if (!is.null(traits)) data <- dplyr::filter(data, trait_name %in% traits)

  add_medoid_centring(data, medoid_lookup)
}

#' Convenience wrapper: the height-only dataset.
build_dataset_height <- function(medoid_lookup = build_medoid_lookup()) {
  build_dataset(traits = "plant_height", sources = "AVH",
                medoid_lookup = medoid_lookup)
}
