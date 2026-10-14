# Build the analysis dataset: trait records + climate + species climate centre,
# with within-group centring on the medoid (Option B; see the decisions log in
# roadmap.md). All trait records, height included, come from one
# austraits.build export.

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

# Common record-level schema of the unified table.
.unified_cols <- c(
  "source", "dataset_id", "taxon_name", "original_name", "family", "genus",
  "trait_name", "value", "log_value", "unit", "value_type",
  "decimal_longitude", "decimal_latitude", "state", "year",
  "annual_precipitation", "average_temperature", "log_annual_precipitation"
)

#' Per-dataset value-type rules (roadmap.md, Invariants).
#'
#' A dataset listed here enters the means with the value types given here
#' instead of the default raw/mean/median.
#'
#' `AVH_2026` (herbarium specimen heights in austraits.build) carries one bound
#' per specimen, a minimum or a maximum read from the label, and no raw values,
#' so the default filter would drop all of it. Each bound is used as one
#' observation, as the earlier AVH CSV loader did.
#'
#' @return tibble(dataset_id, value_type), one row per allowed combination.
dataset_value_type_rules <- function() {
  tibble::tribble(
    ~dataset_id, ~value_type,
    "AVH_2026",  "minimum",
    "AVH_2026",  "maximum"
  )
}

#' Field trait records (uncentred), in the unified schema.
#'
#' Keeps field records with a positive value and climate, and applies the
#' value-type rule: datasets in `rules` keep their listed value types, all
#' others keep `value_types`. Height comes from `AVH_2026` like any other
#' dataset, so `dataset_id` is `AVH_2026` for all its records (the earlier
#' per-collector grouping is not carried by the austraits.build dataset).
#'
#' @param records output of `build_austraits_climate()`, or a parquet path.
#' @param value_types default value types kept.
#' @param rules per-dataset exceptions (`dataset_value_type_rules()`).
load_austraits_records <- function(
    records,
    value_types = c("raw", "mean", "median"),
    rules = dataset_value_type_rules()) {
  if (is.character(records)) records <- arrow::read_parquet(records)
  ruled <- paste(rules$dataset_id, rules$value_type)

  records |>
    dplyr::filter(
      basis_of_record == "field",
      dplyr::if_else(dataset_id %in% rules$dataset_id,
                     paste(dataset_id, value_type) %in% ruled,
                     value_type %in% value_types),
      !is.na(value), value > 0,
      !is.na(annual_precipitation), !is.na(average_temperature)
    ) |>
    dplyr::mutate(
      source = "AusTraits",
      log_value = log(value),
      log_annual_precipitation = log(annual_precipitation),
      year = suppressWarnings(as.integer(substr(collection_date, 1, 4)))
    ) |>
    assign_state_fast() |>
    dplyr::select(dplyr::any_of(.unified_cols))
}

#' Build the unified analysis dataset.
#'
#' Attaches the medoid centre to the trait records with Option-B centring and
#' drops species with no medoid. Climate centring is per species (shared
#' across traits).
#'
#' @param records output of `load_austraits_records()`.
#' @param medoid_lookup output of `build_medoid_lookup()`.
#' @param traits optional character vector to restrict to; NULL keeps all.
build_dataset <- function(records, medoid_lookup = build_medoid_lookup(),
                          traits = NULL) {
  if (!is.null(traits)) records <- dplyr::filter(records, trait_name %in% traits)
  add_medoid_centring(records, medoid_lookup)
}
