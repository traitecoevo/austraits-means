# Build the analysis dataset: trait records + climate + species climate centre,
# with within-group centring on the medoid (Option B; see the decisions log in
# roadmap.md). All trait records, height included, come from one
# austraits.build export.

#' Per-species climate centre (medoid), keyed by trait-data `taxon_name`.
#'
#' The medoids were computed on GBIF names; trait data uses austraits.build
#' names aligned to a different APC release (and rolled up to species by
#' `roll_up_to_species()`). We bridge with the name crosswalk
#' (`build_name_crosswalk()`) and fall back to a direct name match. Where a
#' name has several candidate medoids, the choice is deterministic: a GBIF name
#' identical to the trait name, then one identical to the accepted name, then
#' the first alphabetically; a direct match comes last. `medoid_match` records
#' the route: `name` (same name) or `synonym` (a different GBIF name).
#'
#' @param medoids medoid parquet (one row per GBIF taxon name).
#' @param crosswalk crosswalk table, or a parquet path to one.
#' @return tibble(taxon_name, gbif_taxon_name, medoid_match,
#'   annual_precipitation_centre, average_temperature_centre,
#'   log_annual_precipitation_centre).
build_medoid_lookup <- function(
    medoids = "data-processed/climate_centres.parquet",
    crosswalk) {
  if (is.character(crosswalk)) crosswalk <- arrow::read_parquet(crosswalk)

  med <- arrow::read_parquet(medoids) |>
    dplyr::transmute(
      gbif_taxon_name = taxon_name,
      annual_precipitation_centre = annual_precipitation,
      average_temperature_centre = average_temperature,
      log_annual_precipitation_centre = log(annual_precipitation)
    ) |>
    dplyr::filter(!is.na(annual_precipitation_centre), annual_precipitation_centre > 0)

  centre_cols <- c("taxon_name", "gbif_taxon_name", "medoid_match",
                   "annual_precipitation_centre", "average_temperature_centre",
                   "log_annual_precipitation_centre")

  # Preferred: trait name -> GBIF name (crosswalk) -> medoid
  via_crosswalk <- crosswalk |>
    dplyr::filter(!is.na(gbif_taxon_name)) |>
    dplyr::distinct(taxon_name, gbif_taxon_name, accepted_name) |>
    dplyr::inner_join(med, by = "gbif_taxon_name") |>
    dplyr::mutate(.rank = dplyr::case_when(
      gbif_taxon_name == taxon_name ~ 1L,
      gbif_taxon_name == accepted_name ~ 2L,
      TRUE ~ 3L
    ))

  # Fallback: trait name matches a GBIF/medoid name directly
  direct <- med |>
    dplyr::mutate(taxon_name = gbif_taxon_name, .rank = 4L)

  dplyr::bind_rows(via_crosswalk, direct) |>
    dplyr::arrange(taxon_name, .rank, gbif_taxon_name) |>
    dplyr::distinct(taxon_name, .keep_all = TRUE) |>
    dplyr::mutate(medoid_match = dplyr::if_else(gbif_taxon_name == taxon_name,
                                                "name", "synonym")) |>
    dplyr::select(dplyr::all_of(centre_cols))
}

#' Roll infraspecific taxa up to species (roadmap.md, decision of 14 October
#' 2026).
#'
#' Records of a subspecies, variety or form are relabelled with their species
#' (the build's `binomial`), so the modelled means and the geometric observed
#' means are per species, matching the species-level GBIF names and medoids.
#' The build's own name is kept in `taxon_name_aligned`. Genus-level and other
#' names are left as they are.
#'
#' @param records records with a `taxon_name` column.
#' @param taxa tibble(taxon_name, taxon_rank, binomial): the build's taxa table.
roll_up_to_species <- function(records, taxa) {
  species_of <- taxa |>
    dplyr::filter(taxon_rank %in% c("subspecies", "variety", "form"),
                  !is.na(binomial), binomial != taxon_name) |>
    dplyr::distinct(taxon_name, .species = binomial)

  records |>
    dplyr::left_join(species_of, by = "taxon_name") |>
    dplyr::mutate(taxon_name_aligned = taxon_name,
                  taxon_name = dplyr::coalesce(.species, taxon_name),
                  .after = taxon_name) |>
    dplyr::select(-.species)
}

#' Name-matching rates per trait: species with GBIF records and with a medoid.
#'
#' A species "has GBIF records" if its name reaches a GBIF occurrence name,
#' through the crosswalk or directly; it "has a medoid" if one of those names
#' has a climate centre.
#'
#' @param records trait records (`load_austraits_records()`).
#' @param crosswalk output of `build_name_crosswalk()`.
#' @param gbif_names GBIF occurrence names (`gbif_taxon_names()`).
#' @param medoid_lookup output of `build_medoid_lookup()`.
#' @return one row per trait plus an "all traits" row.
summarise_name_matching <- function(records, crosswalk, gbif_names, medoid_lookup) {
  in_gbif <- union(crosswalk$taxon_name[crosswalk$gbif_taxon_name %in% gbif_names],
                   intersect(records$taxon_name, gbif_names))
  species <- dplyr::distinct(records, trait_name, taxon_name) |>
    dplyr::mutate(has_gbif = taxon_name %in% in_gbif,
                  has_medoid = taxon_name %in% medoid_lookup$taxon_name)
  summarise_rates <- function(d) {
    dplyr::summarise(d,
      n_species = dplyr::n_distinct(taxon_name),
      n_with_gbif = dplyr::n_distinct(taxon_name[has_gbif]),
      n_with_medoid = dplyr::n_distinct(taxon_name[has_medoid]),
      pct_with_medoid = round(100 * n_with_medoid / n_species, 1)
    )
  }
  dplyr::bind_rows(
    species |> dplyr::group_by(trait_name) |> summarise_rates(),
    species |> summarise_rates() |> dplyr::mutate(trait_name = "all traits")
  )
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
