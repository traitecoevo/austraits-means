# Observed means, portal method (`estimate_type = observed_mean`,
# `method = portal_arithmetic`): the species means 7.0.0 product.
#
# A port of the species averages the AusTraits data portal computes offline,
# taken from `austraits.portal` at commit 13479b9 (the commit that built the
# deployed AusTraits 7.0.0 data): `estimate_species_trait_means()` and its
# helpers in `R/species_means.R`, and the steps of `prepare_data_for_portal()`
# and `data-raw/prepare_data.R` that feed it. The logic is copied as it was,
# with only namespaces made explicit and functions renamed. Its quirks are the
# product, so they are kept, not fixed (roadmap.md, Invariants):
#
#   - an arithmetic mean of location means and flora midpoints;
#   - `value_geom_mean` equals `value_mean`, because it is computed after
#     `value_mean` has been overwritten inside the same `summarise()`;
#   - no basis_of_record or life_stage filter on raw/mean/median records;
#   - a lone flora bound (a minimum without its maximum, or the reverse) is
#     used as the midpoint;
#   - no roll-up of infraspecific or genus-level names.
#
# Fixes go in the geometric method (observed_means_geometric.R).

#' Portal commit whose method this file reproduces.
portal_method_commit <- "13479b9"

#' MD5 of each AusTraits release rds on Zenodo, checked before use.
austraits_release_md5 <- c(
  "7.0.0" = "c040dada5776177a1fabaced7ac0afdb"  # 10.5281/zenodo.15718081
)

#' Core-trait list the portal computes means for.
#'
#' @param path the portal's `trait_groups_for_portal.csv`, recovered from
#'   portal commit dc8a965^ (see `data-raw/portal/README.md`).
#' @return character vector of trait names (96 for 7.0.0).
read_portal_core_traits <- function(
    path = "data-raw/portal/trait_groups_for_portal.csv") {
  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::filter(!is.na(core_trait)) |>
    dplyr::pull(trait)
}

#' Read an AusTraits release rds, checking it against the Zenodo checksum.
#'
#' @param path the release rds.
#' @param version release version, used to look up the expected MD5.
read_austraits_release <- function(path, version = "7.0.0") {
  expected <- austraits_release_md5[[version]]
  got <- unname(tools::md5sum(path))
  if (!identical(got, expected)) {
    stop("MD5 of ", path, " is ", got, "; the Zenodo release ", version,
         " has ", expected, ".")
  }
  readRDS(path)
}

#' Provenance written alongside the species means 7.0.0 product.
observed_means_portal_metadata <- function() {
  c(list(
    product = "AusTraits species means",
    version = "7.0.0",
    estimate_type = "observed_mean",
    method = "portal_arithmetic",
    derived_from = "AusTraits 7.0.0 (doi:10.5281/zenodo.15718081)",
    method_source = paste0("traitecoevo/austraits.portal@", portal_method_commit,
                           ", R/species_means.R"),
    trait_list = "the portal's 96 core traits (data-raw/portal/trait_groups_for_portal.csv)",
    documentation = "docs/observed-means-portal-method.md"
  ), repo_provenance())
}

#' Flattened core-trait records, as the portal fed them to the means.
#'
#' Reproduces the deployed build's preparation: the Bryant_2021 fix to
#' `methods$source_primary_key` applied in `data-raw/prepare_data.R`, then
#' `austraits::flatten_database()`, then the core-trait filter. The portal also
#' added a `row_id` and re-encoded `measurement_remarks`; neither column
#' reaches the means, so both steps are omitted.
#'
#' @param austraits a traits.build database object (the 7.0.0 release).
#' @param core_traits output of `read_portal_core_traits()`.
prepare_portal_records <- function(austraits, core_traits) {
  # A small fix for AusTraits v7.0.0: the source_primary_key for the
  # Bryant_2021_2 and Bryant_2021_3 datasets is incorrectly set.
  austraits$methods <- austraits$methods |>
    dplyr::mutate(
      source_primary_key = ifelse(grepl("Bryant_2021", dataset_id),
                                  dataset_id, source_primary_key)
    )

  austraits |>
    austraits::flatten_database() |>
    dplyr::filter(trait_name %in% core_traits)
}

#' Observed means by the portal method (port of
#' `estimate_species_trait_means()`).
#'
#' Numeric traits: records with value_type raw, mean or median become one mean
#' per dataset × location; minimum and maximum records from preserved
#' specimens or literature become one midpoint per observation; the species
#' value is the arithmetic mean of those. Categorical traits: each
#' space-separated value is counted, as `value (count)` in `value_count`.
#'
#' @param austraits flattened records (`prepare_portal_records()`).
#' @return one row per taxon × numeric trait × unit, and one or more rows per
#'   taxon × categorical trait, in the deployed portal's columns.
estimate_observed_means_portal <- function(austraits) {

  traits <-
    austraits |>
    dplyr::select(trait_name, unit) |>
    dplyr::distinct() |>
    dplyr::mutate(
      type = ifelse(is.na(unit), "categorical", "numerical")
    ) |>
    split(~type)

  data_means_numerical <-
    austraits |>
    .portal_means_numerical(traits = traits$numerical$trait_name)

  data_means_categorical <-
    austraits |>
    .portal_summary_categorical(traits = traits$categorical$trait_name)

  # Build species means for database
  austraits_species_averages <-
    data_means_numerical |>
    # combine the two datasets
    dplyr::bind_rows(data_means_categorical) |>
    # add in the trait type
    dplyr::left_join(by = c("trait_name"),
      traits |> dplyr::bind_rows() |> dplyr::select(trait_name, type),
    ) |>
    # add in taxon info
    dplyr::left_join(by = "taxon_name",
      austraits |>
        dplyr::select(taxon_name, taxon_rank:scientific_name_id,
                      -taxon_name_alternatives) |>
        dplyr::distinct()
    ) |>
    dplyr::relocate("value_count", .before = "value_mean") |>
    # Add row_id after aggregation
    dplyr::mutate(row_id = dplyr::row_number())

  # Reorder citations in source columns
  sources <- austraits |>
    dplyr::select(key = source_primary_key, source_primary_citation) |>
    dplyr::distinct()

  ids <- austraits_species_averages |>
    dplyr::select(source_primary_key) |>
    dplyr::distinct() |>
    dplyr::mutate(key = source_primary_key) |>
    tidyr::separate_rows(key, sep = "; ") |>
    dplyr::left_join(sources, by = c("key")) |>
    dplyr::group_by(source_primary_key) |>
    dplyr::summarise(.groups = "drop",
      source_primary_citation = paste(source_primary_citation, collapse = "; ")
    )

  austraits_species_averages |>
    dplyr::left_join(ids, by = "source_primary_key")
}

# Numeric traits: location means and flora midpoints, then the species mean.
.portal_means_numerical <- function(austraits, traits) {

  # any data that is a mean, median or raw, create a site mean
  location_means <- austraits |> .portal_means_locations(traits)

  # any data that is a max or a min (range) and basically from a flora, create
  # a mean value
  flora_means <- austraits |> .portal_means_floras(traits)

  # combine the two, then take means across site and flora replicates
  location_means |>
    dplyr::bind_rows(flora_means) |>
    dplyr::group_by(taxon_name, trait_name, unit) |>
    dplyr::summarise(.groups = "drop",
      value_mean = mean(value_mean),
      value_min = min(value_min),
      value_max = max(value_max),
      value_median = stats::median(value_median),
      # Quirk kept: `value_mean` is already the species mean here, so this is
      # 10^log10(value_mean), i.e. equal to value_mean.
      value_geom_mean = 10^mean(suppressWarnings(log10(value_mean)), na.rm = TRUE),
      all_replicates = sum(all_replicates),
      location_replicates = sum(location_replicates),
      flora_replicates = sum(flora_replicates),
      # record sources
      dataset_id = paste(unique(dataset_id), collapse = "; "),
      observation_id = paste(unique(observation_id), collapse = "; "),
      source_primary_key = paste(unique(stats::na.omit(source_primary_key)),
                                 collapse = "; ")
    ) |>
    dplyr::distinct()
}

# Records with value_type mean, raw or median -> one mean per dataset × location.
.portal_means_locations <- function(austraits, traits) {
  austraits |>
    dplyr::filter(
      trait_name %in% traits,
      value_type %in% c("mean", "raw", "median")
    ) |>
    dplyr::mutate(
      value = as.numeric(value),
      replicates = 1,
      log10_value = suppressWarnings(log10(value))
    ) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::group_by(taxon_name, trait_name, dataset_id, location_id, unit) |>
    dplyr::summarise(
      .groups = "drop",
      dplyr::across(value, list(mean = mean, min = min, max = max,
                                median = stats::median)),
      dplyr::across(c("latitude (deg)", "longitude (deg)", "location_name"),
                    dplyr::first),
      all_replicates = sum(replicates),
      value_geom_mean = 10^mean(log10_value, na.rm = TRUE),
      observation_id = paste(unique(observation_id), collapse = "; "),
      source_primary_key = paste(unique(stats::na.omit(source_primary_key)),
                                 collapse = "; ")
    ) |>
    dplyr::mutate(
      value_type = "location_mean",
      location_replicates = 1,
      flora_replicates = 0
    )
}

# Minimum/maximum records from preserved specimens or literature -> one
# midpoint per observation (a lone bound is its own midpoint).
.portal_means_floras <- function(austraits, traits) {
  flora_data <-
    austraits |>
    dplyr::filter(
      trait_name %in% traits,
      value_type %in% c("minimum", "maximum"),
      basis_of_record %in% c("preserved_specimen", "literature")
    )

  # early exit if no data available
  if (nrow(flora_data) == 0) return(flora_data)

  flora_data |>
    dplyr::mutate(
      value = as.numeric(value),
    ) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::group_by(taxon_name, trait_name, unit, dataset_id, observation_id,
                    original_name) |>
    dplyr::summarise(.groups = "keep",
      dplyr::across(value, list(mean = mean, min = min, max = max)),
      source_primary_key = paste(unique(stats::na.omit(source_primary_key)),
                                 collapse = "; ")
    ) |>
    dplyr::mutate(
      location_replicates = 0,
      flora_replicates = 1,
      all_replicates = 1,
      value_median = value_mean,
      value_type = "flora_mean"
    ) |>
    dplyr::ungroup()
}

# Categorical traits: `value (count)` per taxon × trait. One row is kept per
# distinct value's sources, so a taxon × trait can span several rows.
.portal_summary_categorical <- function(austraits, traits) {
  austraits |>
    .portal_value_summary_categorical(traits) |>
    dplyr::group_by(taxon_name, trait_name) |>
    dplyr::mutate(
      tmp_summary = paste0(value, " (", all_replicates, ")"),
      value_count = paste0(tmp_summary, collapse = "; ")
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-dplyr::all_of(c("tmp_summary", "value"))) |>
    dplyr::distinct()
}

.portal_value_summary_categorical <- function(austraits, traits) {
  austraits |>
    dplyr::filter(trait_name %in% traits) |>
    dplyr::select(dplyr::all_of(c("dataset_id", "taxon_name", "trait_name",
                                  "location_id", "observation_id", "value",
                                  "source_primary_key"))) |>
    dplyr::mutate(value = stringr::str_split(value, " ")) |>
    tidyr::unnest_longer(value) |>
    dplyr::mutate(
      all_replicates = 1
    ) |>
    dplyr::group_by(taxon_name, trait_name, value) |>
    dplyr::summarise(.groups = "drop",
      value = dplyr::first(value),
      all_replicates = sum(all_replicates),
      dataset_id = paste(unique(dataset_id), collapse = "; "),
      observation_id = paste(unique(observation_id), collapse = "; "),
      source_primary_key = paste(unique(stats::na.omit(source_primary_key)),
                                 collapse = "; ")
    ) |>
    dplyr::mutate(
      value_type = "value_count"
    )
}
