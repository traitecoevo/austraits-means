# Observed means, geometric method (`estimate_type = observed_mean`,
# `method = geometric`; #15).
#
# The successor to the portal method (observed_means_portal.R), used from the
# next AusTraits build onward. It keeps the portal's two-stage recipe (one
# mean per dataset × location, one midpoint per flora observation, then the
# species value across those groups), but works on the log scale throughout,
# so the species value is a true geometric mean (terminology.md). Changes from
# the portal method:
#
#   - the species value (`value`) is a geometric mean: location means are
#     means of log values, a flora midpoint is the geometric midpoint
#     sqrt(minimum × maximum), and the species value is exp(mean of those);
#   - the arithmetic mean of the same groups is kept in `arithmetic_mean`;
#   - records with a zero or negative value are dropped (a geometric mean
#     needs positive values), as in the modelled means; traits whose values are
#     all negative, such as leaf_delta13C, drop out;
#   - per-dataset value-type rules (`dataset_value_type_rules()`) apply, as in
#     the modelled means: the AVH height dataset's specimen bounds are used as
#     observations, like raw records;
#   - categorical traits get one row per taxon × trait with the modal value;
#   - subspecies, varieties and forms are rolled up to species
#     (`roll_up_to_species()`), as for the modelled means, so the two join on
#     species.
#
# Kept as in the portal method, so the two stay comparable: no basis_of_record
# or life_stage filter on raw/mean/median records; minimum/maximum records
# from other datasets enter only from preserved specimens and the literature,
# and a lone bound is used as the midpoint; genus-level names are kept as
# their own rows. Each group counts once in the species value, whatever its
# number of records.

#' Traits in the geometric observed means: the portal's core traits, the nine
#' modelled traits, and the woodiness traits the community layers need.
geometric_mean_traits <- function(core_traits) {
  sort(unique(c(core_traits,
                "plant_height", "leaf_mass_per_area", "leaf_area",
                "leaf_N_per_dry_mass", "leaf_N_per_area", "leaf_P_per_area",
                "leaf_P_per_dry_mass", "seed_dry_mass", "wood_density",
                "woodiness", "woodiness_detailed", "plant_growth_form")))
}

#' Flattened records of the given traits from an austraits.build export,
#' with infraspecific taxa rolled up to species.
#'
#' Taxon details (`taxon_rank`, `establishment_means`, `taxon_id`) come from
#' the species' row of the build's taxa table; a rolled-up species absent from
#' that table gets rank "species" and keeps its records' genus and family.
#'
#' @param austraits_rds the export.
#' @param traits trait names to keep.
prepare_records <- function(austraits_rds, traits) {
  austraits <- readRDS(austraits_rds)
  taxa <- austraits$taxa

  austraits |>
    austraits::flatten_database() |>
    dplyr::filter(trait_name %in% traits) |>
    dplyr::select(
      dataset_id, taxon_name, observation_id, trait_name, value, unit,
      value_type, basis_of_record, life_stage, location_id, original_name,
      source_primary_key, source_primary_citation, taxon_rank, genus, family
    ) |>
    roll_up_to_species(dplyr::select(taxa, taxon_name, taxon_rank, binomial)) |>
    dplyr::mutate(taxon_rank = dplyr::if_else(taxon_name != taxon_name_aligned,
                                              "species", taxon_rank)) |>
    dplyr::left_join(
      dplyr::select(taxa, taxon_name, establishment_means, taxon_id),
      by = "taxon_name"
    )
}

#' Observed means by the geometric method.
#'
#' @param records output of `prepare_records()`.
#' @param rules per-dataset value-type rules (`dataset_value_type_rules()`).
#' @return one row per taxon × trait (× unit, for numeric traits).
estimate_observed_means_geometric <- function(records,
                                              rules = dataset_value_type_rules()) {
  numeric <- dplyr::filter(records, !is.na(unit))
  categorical <- dplyr::filter(records, is.na(unit))

  taxa <- records |>
    dplyr::distinct(taxon_name, taxon_rank, genus, family, establishment_means,
                    taxon_id) |>
    dplyr::distinct(taxon_name, .keep_all = TRUE)

  dplyr::bind_rows(
    .geometric_means_numeric(numeric, rules) |> dplyr::mutate(type = "numerical"),
    .geometric_summary_categorical(categorical) |> dplyr::mutate(type = "categorical")
  ) |>
    dplyr::left_join(taxa, by = "taxon_name") |>
    dplyr::arrange(taxon_name, trait_name)
}

.geometric_means_numeric <- function(records, rules) {
  records <- records |>
    dplyr::mutate(value = suppressWarnings(as.numeric(value))) |>
    dplyr::filter(!is.na(value), value > 0)

  ruled <- paste(rules$dataset_id, rules$value_type)
  is_ruled_dataset <- records$dataset_id %in% rules$dataset_id
  location_records <- records[
    (is_ruled_dataset & paste(records$dataset_id, records$value_type) %in% ruled) |
      (!is_ruled_dataset & records$value_type %in% c("raw", "mean", "median")), ]
  flora_records <- records[
    !is_ruled_dataset &
      records$value_type %in% c("minimum", "maximum") &
      records$basis_of_record %in% c("preserved_specimen", "literature"), ]

  group_summary <- function(d) {
    dplyr::summarise(d,
      .groups = "drop",
      log_mean = mean(log(value)),
      arithmetic = mean(value),
      value_min = min(value),
      value_max = max(value),
      n_records = dplyr::n(),
      source_primary_key = paste(unique(stats::na.omit(source_primary_key)),
                                 collapse = "; ")
    )
  }

  locations <- location_records |>
    dplyr::group_by(taxon_name, trait_name, unit, dataset_id, location_id) |>
    group_summary() |>
    dplyr::mutate(kind = "location")

  floras <- flora_records |>
    dplyr::group_by(taxon_name, trait_name, unit, dataset_id, observation_id,
                    original_name) |>
    group_summary() |>
    dplyr::mutate(kind = "flora")

  dplyr::bind_rows(locations, floras) |>
    dplyr::group_by(taxon_name, trait_name, unit) |>
    dplyr::summarise(
      .groups = "drop",
      log_value = mean(log_mean),
      arithmetic_mean = mean(arithmetic),
      value_min = min(value_min),
      value_max = max(value_max),
      n_records = sum(n_records),
      n_locations = sum(kind == "location"),
      n_flora = sum(kind == "flora"),
      n_datasets = dplyr::n_distinct(dataset_id),
      dataset_id = paste(sort(unique(dataset_id)), collapse = "; "),
      source_primary_key = paste(sort(unique(unlist(strsplit(
        source_primary_key[source_primary_key != ""], "; ", fixed = TRUE)))),
        collapse = "; ")
    ) |>
    dplyr::mutate(value = exp(log_value), .after = unit)
}

# Categorical traits: one row per taxon × trait, with value counts and the
# modal value (ties joined by a space, in sorted order).
.geometric_summary_categorical <- function(records) {
  counts <- records |>
    dplyr::mutate(value = stringr::str_split(value, " ")) |>
    tidyr::unnest_longer(value) |>
    dplyr::filter(!is.na(value), value != "") |>
    dplyr::count(taxon_name, trait_name, value, name = "n")

  modes <- counts |>
    dplyr::group_by(taxon_name, trait_name) |>
    dplyr::arrange(value, .by_group = TRUE) |>
    dplyr::summarise(
      .groups = "drop",
      value_mode = paste(value[n == max(n)], collapse = " "),
      value_count = paste0(value, " (", n, ")", collapse = "; ")
    )

  records |>
    dplyr::group_by(taxon_name, trait_name) |>
    dplyr::summarise(
      .groups = "drop",
      n_records = dplyr::n(),
      n_datasets = dplyr::n_distinct(dataset_id),
      dataset_id = paste(sort(unique(dataset_id)), collapse = "; "),
      source_primary_key = paste(sort(unique(stats::na.omit(source_primary_key))),
                                 collapse = "; ")
    ) |>
    dplyr::left_join(modes, by = c("taxon_name", "trait_name"))
}

#' Column dictionary for the geometric observed means.
observed_means_geometric_dictionary <- function() {
  tibble::tribble(
    ~column, ~definition,
    "taxon_name", "Taxon name as aligned in the AusTraits build, with subspecies, varieties and forms rolled up to their species. Genus-level names are kept as their own rows.",
    "trait_name", "AusTraits trait name (APD).",
    "unit", "Unit of the numeric trait values; empty for categorical traits.",
    "estimate_type", "Always observed_mean: computed directly from the records (terminology.md).",
    "method", "Always geometric: the geometric method (docs/observed-means-geometric-method.md).",
    "value", "Numeric traits: the observed mean, a geometric mean. Records are averaged on the log scale within each dataset × location, flora minimum and maximum bounds give a geometric midpoint per observation, and the species value is exp(mean) over those groups, each group counting once. Records with zero or negative values are not used.",
    "log_value", "Natural log of value.",
    "arithmetic_mean", "Arithmetic mean over the same groups and records (each group's arithmetic mean or midpoint, then their mean). Not the observed mean; given for comparison.",
    "value_min", "Smallest record value used.",
    "value_max", "Largest record value used.",
    "n_records", "Number of records used (each minimum and maximum bound counts as a record).",
    "n_locations", "Numeric traits: number of dataset × location groups.",
    "n_flora", "Numeric traits: number of flora observations (minimum/maximum bounds from preserved specimens or the literature).",
    "n_datasets", "Number of datasets contributing.",
    "dataset_id", "Datasets contributing, sorted, joined by '; '.",
    "source_primary_key", "Primary source keys of the contributing records, sorted, joined by '; '.",
    "value_mode", "Categorical traits: the most frequent value; tied values are joined by a space.",
    "value_count", "Categorical traits: each value with the number of records carrying it, as 'value (n)' joined by '; '. A record listing several space-separated values counts once for each.",
    "type", "numerical or categorical.",
    "taxon_rank", "Taxonomic rank of taxon_name.",
    "genus", "Genus.",
    "family", "Family.",
    "establishment_means", "Native, naturalised or both, as given by APC.",
    "taxon_id", "APC taxon identifier."
  )
}

#' Provenance written alongside the geometric observed means.
observed_means_geometric_metadata <- function(austraits_info, traits) {
  c(list(
    product = "AusTraits species means",
    version = austraits_info$austraits_version,
    estimate_type = "observed_mean",
    method = "geometric",
    derived_from = paste0("austraits.build development build ",
                          austraits_info$austraits_version, " (commit ",
                          austraits_info$austraits_build_commit, ")"),
    traits = traits,
    value_type_rules = "AVH_2026 minimum and maximum records used as observations",
    documentation = "docs/observed-means-geometric-method.md"
  ), repo_provenance())
}
