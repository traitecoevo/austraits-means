# Writing the species means products: the data (parquet + gzipped csv), a
# column dictionary and a metadata file, one directory per product version
# (`outputs/species_means/<version>/`). Shared by the observed and modelled
# means, so every export carries `estimate_type` and `method` (terminology.md).

#' Label rows with how they were computed (terminology.md).
#'
#' @param data a species means table.
#' @param estimate_type `observed_mean` or `modelled_mean`.
#' @param method `portal_arithmetic`, `geometric`, `tmb` or `brms`.
#' @param after column after which the labels are placed.
label_estimates <- function(data, estimate_type, method, after = "unit") {
  stopifnot(estimate_type %in% c("observed_mean", "modelled_mean"),
            method %in% c("portal_arithmetic", "geometric", "tmb", "brms"))
  data |>
    dplyr::mutate(estimate_type = estimate_type, method = method) |>
    dplyr::relocate(estimate_type, method, .after = dplyr::all_of(after))
}

#' Write one species means product.
#'
#' Writes `<stem>.parquet`, `<stem>.csv.gz`, `<stem>_dictionary.csv` and
#' `<stem>_metadata.json` to `dir`. The dictionary must define every column.
#'
#' @param data the product table.
#' @param dictionary tibble(column, definition), one row per column of `data`.
#' @param dir output directory.
#' @param stem file name stem.
#' @param metadata named list written as JSON (provenance of the product).
#' @return the four file paths.
write_species_means <- function(data, dictionary, dir, stem, metadata = list()) {
  undefined <- setdiff(names(data), dictionary$column)
  if (length(undefined)) {
    stop("Columns missing from the dictionary: ", paste(undefined, collapse = ", "))
  }
  dictionary <- dictionary |>
    dplyr::filter(column %in% names(data)) |>
    dplyr::mutate(
      type = vapply(column, \(x) class(data[[x]])[1], character(1)),
      .after = column
    ) |>
    dplyr::arrange(match(column, names(data)))

  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, paste0(stem, c(".parquet", ".csv.gz",
                                         "_dictionary.csv", "_metadata.json")))
  arrow::write_parquet(data, paths[1])
  readr::write_csv(data, paths[2], na = "")
  readr::write_csv(dictionary, paths[3])
  metadata <- c(metadata, list(n_rows = nrow(data),
                               written = format(Sys.time(), "%Y-%m-%d")))
  jsonlite::write_json(metadata, paths[4], auto_unbox = TRUE, pretty = TRUE)
  paths
}

#' Provenance of this repository's checkout, for product metadata.
repo_provenance <- function() {
  sha <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                  error = function(e) NA_character_, warning = function(w) NA_character_)
  list(repository = "https://github.com/traitecoevo/austraits-means",
       commit = sha)
}

#' Column dictionary for the portal-method observed means.
observed_means_portal_dictionary <- function() {
  tibble::tribble(
    ~column, ~definition,
    "taxon_name", "Taxon name as aligned in the AusTraits release (APC-accepted where possible). Infraspecific and genus-level names are kept as they are, not rolled up to species.",
    "trait_name", "AusTraits trait name (APD).",
    "unit", "Unit of the numeric trait values; empty for categorical traits.",
    "estimate_type", "Always observed_mean: computed directly from the records (terminology.md).",
    "method", "Always portal_arithmetic: the AusTraits data portal's method for release 7.0.0, an arithmetic mean (see the workflow document).",
    "value_count", "Categorical traits only: each value with the number of records carrying it, as 'value (n)' joined by '; '. A record listing several space-separated values counts once for each.",
    "value_mean", "Numeric traits only: arithmetic mean, over the taxon's location means and flora midpoints, of those values. Each location mean is the arithmetic mean of the raw, mean and median records for one dataset at one location; each flora midpoint is the mean of the minimum and maximum records of one observation from a preserved specimen or the literature.",
    "value_min", "Numeric traits only: smallest record value across the location records and flora bounds used.",
    "value_max", "Numeric traits only: largest record value across the location records and flora bounds used.",
    "value_median", "Numeric traits only: median of the location medians and flora midpoints (a median of group medians, not the median of all records).",
    "value_geom_mean", "Numeric traits only: equal to value_mean where that is positive, and NaN where it is negative (leaf_delta13C). Despite its name this is not a geometric mean; the portal method computes it from the already-averaged value. Kept for fidelity to the deployed portal output.",
    "all_replicates", "Number of records used: each raw, mean or median record counts once and each flora observation counts once (whether it has one bound or two). For categorical traits, the number of records carrying the values on this row.",
    "location_replicates", "Numeric traits only: number of dataset × location groups contributing a location mean.",
    "flora_replicates", "Numeric traits only: number of flora observations contributing a midpoint.",
    "dataset_id", "AusTraits datasets contributing, joined by '; '.",
    "observation_id", "AusTraits observation identifiers contributing, joined by '; '. Identifiers are unique only within a dataset.",
    "source_primary_key", "Primary source keys of the contributing datasets, joined by '; '.",
    "value_type", "value_count for categorical rows; empty for numeric rows.",
    "type", "numerical or categorical.",
    "taxon_rank", "Taxonomic rank of taxon_name, from the AusTraits taxa table.",
    "taxonomic_status", "Taxonomic status of taxon_name in the taxonomic dataset.",
    "taxonomic_dataset", "Taxonomic reference the name was aligned to (APC or APNI).",
    "genus", "Genus, from the AusTraits taxa table.",
    "family", "Family, from the AusTraits taxa table.",
    "binomial", "Binomial of taxon_name, from the AusTraits taxa table.",
    "trinomial", "Trinomial of taxon_name, where applicable.",
    "taxon_distribution", "Distribution by state or territory as given by APC.",
    "establishment_means", "Native, naturalised or both, as given by APC.",
    "scientific_name", "Scientific name with authority.",
    "taxon_id", "APC taxon identifier.",
    "taxon_id_genus", "APC identifier of the genus.",
    "taxon_id_family", "APC identifier of the family.",
    "scientific_name_id", "APNI scientific name identifier.",
    "row_id", "Row number in this table. It depends on row order and is not a stable identifier.",
    "source_primary_citation", "Citations for source_primary_key, joined by '; '."
  )
}
