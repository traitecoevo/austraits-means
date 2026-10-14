# Reliability indicators for the species means (#16): how much evidence sits
# behind each species × trait value, added to all three products.
#
# `reliability` is a PROVISIONAL class; its thresholds are documented as such:
#   high   - at least 3 datasets and at least 10 records;
#   medium - at least 2 datasets or at least 5 records;
#   low    - otherwise.

#' Provisional reliability class from dataset and record counts.
reliability_class <- function(n_datasets, n_records) {
  dplyr::case_when(
    n_datasets >= 3 & n_records >= 10 ~ "high",
    n_datasets >= 2 | n_records >= 5 ~ "medium",
    TRUE ~ "low"
  )
}

#' Column definitions shared by every product's dictionary.
reliability_dictionary <- function() {
  tibble::tribble(
    ~column, ~definition,
    "n_records", "Number of records used.",
    "n_datasets", "Number of AusTraits datasets contributing records.",
    "n_locations", "Numeric traits: number of distinct dataset × location groups (modelled means: distinct record coordinates).",
    "n_flora", "Numeric traits: number of minimum/maximum records from preserved specimens or the literature (floras) used.",
    "n_field", "Number of records used whose basis of record is field.",
    "sd_log", "Numeric traits: standard deviation of the natural log of the record values used; empty with a single record.",
    "reliability", "Provisional class from the evidence behind the value: high (at least 3 datasets and 10 records), medium (at least 2 datasets or 5 records), low (otherwise). Thresholds are provisional."
  )
}

#' Reliability indicators for the portal-method product, from the records the
#' portal method used.
#'
#' Numeric rows: records with value_type raw, mean or median (any basis),
#' plus minimum/maximum records from preserved specimens or the literature,
#' that parse as numbers. Categorical rows: every record of the trait.
#'
#' @param records output of `prepare_portal_records()`.
#' @return list(numeric = per taxon × trait × unit, categorical = per taxon ×
#'   trait).
portal_reliability <- function(records) {
  numeric <- records |>
    dplyr::filter(!is.na(unit)) |>
    dplyr::mutate(value = suppressWarnings(as.numeric(value)),
                  flora = value_type %in% c("minimum", "maximum") &
                    basis_of_record %in% c("preserved_specimen", "literature")) |>
    dplyr::filter(!is.na(value),
                  value_type %in% c("mean", "raw", "median") | flora) |>
    dplyr::group_by(taxon_name, trait_name, unit) |>
    dplyr::summarise(
      .groups = "drop",
      n_records = dplyr::n(),
      n_datasets = dplyr::n_distinct(dataset_id),
      n_locations = dplyr::n_distinct(paste(dataset_id, location_id)[!flora]),
      n_flora = sum(flora),
      n_field = sum(basis_of_record %in% "field"),
      sd_log = if (sum(value > 0) > 1) stats::sd(log(value[value > 0])) else NA_real_
    )
  categorical <- records |>
    dplyr::filter(is.na(unit)) |>
    dplyr::group_by(taxon_name, trait_name) |>
    dplyr::summarise(
      .groups = "drop",
      n_records = dplyr::n(),
      n_datasets = dplyr::n_distinct(dataset_id),
      n_field = sum(basis_of_record %in% "field")
    )
  list(numeric = numeric, categorical = categorical)
}

#' Add reliability indicators to the portal-method product, after its own
#' columns (which stay exactly as the portal deployed them).
add_portal_reliability <- function(means, reliability) {
  num <- means |>
    dplyr::filter(type == "numerical") |>
    dplyr::left_join(reliability$numeric, by = c("taxon_name", "trait_name", "unit"))
  cat <- means |>
    dplyr::filter(type == "categorical") |>
    dplyr::left_join(reliability$categorical, by = c("taxon_name", "trait_name"))
  dplyr::bind_rows(num, cat) |>
    dplyr::arrange(row_id) |>
    dplyr::mutate(reliability = reliability_class(n_datasets, n_records))
}

#' Reliability indicators for modelled means, from the records in the fit.
#'
#' @param data a fit's `data` (field records of species in the model).
modelled_reliability <- function(data) {
  data |>
    dplyr::group_by(taxon_name) |>
    dplyr::summarise(
      .groups = "drop",
      n_datasets = dplyr::n_distinct(dataset_id),
      n_locations = dplyr::n_distinct(paste(decimal_longitude, decimal_latitude)),
      n_flora = 0L,
      n_field = dplyr::n(),
      sd_log = stats::sd(log_value)
    )
}
