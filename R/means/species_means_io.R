# Writing the species means products: the data (parquet + gzipped csv), a
# column dictionary and a metadata file, one directory per product version
# (`outputs/species_means/<version>/`). Shared by the observed and modelled
# means, so every export carries `estimate_type` and `method` (terminology.md).

#' Label rows with how they were computed (terminology.md).
#'
#' @param data a species means table.
#' @param estimate_type `observed_mean` or `modelled_mean`.
#' @param method `geometric`, `arithmetic`, `rtmb`, `tmb` or `brms`.
#' @param after column after which the labels are placed.
label_estimates <- function(data, estimate_type, method, after = "unit") {
  stopifnot(estimate_type %in% c("observed_mean", "modelled_mean"),
            method %in% c("geometric", "arithmetic", "rtmb", "tmb", "brms"))
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
