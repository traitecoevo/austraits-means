# GBIF occurrences: fetch a download by key, clean it, align its names to the
# pinned APC release and attach climate (#32).
#
# Ported from traitecoevo/cornsters_australian_annuals (R/download.R and
# R/clean_gbif.R at bca9cd8), with one addition: only openly licensed records
# (CC0 and CC BY) are kept, so products built on the occurrences can be
# released under CC BY 4.0 (roadmap.md, decision of 19 October 2026). GBIF
# licenses a download as its most restrictive record, so the downloads
# themselves are CC BY-NC; the cleaned extract is not.
#
# The cleaned extract is a Tier B input: build it deliberately with
# R/scripts/build_gbif.R, not in tar_make().

#' GBIF downloads used by this project.
gbif_downloads <- tibble::tribble(
  ~key,                      ~doi,                 ~created,     ~records,
  "0001494-250914085247600", "10.15468/dl.b672m4", "2025-09-14", 24385439,
  "0027765-260928105237408", "10.15468/dl.766erw", "2026-10-15", 27357959
)

#' Licences kept when cleaning (the GBIF `license` column).
gbif_open_licences <- c("CC0_1_0", "CC_BY_4_0")

#' Citation of a GBIF download, as GBIF asks for it.
gbif_citation <- function(key) {
  d <- gbif_downloads[gbif_downloads$key == key, ]
  stopifnot(nrow(d) == 1)
  sprintf("GBIF.org (%s) GBIF Occurrence Download https://doi.org/%s",
          format(as.Date(d$created), "%d %B %Y"), d$doi)
}

#' Path to a download's raw occurrence table, fetching it if absent.
#'
#' An executed download stays retrievable by key without GBIF credentials.
#' The table is GBIF's SIMPLE_CSV format, which is tab-delimited despite its
#' extension (~4 GB zipped, ~17 GB unzipped).
#'
#' @param key GBIF download key.
#' @param dest_dir folder holding (or to hold) `<key>.zip` and `<key>.csv`.
download_gbif <- function(key, dest_dir) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  csv <- file.path(dest_dir, paste0(key, ".csv"))
  if (file.exists(csv)) return(csv)

  zip <- file.path(dest_dir, paste0(key, ".zip"))
  if (!file.exists(zip)) {
    url <- paste0("https://api.gbif.org/v1/occurrence/download/request/", key, ".zip")
    message("Downloading GBIF ", key, " (~4 GB); slow, but only once.")
    withr::with_options(list(timeout = 60 * 60),
                        utils::download.file(url, zip, mode = "wb"))
  }
  utils::unzip(zip, exdir = dest_dir)
  if (!file.exists(csv)) {
    stop("GBIF zip did not contain ", basename(csv), ". Contents: ",
         paste(utils::unzip(zip, list = TRUE)$Name, collapse = ", "))
  }
  csv
}

# Reference tables for CoordinateCleaner: AU country/state centroids and
# capitals (with capitals missing from countryref added), and institutions
# (with CSIRO-Atherton corrected).
gbif_ref_capitals <- function() {
  CoordinateCleaner::countryref |>
    dplyr::filter(iso2 == "AU") |>
    dplyr::bind_rows(
      tibble::tibble(
        capital.lat = c(-33.8688, -35.2809, -27.4698, -31.9505, -12.4634, -34.9285, -37.8136, -27.0000),
        capital.lon = c(151.2093, 149.1300, 153.0251, 115.8605, 130.8456, 138.6007, 144.9631, 151.0000),
        capital = c("Sydney", "Canberra", "Brisbane", "Melbourne", "Darwin", "Adelaide", "Hobart", "Perth"),
        name = c("New South Wales", "Australian Capital Territory", "Queensland", "Victoria",
                 "Northern Territory", "South Australia", "Western Australia", "Tasmania")
      ) |>
        dplyr::mutate(iso3 = "AUS", iso2 = "AU", type = "province")
    )
}

gbif_ref_institutions <- function() {
  CoordinateCleaner::institutions |>
    dplyr::filter(country == "AUS") |>
    dplyr::bind_rows(tibble::tibble(decimalLatitude = -17.2586301,
                                    decimalLongitude = 145.482298,
                                    name = "CSIRO-Atherton"))
}

#' Record-level filters on a raw GBIF table: licence, provenance, coordinate
#' quality, date and rank. Pure and fast; tested in
#' tests/testthat/test-gbif.R.
#'
#' @param x raw SIMPLE_CSV columns (GBIF names, e.g. `decimalLongitude`).
#' @param licences licences to keep.
filter_gbif_records <- function(x, licences = gbif_open_licences) {
  x |>
    dplyr::filter(license %in% licences) |>
    dplyr::filter(countryCode == "AU") |> # already filtered on the GBIF website
    dplyr::filter(occurrenceStatus == "PRESENT") |>
    dplyr::filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) |>
    dplyr::filter(!basisOfRecord %in% c("LIVING_SPECIMEN", "FOSSIL_SPECIMEN")) |>
    dplyr::filter(coordinatePrecision < 0.05 | is.na(coordinatePrecision)) |> # retains most herbarium records
    dplyr::filter(year >= 1900) |> # pre-systematic-botany records have unreliable locations
    dplyr::filter(coordinateUncertaintyInMeters < 10000 | is.na(coordinateUncertaintyInMeters)) |> # ~township scale
    dplyr::filter(!(decimalLatitude == 0 | decimalLongitude == 0)) |>
    dplyr::filter(!grepl("COUNTRY_COORDINATE_MISMATCH", issue),
                  !grepl("RECORDED_DATE_UNLIKELY", issue)) |>
    # subspecies and varieties kept for breadth; names are rolled up later
    dplyr::filter(taxonRank %in% c("SPECIES", "SUBSPECIES", "VARIETY"),
                  !is.na(species), species != "") |>
    dplyr::filter(
      !institutionCode %in% "NHMUK", # digitised herbarium sheets proved unreliable
      !institutionCode %in% "CJBG"   # a small collection with many georeferencing errors
    )
}

#' The CoordinateCleaner passes, applied in row blocks.
#'
#' `cc_cen()` and `cc_cap()` buffer every input point to get a bounding box,
#' so their memory grows with the number of rows; blocks of 1 million rows
#' cut peak memory by more than half and give identical results, because each
#' pass tests every record against a fixed reference set (cornsters
#' R/clean_gbif.R has the measurements; ropensci/CoordinateCleaner#113).
clean_gbif_coordinates <- function(x, chunk_rows = 1e6) {
  ref_capitals <- gbif_ref_capitals()
  ref_institutions <- gbif_ref_institutions()
  run_passes <- function(d) {
    d |>
      CoordinateCleaner::cc_val(verbose = FALSE) |>
      CoordinateCleaner::cc_equ(verbose = FALSE) |>
      CoordinateCleaner::cc_gbif(verbose = FALSE) |>
      CoordinateCleaner::cc_cen(buffer = 5000, ref = ref_capitals, verbose = FALSE) |> # country/state centroids
      CoordinateCleaner::cc_cap(buffer = 2000, ref = ref_capitals, verbose = FALSE) |> # capital-city default coordinates, keeping genuine urban records
      CoordinateCleaner::cc_inst(buffer = 2000, ref = ref_institutions, verbose = FALSE) # herbarium and zoo accessions
  }
  if (is.null(chunk_rows) || nrow(x) <= chunk_rows) return(run_passes(x))
  blocks <- split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / chunk_rows))
  out <- lapply(blocks, function(rows) {
    r <- run_passes(x[rows, ])
    gc()
    r
  })
  data.table::rbindlist(out)
}

#' Clean a raw GBIF download.
#'
#' @param gbif_csv the raw SIMPLE_CSV table (`download_gbif()`).
#' @param licences licences to keep.
#' @param chunk_rows block size for the CoordinateCleaner passes.
#' @return the cleaned records (GBIF column names), with the record count
#'   after each step in `attr(, "steps")`.
clean_gbif <- function(gbif_csv, licences = gbif_open_licences, chunk_rows = 1e6) {
  columns <- c("gbifID", "datasetKey", "occurrenceStatus", "decimalLongitude",
               "decimalLatitude", "basisOfRecord", "coordinatePrecision",
               "coordinateUncertaintyInMeters", "issue", "species", "taxonRank",
               "year", "countryCode", "institutionCode", "license")
  # only data.table::fread reads SIMPLE_CSV reliably (readr and arrow error)
  x <- data.table::fread(gbif_csv, quote = "", select = columns, showProgress = FALSE)
  steps <- c(raw = nrow(x))
  x <- filter_gbif_records(x, licences)
  steps["record_filters"] <- nrow(x)
  x <- clean_gbif_coordinates(x, chunk_rows)
  steps["coordinate_cleaner"] <- nrow(x)
  # one record per location × species × year: keeps resurveys in other years,
  # drops within-year replicates from the same spot
  x <- dplyr::distinct(x, decimalLongitude, decimalLatitude, species, year,
                       .keep_all = TRUE)
  steps["deduplicated"] <- nrow(x)
  x <- dplyr::select(x, -occurrenceStatus, -countryCode, -coordinatePrecision)
  attr(x, "steps") <- steps
  x
}

#' Align GBIF species names to the pinned APC release.
#'
#' @param species unique GBIF `species` names.
#' @param apc_version APCalign release tag (as `apc_version` in `_targets.R`).
align_gbif_names <- function(species, apc_version) {
  resources <- APCalign::load_taxonomic_resources(version = apc_version, quiet = TRUE)
  APCalign::create_taxonomic_update_lookup(
    sort(unique(stats::na.omit(species))), resources = resources,
    # GBIF has already resolved much of the taxonomy; keep its records
    taxonomic_splits = "most_likely_species", quiet = TRUE
  ) |>
    dplyr::mutate(apc_version = apc_version)
}

#' Attach climate to cleaned GBIF records.
#'
#' `taxon_name` is the APC accepted name from `align_gbif_names()`. Records
#' without climate, with non-positive precipitation, or on the known rogue
#' point at the continental centroid (134.1, -25.7) are dropped.
#'
#' @param cleaned output of `clean_gbif()`.
#' @param aligned output of `align_gbif_names()`.
#' @param rasters climate rasters (`load_climate_rasters()`).
build_gbif_climate <- function(cleaned, aligned, rasters = load_climate_rasters()) {
  climate <- terra::extract(
    rasters,
    terra::vect(as.data.frame(cleaned[, c("decimalLongitude", "decimalLatitude")]),
                geom = c("decimalLongitude", "decimalLatitude"),
                crs = terra::crs(rasters)),
    method = "simple", ID = FALSE
  )
  names_lookup <- aligned |>
    dplyr::select(original_name, taxon_name = accepted_name, genus) |>
    dplyr::distinct(original_name, .keep_all = TRUE)

  tibble::as_tibble(cleaned) |>
    dplyr::transmute(
      original_name = species,
      decimal_longitude = decimalLongitude, decimal_latitude = decimalLatitude,
      year, gbif_id = gbifID, dataset_key = datasetKey, license,
      basis_of_record = basisOfRecord, institution_code = institutionCode,
      coordinate_uncertainty_in_meters = coordinateUncertaintyInMeters,
      taxon_rank = taxonRank
    ) |>
    dplyr::mutate(annual_precipitation = climate$wc2.1_2.5m_bio_12,
                  average_temperature = climate$wc2.1_2.5m_bio_1) |>
    dplyr::left_join(names_lookup, by = "original_name") |>
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
        (average_temperature - mean(average_temperature)) / stats::sd(average_temperature)
    ) |>
    dplyr::relocate(original_name, taxon_name, genus)
}
