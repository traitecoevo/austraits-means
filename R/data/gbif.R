# GBIF occurrences: fetch a download by key, clean it, align its names to the
# pinned APC release and attach climate (#32).
#
# Ported from traitecoevo/cornsters_australian_annuals (R/download.R and
# R/clean_gbif.R at bca9cd8), with two additions (roadmap.md, decisions of 19
# October 2026):
#   - only openly licensed records (CC0 and CC BY) are kept, so products built
#     on the occurrences can be released under CC BY 4.0. GBIF licenses a
#     download as its most restrictive record, so the downloads themselves are
#     CC BY-NC; the cleaned extract is not;
#   - names are aligned to APC from the name the data provider supplied
#     (`verbatimScientificName`), falling back to GBIF's interpreted `species`
#     only where the provider name does not align to a species-level name.
#     GBIF interprets names with its own backbone, which since 2026 follows the
#     Catalogue of Life: it places some Australian species in genera APC does
#     not use (Corymbia maculata as Blakella maculata) and lumps some species
#     APC keeps apart, so its names do not line up with AusTraits'.
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
               "coordinateUncertaintyInMeters", "issue", "species",
               "verbatimScientificName", "taxonRank",
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

#' Align GBIF names to the pinned APC release: the provider's name first,
#' GBIF's species name as fallback.
#'
#' Both name sets are aligned with APCalign. For each provider name × GBIF
#' species pair, the provider name's alignment is used if it reaches species
#' rank or below, otherwise the GBIF name's. Infraspecific names are rolled up
#' to their species (`taxon_name`), as on the trait side, keeping the
#' infraspecific name in `taxon_name_infraspecific`.
#'
#' A fuzzy match is accepted only if it keeps the genus, allowing a
#' one-character difference for spelling variants (Leschenaultia /
#' Lechenaultia); fuzzy matches that change the genus are mostly wrong
#' (Acanthus to Amaranthus, Setaria to Scleria) and count as no match.
#'
#' @param cleaned output of `clean_gbif()`.
#' @param apc_version APCalign release tag (as `apc_version` in `_targets.R`).
#' @param lookup optional APCalign lookup from an earlier run (the `lookup`
#'   attribute), reused if it covers every name, to skip the slow alignment.
#' @return one row per provider name × GBIF species pair: `verbatim_name`,
#'   `gbif_name`, `taxon_name`, `taxon_name_infraspecific`, `genus`,
#'   `name_source` (provider, gbif or none), `aligned_reason`, `apc_version`;
#'   plus the full APCalign lookup in `attr(, "lookup")`.
align_gbif_names <- function(cleaned, apc_version, lookup = NULL) {
  pairs <- dplyr::distinct(tibble::as_tibble(cleaned),
                           verbatim_name = verbatimScientificName, gbif_name = species)
  names_in <- sort(unique(stats::na.omit(c(pairs$verbatim_name, pairs$gbif_name))))
  names_in <- names_in[names_in != ""]
  if (is.null(lookup) || !all(names_in %in% lookup$original_name)) {
    resources <- APCalign::load_taxonomic_resources(version = apc_version, quiet = TRUE)
    lookup <- APCalign::create_taxonomic_update_lookup(
      names_in, resources = resources,
      taxonomic_splits = "most_likely_species", quiet = TRUE
    )
  }
  species_level <- c("species", "subspecies", "variety", "form")
  first_word <- function(x) sub(" .*", "", x)
  genus_kept <- function(from, to) {
    diag(utils::adist(first_word(from), first_word(to))) <= 1
  }
  hit <- lookup |>
    dplyr::filter(taxon_rank %in% species_level, !is.na(accepted_name)) |>
    dplyr::filter(!grepl("fuzzy", aligned_reason, ignore.case = TRUE) |
                    genus_kept(original_name, aligned_name)) |>
    dplyr::select(original_name, accepted_name, taxon_rank, genus, aligned_reason)

  out <- pairs |>
    dplyr::left_join(dplyr::rename_with(hit, ~ paste0(.x, "_v"), -original_name),
                     by = c(verbatim_name = "original_name")) |>
    dplyr::left_join(dplyr::rename_with(hit, ~ paste0(.x, "_g"), -original_name),
                     by = c(gbif_name = "original_name")) |>
    dplyr::mutate(
      name_source = dplyr::case_when(!is.na(accepted_name_v) ~ "provider",
                                     !is.na(accepted_name_g) ~ "gbif",
                                     TRUE ~ "none"),
      taxon_name_infraspecific = dplyr::coalesce(accepted_name_v, accepted_name_g),
      rank = dplyr::coalesce(taxon_rank_v, taxon_rank_g),
      genus = dplyr::coalesce(genus_v, genus_g),
      aligned_reason = dplyr::coalesce(aligned_reason_v, aligned_reason_g),
      # roll infraspecific names up to their species (first two words; APC
      # infraspecific names are "Genus epithet rank name")
      taxon_name = dplyr::if_else(
        rank %in% c("subspecies", "variety", "form"),
        vapply(strsplit(dplyr::coalesce(taxon_name_infraspecific, ""), " ", fixed = TRUE),
               \(w) paste(utils::head(w, 2), collapse = " "), character(1)),
        taxon_name_infraspecific),
      apc_version = apc_version
    ) |>
    dplyr::select(verbatim_name, gbif_name, taxon_name, taxon_name_infraspecific,
                  genus, name_source, aligned_reason, apc_version)
  attr(out, "lookup") <- lookup
  out
}

#' Attach climate to cleaned GBIF records.
#'
#' `taxon_name` is the APC accepted species from `align_gbif_names()`. Records
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

  tibble::as_tibble(cleaned) |>
    dplyr::transmute(
      original_name = verbatimScientificName, gbif_name = species,
      decimal_longitude = decimalLongitude, decimal_latitude = decimalLatitude,
      year, gbif_id = gbifID, dataset_key = datasetKey, license,
      basis_of_record = basisOfRecord, institution_code = institutionCode,
      coordinate_uncertainty_in_meters = coordinateUncertaintyInMeters,
      gbif_taxon_rank = taxonRank
    ) |>
    dplyr::mutate(annual_precipitation = climate$wc2.1_2.5m_bio_12,
                  average_temperature = climate$wc2.1_2.5m_bio_1) |>
    dplyr::left_join(
      dplyr::select(aligned, verbatim_name, gbif_name, taxon_name,
                    taxon_name_infraspecific, genus, name_source),
      by = c(original_name = "verbatim_name", "gbif_name")) |>
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
    dplyr::relocate(original_name, gbif_name, taxon_name, genus, name_source)
}
