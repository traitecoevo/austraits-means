# Builders for the upstream source artefacts.
#
# The GBIF occurrences are expensive and rarely rebuilt, so in the targets
# pipeline (`_targets.R`) the cleaned extract is tracked as a *file input*
# rather than recomputed on every `tar_make()`. Its builders (download,
# cleaning, name alignment, climate) are in R/data/gbif.R and run by
# R/scripts/build_gbif.R. Medoids are built by `climate_centres.qmd`
# (process_taxon_medoid over all taxa). The name crosswalk
# (`build_name_crosswalk()`) runs inside the pipeline, with the APC release
# pinned in `_targets.R`.
#
# Trait records, height included, come from an austraits.build export
# (`data-raw/austraits/README.md`); `build_austraits_climate()` runs inside the
# pipeline, because flattening a build takes seconds.

#' Name crosswalk from trait-record names to GBIF names, via one APC release.
#'
#' The trait records (austraits.build) and the GBIF occurrences, and so the
#' medoids, were name-aligned against different APC releases, so a species can
#' carry different names in the two (#1). Both name sets are aligned to the
#' same pinned APC release with APCalign and matched on the accepted name.
#' Names that do not align to an accepted name are left out; such names can
#' still match a GBIF name directly in `build_medoid_lookup()`.
#'
#' @param taxon_names trait-record taxon names (all taxa in the build).
#' @param gbif_names taxon names in the GBIF occurrences (the medoid keys).
#' @param apc_version APC/APNI release to align against (an APCalign release
#'   tag, e.g. "2026-07-30").
#' @return tibble(taxon_name, accepted_name, gbif_taxon_name, apc_version); a
#'   trait name can map to more than one GBIF name (synonyms in GBIF).
build_name_crosswalk <- function(taxon_names, gbif_names, apc_version) {
  resources <- APCalign::load_taxonomic_resources(version = apc_version,
                                                  quiet = TRUE)
  align_accepted <- function(names) {
    APCalign::create_taxonomic_update_lookup(
      sort(unique(stats::na.omit(names))), resources = resources,
      quiet = TRUE
    ) |>
      dplyr::select(original_name, accepted_name = suggested_name) |>
      dplyr::filter(!is.na(accepted_name))
  }

  trait <- align_accepted(taxon_names) |> dplyr::rename(taxon_name = original_name)
  gbif <- align_accepted(gbif_names) |> dplyr::rename(gbif_taxon_name = original_name)

  trait |>
    dplyr::inner_join(gbif, by = "accepted_name", relationship = "many-to-many") |>
    dplyr::mutate(apc_version = apc_version) |>
    dplyr::arrange(taxon_name, gbif_taxon_name)
}

#' Taxon names in the GBIF occurrence parquet.
gbif_taxon_names <- function(gbif_path) {
  arrow::open_dataset(gbif_path) |>
    dplyr::distinct(taxon_name) |>
    dplyr::collect() |>
    dplyr::pull(taxon_name) |>
    stats::na.omit() |>
    sort()
}

#' GBIF download the cleaned occurrence extract in use was made from (a key
#' in `gbif_downloads`, R/data/gbif.R).
gbif_download_key <- "0027765-260928105237408"

#' Data basis of the modelled means: where each input came from and when.
#'
#' @param austraits_info output of `austraits_build_info()`.
#' @param apc_version APC release of the name crosswalk.
#' @param gbif_file,medoids_file the cleaned GBIF parquet and the medoids.
#' @return named list, written into product metadata.
data_basis <- function(austraits_info, apc_version, gbif_file, medoids_file) {
  c(austraits_info, list(
    apc_version = apc_version,
    apcalign_package = as.character(utils::packageVersion("APCalign")),
    gbif_download = gbif_citation(gbif_download_key),
    gbif_extract_built = format(file.mtime(gbif_file), "%Y-%m-%d"),
    medoids_built = format(file.mtime(medoids_file), "%Y-%m-%d")
  ))
}

#' Version string and commit of an austraits.build export.
#'
#' @param austraits_rds an rds written by austraits.build (`build.R`).
#' @return list(austraits_version, austraits_build_commit).
austraits_build_info <- function(austraits_rds) {
  info <- readRDS(austraits_rds)$build_info
  list(austraits_version = info$version,
       austraits_build_commit = info$git_SHA)
}

#' Flatten AusTraits and attach climate for the analysis traits.
#'
#' Reads an austraits.build export (a development build since #12; it carries
#' the AVH height dataset `AVH_2026`), keeps field records of the analysis
#' traits that have coordinates, and attaches WorldClim climate. Value-type
#' rules are applied later, in `load_austraits_records()`.
#'
#' @param austraits_rds the austraits.build export.
#' @param traits trait names to keep.
#' @param rasters climate rasters (`load_climate_rasters()`).
#' @return tibble of field, georeferenced trait records with climate attached.
build_austraits_climate <- function(
    austraits_rds,
    traits = c("plant_height", "leaf_mass_per_area", "leaf_area",
               "leaf_N_per_dry_mass", "leaf_N_per_area", "leaf_P_per_area",
               "leaf_P_per_dry_mass", "seed_dry_mass", "wood_density"),
    rasters = load_climate_rasters()) {

  austraits <- readRDS(austraits_rds) |>
    austraits::flatten_database()

  austraits |>
    dplyr::filter(
      trait_name %in% traits,
      basis_of_record == "field",
      !is.na(`latitude (deg)`), !is.na(`longitude (deg)`)
    ) |>
    dplyr::select(
      dataset_id, taxon_name, observation_id, trait_name, value, unit,
      value_type, basis_of_record, life_stage, location_id, collection_date,
      original_name, taxon_rank, genus, family,
      `latitude (deg)`, `longitude (deg)`
    ) |>
    dplyr::mutate(dplyr::across(
      c("latitude (deg)", "longitude (deg)", value),
      ~ suppressWarnings(as.numeric(.x))
    )) |>
    dplyr::filter(!is.na(`latitude (deg)`), !is.na(`longitude (deg)`),
                  !is.na(value)) |>
    dplyr::rename(decimal_longitude = `longitude (deg)`,
                  decimal_latitude = `latitude (deg)`) |>
    extract_climate(rasters = rasters) |>
    dplyr::filter(!is.na(annual_precipitation), !is.na(average_temperature))
}

#' Climate medoid of every taxon in a GBIF extract (#33).
#'
#' Ported from `climate_centres.qmd` / `process_taxon_medoid()`, without the
#' per-taxon plots. For each taxon: distinct points, a random subsample of
#' `max_n` points for very common taxa (the distance matrix is O(n^2)), then
#' `locate_medoid_2d()`, which drops the points furthest from all others in
#' standardised climate space and returns the point closest to the rest.
#' Taxa run in parallel (furrr, multisession) in blocks balanced by size.
#'
#' @param gbif_path a cleaned extract (`R/scripts/build_gbif.R`).
#' @param workers number of parallel workers.
#' @param max_n largest number of points per taxon.
#' @param cull_percentage share of outlying points dropped (at most 30).
#' @param seed subsampling seed, set per taxon so results do not depend on
#'   which worker runs it.
#' @return one row per taxon: the medoid record, `n_points` (distinct points),
#'   `n_used` (after subsampling and culling).
build_climate_medoids <- function(gbif_path, workers = 10, max_n = 20000,
                                  cull_percentage = 0.01, seed = 123) {
  points <- arrow::open_dataset(gbif_path) |>
    dplyr::filter(!is.na(taxon_name), !is.na(annual_precipitation),
                  !is.na(average_temperature)) |>
    dplyr::select(taxon_name, decimal_longitude, decimal_latitude,
                  annual_precipitation, average_temperature,
                  log_annual_precipitation_centered, average_temperature_centered) |>
    dplyr::collect() |>
    dplyr::distinct()

  one_taxon <- function(d) {
    n_points <- nrow(d)
    if (n_points > max_n) {
      set.seed(seed)
      d <- dplyr::slice_sample(d, n = max_n)
    }
    med <- locate_medoid_2d(d, cull_percentage = cull_percentage)
    dplyr::mutate(med$medoid, n_points = n_points, n_used = nrow(med$data_taxon))
  }

  by_taxon <- split(points, points$taxon_name)
  rm(points)
  # blocks of similar total work (n^2, capped at max_n), largest first
  work <- pmin(vapply(by_taxon, nrow, integer(1)), max_n)^2
  block <- integer(length(by_taxon))
  load <- numeric(workers * 4)
  for (i in order(work, decreasing = TRUE)) {
    b <- which.min(load); block[i] <- b; load[b] <- load[b] + work[i]
  }
  blocks <- split(by_taxon, block)

  old_plan <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  out <- furrr::future_map(blocks, function(taxa) {
    dplyr::bind_rows(lapply(taxa, one_taxon))
  }, .options = furrr::furrr_options(
    seed = NULL, packages = "dplyr",
    globals = list(one_taxon = one_taxon, locate_medoid_2d = locate_medoid_2d,
                   max_n = max_n, cull_percentage = cull_percentage, seed = seed)))
  dplyr::bind_rows(out) |> dplyr::arrange(taxon_name)
}
