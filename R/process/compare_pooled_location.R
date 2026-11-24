# The pooled models' one remaining role (roadmap.md, 25 November 2026): to
# show what the location model fixes. Three nested models, named by what they
# pool (terminology.md):
#
#   pooled            growth forms pooled, locations pooled  (`tmb_fit`)
#   locations pooled  growth forms apart,  locations pooled  (`woody_fit`)
#   location model    growth forms apart,  locations apart   (`location_fit`)
#
# Per trait and growth form: how far the species means move at each step,
# and the community SD among species within a cell that the two
# per-growth-form models imply against what is observed (the benchmark of
# docs/local-range-model.qmd). The pooled model has one species variance for
# both growth forms, so it has no per-growth-form community SD.

#' @param m a fit_location_target() result.
#' @param pooled the trait's pooled `trait_climate_fit` (with centres).
#' @param locations_pooled the trait's woodiness model without location
#'   effects (fit_woodiness_model(), the `woody_fit` target).
#' @return list(trait, shift, benchmark) or NULL. `shift`: per growth form,
#'   the species in all three models and, for each step, the SD and median
#'   absolute change in log mean: `forms` (pooled to locations pooled),
#'   `locations` (locations pooled to location model) and `total`, with the
#'   correlation of the pooled and location-model means. `benchmark`: per
#'   model and growth form, the ratio of modelled to observed community SD
#'   (summarise_community_benchmark()); with `sorting`, also the location
#'   model with its species spread scaled by the sorting factor.
compare_pooled_location <- function(m, pooled, locations_pooled, sorting = NULL,
                                    min_species = 5) {
  if (is.null(m) || is.null(pooled) || is.null(locations_pooled)) return(NULL)
  means <- dplyr::transmute(m$centres, taxon_name,
                            woodiness_group = as.character(woodiness_group),
                            log_location = log_estimate) |>
    dplyr::inner_join(dplyr::transmute(species_centre_values(locations_pooled), taxon_name,
                                       log_forms = log_centre),
                      by = "taxon_name") |>
    dplyr::inner_join(dplyr::transmute(pooled$centres, taxon_name, log_pooled = log_estimate),
                      by = "taxon_name")
  shift <- means |>
    dplyr::group_by(woodiness_group) |>
    dplyr::summarise(
      .groups = "drop", n_species = dplyr::n(),
      sd_change_forms = stats::sd(log_forms - log_pooled),
      median_abs_change_forms = stats::median(abs(log_forms - log_pooled)),
      sd_change_locations = stats::sd(log_location - log_forms),
      median_abs_change_locations = stats::median(abs(log_location - log_forms)),
      sd_change_total = stats::sd(log_location - log_pooled),
      median_abs_change_total = stats::median(abs(log_location - log_pooled)),
      correlation_total = stats::cor(log_location, log_pooled)) |>
    dplyr::mutate(trait = m$trait, .before = 1)

  # the community benchmark needs each record's cell
  locations_pooled$data <- add_location_key(locations_pooled$data, "cell",
                                            cell_km = location_cell_km())
  noise <- tidy_varcomp(m)
  benchmark <- dplyr::bind_rows(
    `location model` = summarise_community_benchmark(
      benchmark_community_variance(m, noise_vc = noise, min_species = min_species)),
    `location model, sorting` = if (!is.null(sorting)) summarise_community_benchmark(
      benchmark_community_variance(.with_sorting(m, sorting), noise_vc = noise,
                                   min_species = min_species)),
    `locations pooled` = summarise_community_benchmark(
      benchmark_community_variance(locations_pooled, noise_vc = noise,
                                   min_species = min_species)),
    .id = "model") |>
    dplyr::mutate(trait = m$trait, .before = 1)
  list(trait = m$trait, shift = shift, benchmark = benchmark)
}

# A location-model fit with each species' variance scaled by its growth
# form's sorting factor (community_sorting.R), for the benchmark.
.with_sorting <- function(m, sorting) {
  g <- as.character(m$centres$woodiness_group[match(m$species_variance$taxon_name,
                                                    m$centres$taxon_name)])
  k <- vapply(g, function(x) .sorting_k(sorting, m$trait, x), numeric(1))
  m$species_variance$species_variance <- m$species_variance$species_variance * k
  m
}

#' Write the comparison across traits to `dir`.
#' @return the file paths.
write_pooled_location_comparison <- function(cmp, dir = "outputs/comparisons") {
  cmp <- purrr::compact(cmp)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, c("pooled_vs_location_species_means.csv",
                            "pooled_vs_location_community_sd.csv"))
  readr::write_csv(purrr::map(cmp, "shift") |> purrr::list_rbind(), paths[1])
  readr::write_csv(purrr::map(cmp, "benchmark") |> purrr::list_rbind(), paths[2])
  paths
}
