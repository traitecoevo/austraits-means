# Round-2 simple species means: a fallback estimate for every species x trait,
# independent of the climate-centre model. Covers species/traits that round 1
# can't handle — too few records to fit, or no medoid (dropped from the model).
#
# These run on the *uncentred, un-medoid-filtered* record pool
# (load_height_records() + load_austraits_records()), so coverage is maximal.
# Traits are log-normal, so the central estimate is the geometric mean.

#' Per-species, per-trait simple summary (geometric mean + spread).
#'
#' @param records record-level data with `taxon_name`, `trait_name`, `value`,
#'   `log_value` (e.g. `bind_rows(load_height_records(), load_austraits_records())`).
#' @return tibble, one row per species x trait.
estimate_simple_means <- function(records) {
  records |>
    dplyr::filter(!is.na(log_value), is.finite(log_value)) |>
    dplyr::group_by(taxon_name, trait_name) |>
    dplyr::summarise(
      .groups = "drop",
      dplyr::across(dplyr::any_of(c("family", "genus", "unit")), dplyr::first),
      n = dplyr::n(),
      geom_mean = exp(mean(log_value)),
      sd_log = stats::sd(log_value),
      value_min = min(value),
      value_max = max(value)
    )
}

#' Read the round-1 climate-centre exports back into one long table.
#'
#' @param traits trait names to look for.
#' @param dir directory holding `<trait>_climate_response.csv.gz`.
load_centre_estimates <- function(traits, dir = "outputs/traits") {
  purrr::map_dfr(traits, function(tr) {
    f <- file.path(dir, paste0(clean_species_filename(tr), "_climate_response.csv.gz"))
    if (!file.exists(f)) return(tibble::tibble())
    readr::read_csv(f, show_col_types = FALSE) |>
      dplyr::transmute(
        taxon_name, trait_name = tr,
        centre_estimate = estimate,
        centre_conf.low = conf.low, centre_conf.high = conf.high,
        n_centre = n
      )
  })
}

#' Integrated best-estimate product: climate-centre where available, else simple.
#'
#' @param simple_means output of `estimate_simple_means()` (all species x trait).
#' @param centre_estimates output of `load_centre_estimates()`.
#' @return one row per species x trait with `estimate`, `method`, and both inputs.
combine_estimates <- function(simple_means, centre_estimates) {
  simple_means |>
    dplyr::left_join(centre_estimates, by = c("taxon_name", "trait_name")) |>
    dplyr::mutate(
      method = dplyr::if_else(!is.na(centre_estimate), "climate_centre", "simple_mean"),
      estimate = dplyr::coalesce(centre_estimate, geom_mean)
    ) |>
    dplyr::relocate(estimate, method, .after = trait_name)
}
