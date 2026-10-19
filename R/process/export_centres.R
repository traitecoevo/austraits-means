# Export the modelled means in the product schema of terminology.md: per
# trait as each trait is fitted (`<trait>_climate_response.csv.gz` +
# `_meta.csv` in outputs/traits/), and as one product table across traits
# (`modelled_means_table()`, written to outputs/species_means/<version>/).
# The species centres inside a fit keep their working names (`estimate`,
# `log_estimate`, ...); only the exports use the product names.

#' Species centres of one fit in the product schema.
#'
#' @param centres output of `species_centres()`.
#' @param trait trait name.
#' @param unit the trait's unit.
#' @param backend "tmb" or "brms" (becomes `method`).
#' @param data the fit's records, for the reliability indicators (#16).
centres_to_modelled_means <- function(centres, trait, unit, backend, data) {
  centres |>
    dplyr::transmute(
      taxon_name, trait_name = trait, unit = unit,
      estimate_type = "modelled_mean", method = backend,
      value = estimate, value_lower = conf.low, value_upper = conf.high,
      log_value = log_estimate, log_value_lower = log_conf.low,
      log_value_upper = log_conf.high,
      n_records = n,
      annual_precipitation_centre, average_temperature_centre,
      species_intercept = g_0,
      species_slope_precipitation = g_log_annual_precipitation_c,
      species_slope_temperature = g_average_temperature_c,
      family, genus
    ) |>
    dplyr::left_join(modelled_reliability(data), by = "taxon_name") |>
    dplyr::relocate(n_datasets, n_locations, n_flora, n_field, sd_log,
                    .after = n_records) |>
    dplyr::mutate(reliability = reliability_class(n_datasets, n_records),
                  .after = sd_log)
}

#' Modelled means of every fitted trait, in the product schema.
#'
#' @param fits list of `trait_climate_fit` (e.g. the `tmb_fits` target).
modelled_means_table <- function(fits) {
  purrr::map_dfr(purrr::compact(fits), \(f) centres_to_modelled_means(
    f$centres, f$trait, unit = f$data$unit[1], backend = f$backend,
    data = f$data))
}

#' Column dictionary for the modelled means.
modelled_means_dictionary <- function() {
  tibble::tribble(
    ~column, ~definition,
    "taxon_name", "Species name as aligned in the AusTraits build, with subspecies, varieties and forms rolled up to their species.",
    "trait_name", "AusTraits trait name (APD).",
    "unit", "Unit of value.",
    "estimate_type", "Always modelled_mean: the model-estimated geometric mean at the centre of the species' climate distribution (terminology.md).",
    "method", "Fitting backend: tmb (glmmTMB, maximum likelihood) or brms (Stan, full posterior).",
    "value", "The modelled mean: exp of the model's prediction on the log scale for the species at its climate centre (medoid), where within-species climate deviations are zero.",
    "value_lower", "Lower bound of the 94% interval of value (TMB: exp(fit - z se), delta method).",
    "value_upper", "Upper bound of the 94% interval of value.",
    "log_value", "Natural log of value: the prediction on the model scale.",
    "log_value_lower", "Lower bound of the 94% interval on the log scale.",
    "log_value_upper", "Upper bound of the 94% interval on the log scale.",
    "annual_precipitation_centre", "Mean annual precipitation (mm, WorldClim BIO12) at the species' climate centre, the medoid of its GBIF occurrences.",
    "average_temperature_centre", "Mean annual temperature (degC, WorldClim BIO1) at the species' climate centre.",
    "species_intercept", "The species' random intercept (log scale): its departure from the across-species prediction at its climate centre.",
    "species_slope_precipitation", "The species' random slope for within-species log precipitation (departure from the mean within-species slope).",
    "species_slope_temperature", "The species' random slope for within-species temperature (departure from the mean within-species slope).",
    "family", "Family.",
    "genus", "Genus."
  ) |>
    dplyr::bind_rows(
      reliability_dictionary() |>
        dplyr::mutate(definition = dplyr::if_else(
          column == "n_records",
          "Number of field records of the species used in the fit (species need more than 5).",
          definition))
    )
}

#' Combined species means: for each modelled trait, the modelled mean where
#' there is one, otherwise the geometric observed mean. `estimate_type` and
#' `method` say which each row is.
combine_species_means <- function(fits, observed_means) {
  modelled <- modelled_means_table(fits) |>
    dplyr::select(taxon_name, trait_name, unit, estimate_type, method, value,
                  value_lower, value_upper, log_value, n_records)
  observed <- observed_means |>
    dplyr::filter(type == "numerical", trait_name %in% unique(modelled$trait_name)) |>
    dplyr::transmute(taxon_name, trait_name, unit, estimate_type = "observed_mean",
                     method = "geometric", value, log_value, n_records)

  dplyr::bind_rows(
    modelled,
    dplyr::anti_join(observed, modelled, by = c("taxon_name", "trait_name"))
  ) |>
    dplyr::arrange(trait_name, taxon_name)
}

#' Column dictionary for the combined species means.
species_means_combined_dictionary <- function() {
  dplyr::bind_rows(
    modelled_means_dictionary() |>
      dplyr::filter(column %in% c("taxon_name", "trait_name", "unit", "value_lower",
                                  "value_upper", "log_value")),
    tibble::tribble(
      ~column, ~definition,
      "estimate_type", "modelled_mean where the species has one, otherwise observed_mean (terminology.md).",
      "method", "tmb for modelled means; geometric for observed means.",
      "value", "The species mean: the modelled mean where there is one, otherwise the geometric observed mean. Both are geometric means.",
      "n_records", "Records behind the value: field records in the model (modelled means) or all records used (observed means)."
    )
  )
}

#' Provenance written with the modelled means.
modelled_means_metadata <- function(data_basis) {
  c(list(
    product = "AusTraits species means",
    version = data_basis$austraits_version,
    estimate_type = "modelled_mean",
    method = "tmb",
    model = paste(deparse(build_model_formula()), collapse = " "),
    centring = "climate predictors centred on each species' medoid: x_mn = medoid, x_c = record - medoid",
    species_included = "species with more than 5 field records and a climate medoid",
    interval = "94%, delta method on the log scale",
    documentation = "terminology.md; analysis.qmd"
  ), data_basis, repo_provenance())
}

#' Provenance written with the combined species means.
species_means_combined_metadata <- function(data_basis) {
  c(list(
    product = "AusTraits species means",
    version = data_basis$austraits_version,
    content = "one species mean per species and trait: the modelled mean (TMB) where there is one, otherwise the geometric observed mean",
    documentation = "terminology.md; docs/observed-means-geometric-method.md"
  ), data_basis, repo_provenance())
}

#' Write one trait's modelled means (`<trait>_climate_response.csv.gz`) and
#' its column dictionary (`_meta.csv`) to `out_dir`.
#'
#' @return (invisibly) the data CSV path.
export_species_centres <- function(centres, trait, out_dir = "outputs/traits",
                                   unit = NA_character_, backend = "tmb", data) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  table <- centres_to_modelled_means(centres, trait, unit, backend, data)
  stem <- file.path(out_dir, paste0(clean_species_filename(trait), "_climate_response"))
  readr::write_csv(table, paste0(stem, ".csv.gz"))
  readr::write_csv(
    modelled_means_dictionary() |> dplyr::filter(column %in% names(table)),
    paste0(stem, "_meta.csv"))
  invisible(paste0(stem, ".csv.gz"))
}
