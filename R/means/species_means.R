# Numeric traits


#' Estimate Species Trait Means
#'
#' Calculates mean trait values for each species in the provided Austraits dataset.
#'
#' @param austraits A data frame or list containing trait data, typically in the Austraits format.
#'
#' @return A data frame with estimated mean trait values for each species.
#'
#' @details
#' This function processes the input Austraits dataset to compute mean values of traits for each species.
#'
#' @examples
#' \dontrun{
#' means <- estimate_species_trait_means(austraits)
#' }
#'
#' @export
estimate_species_trait_means <- function(austraits) {

  traits <-
    austraits |>
    select(trait_name, unit) |>
    dplyr::distinct() |>
    dplyr::mutate(
      type = ifelse(is.na(unit), "categorical", "numerical")
    ) |>
    split(~type)
  
  data_means_numerical <- 
    austraits |>
    estimate_species_trait_means_numerical(
      traits = traits$numerical$trait_name) # |>
      # dplyr::mutate(
      #   dplyr::across(tidyselect::where(is.numeric), as.character)
      # )
  
  data_means_categorical <-
    austraits |> estimate_species_trait_summary_categorical(
      traits = traits$categorical$trait_name)

  
  # Build species means for database
  data_means_numerical |>
    # combine the two datasets
    dplyr::bind_rows(data_means_categorical) |>
    # add in the trait type
    dplyr::left_join(by = c("trait_name"),
      traits |> dplyr::bind_rows() |> dplyr::select(trait_name, type),
    ) |>
    # add in taxon info
    dplyr::left_join(by = "taxon_name",
      austraits |>
        dplyr::select(taxon_name, taxon_rank:scientific_name_id, -taxon_name_alternatives) |>
        dplyr::distinct()
    )
}

#' @keywords internal
#' @noRd
estimate_species_trait_means_numerical <- function(austraits, traits) {

  # any data that is a mean, median or raw, create a site mean
  location_means <- austraits |> estimate_species_trait_means_locations(traits)

  # any data that is a max or a min (range) and basically from a flora, create a mean value
  flora_means <- austraits |> estimate_species_trait_means_floras(traits)
  
  # combine the two, then take means across site and flora replicates
  means <- location_means |>
    dplyr::bind_rows(flora_means)  |>
    dplyr::group_by(taxon_name, trait_name, unit) |>
    dplyr::summarise(.groups = "drop",
      value_mean = mean(value_mean),
      value_min = min(value_min),
      value_max = max(value_max),
      value_median = median(value_median),
      value_geom_mean = 10^mean(suppressWarnings(log10(value_mean)), na.rm= TRUE),
      all_replicates = sum(all_replicates),
      location_replicates = sum(location_replicates),
      flora_replicates = sum(flora_replicates),
      # record sources
      dataset_id = paste(unique(dataset_id), collapse = "; "),
      observation_id = paste(unique(observation_id), collapse = "; "),
    ) |>
    dplyr::distinct()
  
  means
}

#' @keywords internal
#' @noRd
estimate_species_trait_means_locations <- function(austraits, traits) {
  # any data that is a mean, median or raw, create a location mean
  x <- austraits |>
  dplyr::filter(
    trait_name %in% traits,
    value_type %in% c("mean", "raw", "median")
  ) |>
  dplyr::mutate(
    value = as.numeric(value),
    replicates = 1,
    log10_value = suppressWarnings(log10(value))
  ) |>
  dplyr::group_by(taxon_name, trait_name, dataset_id, location_id, unit) |>
  dplyr::summarise(
    .groups = "drop",
    dplyr::across(value, list(mean = mean, min = min, max = max, median = median), na.rm = TRUE),
    dplyr::across(c("latitude (deg)", "longitude (deg)", "location_name"), dplyr::first),
    all_replicates = sum(replicates),
    value_geom_mean = 10^mean(log10_value, na.rm= TRUE),
    observation_id = paste(unique(observation_id), collapse = "; "),
  ) |>
  dplyr::mutate(
    value_type = "location_mean",
    location_replicates = 1,
    flora_replicates = 0
  )
}

#' @keywords internal
#' @noRd
estimate_species_trait_means_floras <- function(austraits, traits) {
  flora_data <- 
    austraits |>
    dplyr::filter(
      trait_name %in% traits,
      value_type %in% c("minimum", "maximum"), 
      basis_of_record %in% c("preserved_specimen", "literature")
    )
  
  # early exit if no data available
  if(nrow(flora_data) == 0) {return(flora_data)}
  
  flora_data |>
  dplyr::mutate(
    value = as.numeric(value),
  )|>
  dplyr::group_by(taxon_name, trait_name, unit, dataset_id, observation_id, original_name) |> 
  dplyr::summarise(
    dplyr::across(value, list(mean = mean, min = min, max = max), na.rm = TRUE), 
  ) |>
  dplyr::mutate(
    location_replicates = 0,
    flora_replicates = 1,
    all_replicates = 1,
    value_median = value_mean,
    value_type = "flora_mean"
  ) |> 
  dplyr::ungroup()

}

#' @keywords internal
#' @noRd
estimate_species_trait_summary_categorical <- function(austraits, 
traits) {

  austraits |>
  estimate_species_trait_value_summary_categorical(traits) |>
  dplyr::group_by(taxon_name, trait_name) |>
  dplyr:: mutate(
      tmp_summary = paste0(value, " (", all_replicates, ")"),
      value_count = paste0(tmp_summary, collapse = "; ")
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-dplyr::all_of(c("tmp_summary", "value"))) |>
  dplyr::distinct() #|>
  # dplyr::mutate(
  #   dplyr::across(tidyselect::where(is.numeric), as.character)
  # )
}

#' @keywords internal
#' @noRd
estimate_species_trait_value_summary_categorical <- function(austraits, traits) {

  austraits |>
  dplyr::filter(trait_name %in% traits) |>
  dplyr::select(dplyr::all_of(c("dataset_id", "taxon_name", "trait_name", "location_id", "observation_id", "value")))|>
  dplyr::mutate(value = stringr::str_split(value, " ")) |>
  tidyr::unnest_longer(value) |>
  dplyr::mutate(
    all_replicates = 1
  ) |>
  dplyr::group_by(taxon_name, trait_name, value) |>
  summarise(.groups = "drop",
    value = first(value),
    all_replicates = sum(all_replicates),
    dataset_id = paste(unique(dataset_id), collapse = "; "),
    observation_id = paste(unique(observation_id), collapse = "; ")
  ) |>
  mutate(
    value_type = "value_count"
  )
}
