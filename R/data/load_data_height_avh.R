assign_state_fast <- function(data) {
  data |>
    dplyr::mutate(
      vic_nsw_border_lat = -34 - 0.25 * (decimal_longitude - 141),
      state = dplyr::case_when(
        is.na(decimal_longitude) | is.na(decimal_latitude) ~ NA_character_,
        decimal_latitude < -39 ~ "TAS",
        dplyr::between(decimal_longitude, 148.75, 149.40) & dplyr::between(decimal_latitude, -35.95, -35.10) ~ "ACT",
        decimal_longitude < 129 ~ "WA",
        decimal_longitude < 138 & decimal_latitude > -26 ~ "NT",
        decimal_longitude < 141 ~ "SA",
        decimal_latitude > -29 ~ "QLD",
        decimal_latitude < vic_nsw_border_lat ~ "VIC",
        TRUE ~ "NSW"
      )
    ) |>
    dplyr::select(-vic_nsw_border_lat) |>
    dplyr::relocate(state, .after = decimal_latitude)
}

load_data <- function(avh_csv_path = NULL) {
  default_csv <- "data/AVH_height_records_20260311.csv"
  env_csv <- Sys.getenv("AVH_HEIGHT_CSV", unset = "")

  # Path priority: function argument > env var > project default
  avh_csv <- dplyr::coalesce(avh_csv_path, if (nzchar(env_csv)) env_csv else NA_character_, default_csv)

  if (!file.exists(avh_csv)) {
    stop(
      "AVH file not found: ", avh_csv, "\n",
      "Set `avh_csv_path=` in `load_data()` or env var `AVH_HEIGHT_CSV`."
    )
  }

  if (file.access(avh_csv, 4) != 0) {
    stop(
      "AVH file exists but is not readable: ", avh_csv, "\n",
      "On macOS this is often a Files/OneDrive permission issue. ",
      "Move/copy the file into the project data folder or grant terminal/R full disk access."
    )
  }

  data <-
    tryCatch(
      readr::read_csv(avh_csv, show_col_types = FALSE),
      error = function(e) {
        stop(
          "Failed to read AVH CSV: ", avh_csv, "\n",
          "Original error: ", conditionMessage(e), "\n",
          "If this points to OneDrive/CloudStorage, copy the CSV locally and pass `avh_csv_path=`."
        )
      }
    ) |>
    dplyr::mutate(trait_name = "height", log_value = log(height_final)) |>
    dplyr::select(family, genus, taxon_name,
      recordNumber, row_number, identifier, institutionCode, year,
      trait_name,
      value_type,
      value = height_final,
      log_value,
      decimal_longitude = decimalLongitude, decimal_latitude = decimalLatitude,
      original_name = species,
    ) |>
    dplyr::filter(dplyr::if_any(where(is.double), ~ !(is.na(.x))), value > 0)

  vars <- c("wc2.1_2.5m_bio_1.tif", "wc2.1_2.5m_bio_12.tif")

  australia <- terra::rast("data/australia.tif")

  stacked_rasters <-
    file.path("data/climate/", vars) |>
    terra::rast() |>
    terra::crop(australia)

  austraits_geo <-
    data |>
    dplyr::filter(
      !is.na(decimal_latitude), !is.na(decimal_longitude)
    )

  trait_data_vect <-
    austraits_geo |>
    dplyr::select("decimal_longitude", "decimal_latitude") |>
    terra::vect(
      geom = c("decimal_longitude", "decimal_latitude"),
      crs = terra::crs(stacked_rasters)
    )

  data_traits_climate <-
    terra::extract(x = stacked_rasters, y = trait_data_vect, method = "simple", bind = TRUE, cells = FALSE, ID = FALSE, dataframe = TRUE) |>
    tibble::as_tibble() |>
    dplyr::bind_cols(austraits_geo) |>
    dplyr::rename(
      annual_precipitation = wc2.1_2.5m_bio_12, average_temperature = wc2.1_2.5m_bio_1
    ) |>
    dplyr::select(decimal_longitude, decimal_latitude, annual_precipitation, average_temperature, dplyr::everything()) |>
    dplyr::mutate(
      log_annual_precipitation = log(annual_precipitation),
    ) |>
    dplyr::filter(!is.na(log_annual_precipitation) & !is.na(average_temperature))

  data_traits_climate <- assign_state_fast(data_traits_climate)

  climate_centres <- data_traits_climate |>
    dplyr::group_by(taxon_name) |>
    dplyr::summarise(
      .groups = "drop",
      dplyr::across(c(log_annual_precipitation, average_temperature, decimal_latitude, decimal_longitude), mean, .names = "{col}_mn")
    )

  data_traits_climate |>
    dplyr::left_join(climate_centres, by = "taxon_name") |>
    dplyr::filter(!is.na(log_annual_precipitation_mn)) |>
    dplyr::group_by(taxon_name) |>
    dplyr::mutate(
      log_annual_precipitation_c = log_annual_precipitation - log_annual_precipitation_mn,
      average_temperature_c = average_temperature - average_temperature_mn
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(taxon_name, trait_name) |>
    dplyr::mutate(
      log_value_c = log_value - mean(log_value)
    ) |>
    dplyr::ungroup()
}
