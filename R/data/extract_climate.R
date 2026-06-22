# Shared WorldClim climate extraction.
#
# Two bioclimatic layers are used throughout the project:
#   BIO_1  -> average_temperature   (Annual Mean Temperature, degC)
#   BIO_12 -> annual_precipitation  (Annual Precipitation, mm/yr)
# Both at 2.5 arc-minute resolution, cropped to Australia.

#' Load the stacked WorldClim rasters, cropped to Australia.
#'
#' @return a terra SpatRaster with layers wc2.1_2.5m_bio_1 and wc2.1_2.5m_bio_12.
load_climate_rasters <- function(
    vars = c("wc2.1_2.5m_bio_1.tif", "wc2.1_2.5m_bio_12.tif"),
    climate_dir = "data/climate",
    australia_tif = "data/australia.tif") {
  australia <- terra::rast(australia_tif)
  file.path(climate_dir, vars) |>
    terra::rast() |>
    terra::crop(australia)
}

#' Extract average_temperature + annual_precipitation at point coordinates.
#'
#' Returns the input `data` with `annual_precipitation`, `average_temperature`
#' and `log_annual_precipitation` columns prepended. Rows with missing climate
#' (e.g. offshore points) are kept; filter downstream as needed.
#'
#' @param data data frame with longitude/latitude columns.
#' @param lon,lat column names holding longitude/latitude.
#' @param rasters a SpatRaster from `load_climate_rasters()` (loaded by default).
extract_climate <- function(data,
    lon = "decimal_longitude", lat = "decimal_latitude",
    rasters = load_climate_rasters()) {

  pts <- data |>
    dplyr::select(dplyr::all_of(c(lon, lat))) |>
    terra::vect(geom = c(lon, lat), crs = terra::crs(rasters))

  climate <-
    terra::extract(x = rasters, y = pts, method = "simple",
                   ID = FALSE, xy = FALSE, dataframe = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(
      annual_precipitation = wc2.1_2.5m_bio_12,
      average_temperature = wc2.1_2.5m_bio_1
    ) |>
    dplyr::mutate(log_annual_precipitation = log(annual_precipitation))

  dplyr::bind_cols(climate, data)
}
