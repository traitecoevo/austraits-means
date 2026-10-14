# Modelled community mean surface (layer variant `model_surface`).
#
# The population-level prediction of each trait's hierarchical model at a
# cell's climate: the expected log trait of a species whose climate centre is
# the cell's climate, averaging over species (species effects and
# within-species terms set to zero):
#
#   log m(cell) = b0 + b_P * log P(cell) + b_T * T(cell)
#
# One surface per trait, for all species together, because the per-trait
# models are not split into woody and non-woody. Under the model the
# community variance is the between-species intercept variance, the same in
# every cell, so it is recorded rather than mapped. Cells whose climate lies
# outside the range of the species centres the model was fitted on are
# flagged as extrapolation.

#' Mean climate of each grid cell (WorldClim, averaged onto the grid).
#'
#' @param rasters climate rasters (`load_climate_rasters()`).
#' @return tibble(cell_x, cell_y, longitude, latitude, annual_precipitation,
#'   log_annual_precipitation, average_temperature), for cells touching
#'   Australian land.
cell_climate <- function(rasters, settings = community_layer_settings()) {
  e <- settings$extent
  template <- terra::rast(xmin = e[["xmin"]], xmax = e[["xmax"]],
                          ymin = e[["ymin"]], ymax = e[["ymax"]],
                          resolution = settings$resolution, crs = settings$crs)
  clim <- terra::resample(rasters, template, method = "average")
  names(clim) <- c("average_temperature", "annual_precipitation")
  # The rasters are cropped to Australia's bounding box, not its coastline:
  # keep cells touching Australian land only.
  australia <- terra::vect(ozmaps::ozmap_country) |> terra::project(settings$crs)
  clim <- terra::mask(clim, australia, touches = TRUE)
  terra::as.data.frame(clim, xy = TRUE, na.rm = TRUE) |>
    tibble::as_tibble() |>
    dplyr::filter(annual_precipitation > 0) |>
    dplyr::transmute(
      cell_x = floor(x / settings$resolution),
      cell_y = floor(y / settings$resolution),
      longitude = x, latitude = y,
      annual_precipitation,
      log_annual_precipitation = log(annual_precipitation),
      average_temperature
    )
}

#' Population-level prediction of each trait model at every cell's climate.
#'
#' @param fits list of TMB `trait_climate_fit` (the `tmb_fits` target).
#' @param climate output of `cell_climate()`.
#' @param level interval mass (delta method on the log scale).
#' @return one row per trait × cell: `community_mean` with `lower`/`upper`,
#'   `log_mean`, `extrapolated`, and the model's between-species SD
#'   (`sd_species`, log scale).
predict_model_surface <- function(fits, climate, level = 0.94) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  purrr::map_dfr(purrr::compact(fits), function(f) {
    stopifnot(f$backend == "tmb")
    b <- glmmTMB::fixef(f$fit)$cond
    V <- stats::vcov(f$fit)$cond[names(b), names(b)]
    X <- matrix(0, nrow(climate), length(b), dimnames = list(NULL, names(b)))
    X[, "(Intercept)"] <- 1
    X[, "log_annual_precipitation_mn"] <- climate$log_annual_precipitation
    X[, "average_temperature_mn"] <- climate$average_temperature
    eta <- as.numeric(X %*% b)
    se <- sqrt(rowSums((X %*% V) * X))

    p_range <- range(f$centres$log_annual_precipitation_mn)
    t_range <- range(f$centres$average_temperature_mn)
    sd_species <- attr(glmmTMB::VarCorr(f$fit)$cond[[f$group]], "stddev")[[1]]

    climate |>
      dplyr::mutate(
        trait_name = f$trait,
        log_mean = eta,
        community_mean = exp(eta),
        lower = exp(eta - z * se),
        upper = exp(eta + z * se),
        extrapolated = log_annual_precipitation < p_range[1] |
          log_annual_precipitation > p_range[2] |
          average_temperature < t_range[1] | average_temperature > t_range[2],
        sd_species = sd_species,
        .before = 1
      )
  })
}

#' Metadata written with the model surface.
model_surface_metadata <- function(data_basis, settings = community_layer_settings()) {
  c(list(
    product = "AusTraits community layers, modelled surface (first draft)",
    status = "provisional; one of several layer variants built for review (#20)",
    statistic = "community_mean: exp(b0 + b_P log P + b_T T) from each trait's hierarchical (TMB) model at the cell's mean climate, i.e. the expected trait of a species centred on that climate, averaging over species; lower/upper: 94% interval (delta method, fixed-effect uncertainty only)",
    community_variance = "under the model, the between-species variance sd_species^2 (log scale), the same in every cell; not mapped",
    groups = "all species together; the trait models are not split into woody and non-woody",
    climate = "WorldClim 2.1 BIO1 and BIO12 at 2.5 arc-minutes, averaged onto the grid",
    spatial_unit = paste0(settings$resolution, " degree grid, ", settings$crs,
                          "; longitude and latitude give cell centres"),
    extrapolated = "TRUE where the cell's climate lies outside the range of species climate centres in the fit"
  ), data_basis, repo_provenance())
}

#' Write the surface: parquet, one GeoTIFF per trait (community_mean), and
#' metadata. Returns the file paths.
write_model_surface <- function(surface, dir, metadata,
                                settings = community_layer_settings()) {
  dir.create(file.path(dir, "geotiff"), recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, c("model_surface.parquet", "model_surface_metadata.json"))
  arrow::write_parquet(surface, paths[1])
  sds <- dplyr::distinct(surface, trait_name, sd_species)
  jsonlite::write_json(c(metadata, list(sd_species = stats::setNames(
    as.list(sds$sd_species), sds$trait_name))), paths[2],
    auto_unbox = TRUE, pretty = TRUE)

  e <- settings$extent
  template <- terra::rast(xmin = e[["xmin"]], xmax = e[["xmax"]],
                          ymin = e[["ymin"]], ymax = e[["ymax"]],
                          resolution = settings$resolution, crs = settings$crs)
  tifs <- vapply(sort(unique(surface$trait_name)), function(tr) {
    d <- surface[surface$trait_name == tr, ]
    v <- rep(NA_real_, terra::ncell(template))
    v[terra::cellFromXY(template, cbind(d$longitude, d$latitude))] <- d$community_mean
    r <- template
    terra::values(r) <- v
    names(r) <- "community_mean"
    f <- file.path(dir, "geotiff", paste0(tr, "_all_community_mean_modelled.tif"))
    terra::writeRaster(r, f, overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
    f
  }, character(1))
  c(paths, unname(tifs))
}

#' The model surface beside the observed community mean of all species
#' (combined species means), on one colour scale. Extrapolated cells are grey.
plot_model_surface <- function(surface, observed_layers, trait, unit = NULL) {
  s <- dplyr::filter(surface, trait_name == trait)
  o <- dplyr::filter(observed_layers, trait_name == trait, woodiness_group == "all",
                     !is.na(community_mean))
  limits <- range(c(s$community_mean[!s$extrapolated], o$community_mean))
  oz <- ozmaps::ozmap_states
  blue <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
  panel <- function(d, title, grey = NULL) {
    p <- ggplot2::ggplot(d, ggplot2::aes(longitude, latitude))
    if (!is.null(grey)) {
      p <- p + ggplot2::geom_tile(data = grey, fill = "#e3e2de", width = 0.25, height = 0.25)
    }
    p +
      ggplot2::geom_tile(ggplot2::aes(fill = community_mean), width = 0.25, height = 0.25) +
      ggplot2::geom_sf(data = oz, inherit.aes = FALSE, fill = NA,
                       colour = "#6b6a66", linewidth = 0.2) +
      ggplot2::scale_fill_gradientn(colours = blue, transform = "log10",
                                    limits = limits, oob = scales::squish, name = NULL) +
      ggplot2::coord_sf(xlim = c(112, 154), ylim = c(-44, -9.5), expand = FALSE) +
      ggplot2::labs(title = title, x = NULL, y = NULL) +
      ggplot2::theme_void(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                     legend.key.height = grid::unit(1.2, "cm"))
  }
  unit_label <- if (is.null(unit) || is.na(unit)) "" else paste0(" (", unit, ")")
  patchwork::wrap_plots(
    panel(dplyr::filter(s, !extrapolated), "Modelled: model surface at cell climate",
          grey = dplyr::filter(s, extrapolated)),
    panel(o, "Observed: community mean of species present (combined means)"),
    ncol = 2, guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = paste0(trait, ", all species, community mean", unit_label),
      subtitle = "First draft, provisional choices; grey: climate outside the range of species centres in the model",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                             plot.subtitle = ggplot2::element_text(size = 9, colour = "#6b6a66")))
}

#' Save one model-surface map per trait; returns the file paths.
write_model_surface_maps <- function(surface, observed_layers, means, dir) {
  units <- dplyr::distinct(means, trait_name, unit)
  unlist(lapply(sort(unique(surface$trait_name)), function(tr) {
    p <- plot_model_surface(surface, observed_layers, tr,
                            units$unit[units$trait_name == tr])
    save_figure(p, file.path(dir, "figures", paste0(tr, "_model_surface")),
                width = 10, height = 5)
  }), use.names = FALSE)
}
