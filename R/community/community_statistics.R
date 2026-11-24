# Community statistics (roadmap.md, 25 November 2026): for a location and
# its climate, the woody and non-woody community mean and community range
# (terminology.md), from the location model (R/model/spread_rtmb.R). They
# replace the three 10 km community layers (observed local, observed
# presence, modelled) and their comparisons, dropped on 25 November 2026.
#
# A location is a cell of the location model's grid (location_cell_km(), 50
# km). The species of a cell are those with GBIF records there, each placed
# in its growth form (#19) and at its climate centre (the medoid of its
# occurrences); they need no trait data. For species s of growth form g with
# centre m_s, at a cell of climate x, the location model gives the expected
# species mean at an average location
#   mu_s = b0_g + b_between,g . m_s + b_within,g . (x - m_s)
# and the spread of the species means there (local_range.R, header) is
#   var_s(mu_s) + mean_s(sd_u(s)^2) + sd_w^2 + mean_s(sum_k sd_k^2 (x_k - m_sk)^2)
# with sd_u(s) the species SD at its centre (shrunk-linear in climate, kept
# within the range of SDs the model gives its fitted species of that growth
# form, so the slopes in log SD are not extrapolated), sd_w the species x location SD and sd_k the
# SDs of the species' climate slopes. The community mean is exp(mean_s mu_s),
# the community range the central 90% interval of a normal on the log scale
# with that variance.
#
# community_statistics() takes the species present from each cell's own
# GBIF records; it is the check on the product, community_statistics_climate()
# (community_climate.R), which predicts from climate alone.

#' Grid settings: the location model's cells, in an equal-area projection.
community_layer_settings <- function() {
  list(
    resolution = location_cell_km() * 1000, min_species = 5,
    crs = "+proj=laea +lat_0=-25.2744 +lon_0=133.7751 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
    # Australia's occurrences in this projection, a multiple of 10 and 50 km
    extent = c(xmin = -2300000, xmax = 2200000, ymin = -2200000, ymax = 1800000)
  )
}

#' Empty raster on the grid.
community_grid <- function(settings = community_layer_settings()) {
  e <- settings$extent
  terra::rast(xmin = e[["xmin"]], xmax = e[["xmax"]], ymin = e[["ymin"]],
              ymax = e[["ymax"]], resolution = settings$resolution, crs = settings$crs)
}

#' Grid cell of each longitude / latitude point.
#' @return tibble(cell_x, cell_y): cell (i, j) spans [i, i + 1) x resolution
#'   metres east and likewise north in the grid projection.
lonlat_to_cell <- function(lon, lat, settings = community_layer_settings()) {
  xy <- sf::sf_project("EPSG:4326", settings$crs, cbind(lon, lat))
  tibble::tibble(cell_x = as.integer(floor(xy[, 1] / settings$resolution)),
                 cell_y = as.integer(floor(xy[, 2] / settings$resolution)))
}

#' Add cell-centre coordinates: `x`, `y` (metres, grid projection) and
#' `longitude`, `latitude`.
add_cell_centres <- function(d, settings = community_layer_settings()) {
  x <- (d$cell_x + 0.5) * settings$resolution
  y <- (d$cell_y + 0.5) * settings$resolution
  ll <- sf::sf_project(settings$crs, "EPSG:4326", cbind(x, y))
  dplyr::mutate(d, x = x, y = y, longitude = ll[, 1], latitude = ll[, 2],
                .after = cell_y)
}

#' Species present per grid cell, from the GBIF occurrences.
#' @return tibble(gbif_taxon_name, cell_x, cell_y).
gbif_presence_cells <- function(gbif_path, settings = community_layer_settings()) {
  occ <- arrow::open_dataset(gbif_path) |>
    dplyr::select(taxon_name, decimal_longitude, decimal_latitude) |>
    dplyr::filter(!is.na(taxon_name)) |>
    dplyr::collect()
  dplyr::bind_cols(dplyr::select(occ, gbif_taxon_name = taxon_name),
                   lonlat_to_cell(occ$decimal_longitude, occ$decimal_latitude, settings)) |>
    dplyr::distinct()
}

#' One species name (as used by the species means) per GBIF name.
#'
#' A GBIF name that is itself a species-means name maps to itself; otherwise
#' the crosswalk, with its trait-side names rolled up to species, is used,
#' preferring the same name, then the accepted name, then the first
#' alphabetically.
#'
#' @param gbif_names GBIF occurrence names.
#' @param crosswalk output of `build_name_crosswalk()`.
#' @param taxa the build's taxa table (taxon_name, taxon_rank, binomial).
#' @param species_names the names to map to (those with a woody class).
gbif_to_species <- function(gbif_names, crosswalk, taxa, species_names) {
  via_crosswalk <- crosswalk |>
    roll_up_to_species(taxa) |>
    dplyr::filter(taxon_name %in% species_names) |>
    dplyr::mutate(.rank = dplyr::case_when(
      taxon_name_aligned == gbif_taxon_name ~ 1L,
      taxon_name == accepted_name ~ 2L,
      TRUE ~ 3L)) |>
    dplyr::select(gbif_taxon_name, taxon_name, .rank)
  direct <- tibble::tibble(gbif_taxon_name = intersect(gbif_names, species_names)) |>
    dplyr::mutate(taxon_name = gbif_taxon_name, .rank = 0L)

  dplyr::bind_rows(direct, via_crosswalk) |>
    dplyr::arrange(gbif_taxon_name, .rank, taxon_name) |>
    dplyr::distinct(gbif_taxon_name, .keep_all = TRUE) |>
    dplyr::select(gbif_taxon_name, taxon_name)
}

#' Mean climate of each grid cell (WorldClim averaged onto the grid, cells
#' touching Australian land only).
#'
#' @param rasters climate rasters (`load_climate_rasters()`).
#' @return tibble(cell_x, cell_y, annual_precipitation,
#'   log_annual_precipitation, average_temperature).
cell_climate <- function(rasters, settings = community_layer_settings()) {
  template <- community_grid(settings)
  clim <- terra::project(rasters, template, method = "average")
  names(clim) <- c("average_temperature", "annual_precipitation")
  australia <- terra::vect(sf::st_transform(ozmaps::ozmap_country, settings$crs))
  clim <- terra::mask(clim, australia, touches = TRUE)
  terra::as.data.frame(clim, xy = TRUE, na.rm = TRUE) |>
    tibble::as_tibble() |>
    dplyr::filter(annual_precipitation > 0) |>
    dplyr::transmute(
      cell_x = as.integer(floor(x / settings$resolution)),
      cell_y = as.integer(floor(y / settings$resolution)),
      annual_precipitation,
      log_annual_precipitation = log(annual_precipitation),
      average_temperature
    )
}

#' The species of each cell, with growth form and climate centre.
#'
#' @param presence output of `gbif_presence_cells()` (GBIF names x cells).
#' @param medoids_path medoid parquet, keyed by GBIF name.
#' @param name_map,groups to place each GBIF species in woody / non-woody;
#'   species without a class are left out.
#' @return tibble(cell_x, cell_y, gbif_taxon_name, woodiness_group,
#'   centre_log_annual_precipitation, centre_average_temperature).
cell_species <- function(presence, medoids_path, name_map, groups) {
  centres <- arrow::read_parquet(medoids_path) |>
    dplyr::filter(!is.na(annual_precipitation), annual_precipitation > 0) |>
    dplyr::transmute(gbif_taxon_name = taxon_name,
                     centre_log_annual_precipitation = log(annual_precipitation),
                     centre_average_temperature = average_temperature)
  grp <- name_map |>
    dplyr::inner_join(dplyr::filter(groups, woodiness_group %in% c("woody", "non_woody")) |>
                        dplyr::select(taxon_name, woodiness_group), by = "taxon_name") |>
    dplyr::distinct(gbif_taxon_name, woodiness_group)
  presence |>
    dplyr::inner_join(centres, by = "gbif_taxon_name") |>
    dplyr::inner_join(grp, by = "gbif_taxon_name") |>
    dplyr::relocate(woodiness_group, .after = gbif_taxon_name)
}

#' Community mean and range per trait x cell x growth form, from the
#' location model (see the file header).
#'
#' @param fits list of fit_location_target() results.
#' @param climate output of `cell_climate()`.
#' @param species output of `cell_species()`.
#' @param level mass of the community range.
#' @return long tibble: trait_name, cell, woodiness_group, n_species, cell
#'   climate, community_mean, community_sd (natural-log units), range_lower,
#'   range_upper (empty with fewer than `min_species` species), and
#'   `extrapolated` (the cell's climate outside the records' for that trait
#'   and growth form).
community_statistics <- function(fits, climate, species, level = 0.90,
                                 settings = community_layer_settings()) {
  cells <- dplyr::inner_join(species, climate, by = c("cell_x", "cell_y"))
  purrr::map(purrr::compact(fits), function(m) {
    vc <- m$vc
    sd_slope <- function(term) vc$sd[vc$grouping == "taxon_name" & vc$term == term]
    sP <- sd_slope("log_annual_precipitation_c")
    sT <- sd_slope("average_temperature_c")
    sw <- vc$sd[vc$grouping == "taxon_name:location"]
    s <- m$centre_scaling
    fitted_sp <- dplyr::distinct(m$data, taxon_name, woodiness_group)
    purrr::map(m$groups, function(g) {
      b <- m$fixed[.group_coef_names(m, g)]
      names(b) <- c("b0", "bP", "bT", "wP", "wT")
      d <- dplyr::filter(cells, woodiness_group == g)
      if (nrow(d) == 0) return(NULL)
      mP <- d$centre_log_annual_precipitation
      mT <- d$centre_average_temperature
      xP <- d$log_annual_precipitation
      xT <- d$average_temperature
      mu <- b[["b0"]] + b[["bP"]] * mP + b[["bT"]] * mT + b[["wP"]] * (xP - mP) + b[["wT"]] * (xT - mT)
      # species SD at the centre, kept within the range of the fitted
      # species' SDs: clamping each climate axis alone still lets the two
      # slopes compound at a corner no fitted species occupies (woody seed
      # mass reached SD 5.2 at the warm, wet corner against 3.5 at most)
      p <- (mP - s$p_mean) / s$p_sd
      t <- (mT - s$t_mean) / s$t_sd
      th <- m$spread$theta
      log_sd <- th[[paste0(g, ":(Intercept)")]] +
        Reduce(`+`, lapply(m$spread$terms, function(k)
          th[[paste0(g, ":", k)]] * .spread_term_cols(p, t, k)[, 1]), 0)
      fitted_log_sd <- 0.5 * log(m$species_variance$species_variance[
        m$species_variance$taxon_name %in% fitted_sp$taxon_name[fitted_sp$woodiness_group == g]])
      log_sd <- pmin(pmax(log_sd, min(fitted_log_sd)), max(fitted_log_sd))
      rec <- dplyr::filter(m$data, woodiness_group == g)
      outside <- function(v, r) v < min(r) | v > max(r)
      sm <- d |>
        dplyr::mutate(mu = mu, var_u = exp(2 * log_sd),
                      var_slopes = sP^2 * (xP - mP)^2 + sT^2 * (xT - mT)^2) |>
        dplyr::group_by(cell_x, cell_y, annual_precipitation, log_annual_precipitation,
                        average_temperature) |>
        dplyr::summarise(.groups = "drop", n_species = dplyr::n(),
                         mean_log = mean(mu),
                         variance_log = mean((mu - mean(mu))^2) + mean(var_u) + sw^2 +
                           mean(var_slopes)) |>
        dplyr::mutate(
          trait_name = m$trait, woodiness_group = g, .before = 1,
          extrapolated = outside(log_annual_precipitation, rec$log_annual_precipitation) |
            outside(average_temperature, rec$average_temperature))
      sm |>
        dplyr::bind_cols(community_range(sm$mean_log, sm$variance_log, level)) |>
        dplyr::mutate(dplyr::across(c(community_mean, community_sd, range_lower, range_upper),
                                    \(v) dplyr::if_else(n_species >= settings$min_species,
                                                        v, NA_real_))) |>
        dplyr::select(-mean_log, -variance_log)
    }) |> purrr::list_rbind()
  }) |>
    purrr::list_rbind() |>
    add_cell_centres(settings) |>
    dplyr::arrange(trait_name, woodiness_group, cell_y, cell_x)
}

#' Write the community statistics: a long parquet, one GeoTIFF per trait x
#' growth form x statistic, and a metadata file.
#' @param statistics columns written as GeoTIFFs.
#' @return the file paths.
write_community_layers <- function(layers, dir, metadata,
                                   statistics = c("community_mean", "community_sd",
                                                  "range_lower", "range_upper", "n_species"),
                                   settings = community_layer_settings()) {
  dir.create(file.path(dir, "geotiff"), recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, c("community_statistics.parquet", "community_statistics_metadata.json"))
  arrow::write_parquet(layers, paths[1])
  jsonlite::write_json(c(metadata, list(n_rows = nrow(layers))), paths[2],
                       auto_unbox = TRUE, pretty = TRUE)

  template <- community_grid(settings)
  combos <- dplyr::distinct(layers, trait_name, woodiness_group)
  tifs <- unlist(purrr::pmap(combos, function(trait_name, woodiness_group) {
    d <- layers[layers$trait_name == trait_name & layers$woodiness_group == woodiness_group, ]
    cells <- terra::cellFromXY(template, cbind(d$x, d$y))
    vapply(statistics, function(st) {
      r <- template
      v <- rep(NA_real_, terra::ncell(r))
      v[cells] <- as.numeric(d[[st]])
      terra::values(r) <- v
      names(r) <- st
      f <- file.path(dir, "geotiff", paste0(trait_name, "_", woodiness_group, "_", st, ".tif"))
      terra::writeRaster(r, f, overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
      f
    }, character(1))
  }))
  c(paths, unname(tifs))
}

#' Metadata of the community statistics.
community_statistics_metadata <- function(data_basis, settings = community_layer_settings()) {
  c(list(
    product = "AusTraits community statistics (unpublished draft)",
    description = paste(
      "For each cell and growth form, from the cell's climate alone, the community mean and",
      "community range from the location model (R/model/spread_rtmb.R,",
      "R/community/community_climate.R): the species' expected means at the cell's climate",
      "and an average location, averaged over the climate centres of the species that",
      "co-occur at that climate (a smooth surface learned from GBIF records); the range from",
      "the model's variance among species within a location."),
    weighting = "presence only: every species counts equally (#18)",
    groups = "woody / non-woody by the rule of #19; species without a class are left out",
    spatial_unit = paste0(location_cell_km(), " x ", location_cell_km(),
                          " km cells of an equal-area grid (", settings$crs,
                          "), the location model's locations; x and y give cell centres in",
                          " metres, longitude and latitude in degrees"),
    extrapolation = paste("`extrapolated`: the cell's climate is outside that of the trait's records",
                          "for the growth form; `centres_extrapolated`: outside the cells the centre",
                          "surfaces were fitted to"),
    definitions = list(
      community_mean = "exp(mean over the species of their expected log mean); a geometric mean, in trait units",
      community_sd = "SD among the species' means within the cell, natural-log units",
      range_lower = "lower end of the central 90% of the species' means, trait units",
      range_upper = "upper end of the central 90% of the species' means, trait units",
      extrapolated = "the cell's climate lies outside that of the trait's records for the growth form")
  ), data_basis, repo_provenance())
}

#' A map panel of one statistic on the grid, with state outlines. Cells with
#' the statistic empty (species present, too few) are drawn grey.
map_cells <- function(d, stat, colours, title, transform = "identity",
                      settings = community_layer_settings(), facet = NULL,
                      limits = NULL) {
  oz <- sf::st_transform(ozmaps::ozmap_states, settings$crs)
  r <- settings$resolution
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_tile(data = dplyr::filter(d, is.na(.data[[stat]])),
                       fill = "#e3e2de", width = r, height = r) +
    ggplot2::geom_tile(data = dplyr::filter(d, !is.na(.data[[stat]])),
                       ggplot2::aes(fill = .data[[stat]]), width = r, height = r) +
    ggplot2::geom_sf(data = oz, inherit.aes = FALSE, fill = NA,
                     colour = "#6b6a66", linewidth = 0.2) +
    ggplot2::scale_fill_gradientn(colours = colours, transform = transform,
                                  limits = limits, oob = scales::squish, name = NULL) +
    ggplot2::coord_sf(crs = settings$crs, datum = NA, expand = FALSE,
                      xlim = c(-2100000, 2100000), ylim = c(-2200000, 1700000)) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA),
                   plot.title = ggplot2::element_text(size = 10),
                   strip.text = ggplot2::element_text(hjust = 0, size = 9),
                   legend.key.height = grid::unit(1.2, "cm"))
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(facet, ncol = 1)
  p
}

layer_colours <- list(
  mean = c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"),
  variance = c("#fde3d6", "#f6a57f", "#eb6834", "#b8461d", "#7a2c10")
)

#' Map of one trait: community mean and community SD, woody and non-woody.
plot_community_statistics <- function(layers, trait, unit = NULL) {
  d <- dplyr::filter(layers, trait_name == trait) |>
    dplyr::mutate(woodiness_group = factor(woodiness_group, c("woody", "non_woody"),
                                           c("Woody", "Non-woody")))
  unit_label <- if (is.null(unit) || is.na(unit)) "" else paste0(" (", unit, ")")
  patchwork::wrap_plots(
    map_cells(d, "community_mean", layer_colours$mean,
              paste0("Community mean", unit_label), "log10", facet = ~woodiness_group),
    map_cells(d, "community_sd", layer_colours$variance,
              "Community SD (natural-log units)", facet = ~woodiness_group),
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = paste0(trait, ": community statistics, location model"),
      subtitle = paste0(location_cell_km(), " km cells; from each cell's climate alone ",
                        "(co-occurring species' centres learned from GBIF)"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                             plot.subtitle = ggplot2::element_text(size = 9, colour = "#6b6a66")))
}

#' Save one map per trait (`figures/<trait>_community_statistics`); returns
#' the file paths.
write_community_maps <- function(layers, units, dir) {
  unlist(lapply(sort(unique(layers$trait_name)), function(tr) {
    p <- plot_community_statistics(layers, tr, units$unit[units$trait_name == tr])
    save_figure(p, file.path(dir, "figures", paste0(tr, "_community_statistics")),
                width = 9, height = 8)
  }), use.names = FALSE)
}

#' Units of the modelled traits.
trait_units <- function(observed_means, traits) {
  observed_means |>
    dplyr::filter(type == "numerical", trait_name %in% traits) |>
    dplyr::distinct(trait_name, unit)
}
