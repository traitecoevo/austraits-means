# Community layers (#28): community statistics of trait values per cell of a
# 10 x 10 km equal-area grid, for all species and for woody and non-woody
# species separately. Decisions (roadmap.md, 19 October 2026): presence only,
# every species counting equally (#18); the woody rule of #19; the grid of
# Andrew et al. 2025, at least 10 species per cell and group, GBIF points for
# presence (#20).
#
# Three products:
#   - observed_local:    from the trait records collected in the cell: each
#                        species' records there are averaged, then the
#                        species. No model, no occurrence data.
#   - observed_presence: from the species GBIF records in the cell, each
#                        given its observed (geometric) species mean.
#   - modelled:          the fitted trait-climate model's net average for the
#                        cell (community_modelled.R).
#
# Definitions (terminology.md): community mean = exp(mean of the species' log
# values); community variance = variance of those log values, between species
# only, in natural-log units.

#' Grid settings.
community_layer_settings <- function() {
  list(
    resolution = 10000, min_species = 10,
    crs = "+proj=laea +lat_0=-25.2744 +lon_0=133.7751 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
    # Australia's occurrences in this projection, snapped to 10 km
    extent = c(xmin = -2300000, xmax = 2200000, ymin = -2200000, ymax = 1800000)
  )
}

#' Empty raster on the layer grid.
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
#' @param species_names the names used by the species means and the lookup.
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

# Species values per cell, duplicated into the woody / non-woody group of
# each species and into "all" (unassigned species count only in "all").
.with_groups <- function(d, groups) {
  dplyr::bind_rows(
    d |>
      dplyr::inner_join(groups |>
                          dplyr::filter(woodiness_group != "unassigned") |>
                          dplyr::select(taxon_name, woodiness_group),
                        by = "taxon_name"),
    dplyr::mutate(d, woodiness_group = "all")
  )
}

# Community statistics from one log value per species x trait x cell x group.
.community_statistics <- function(values, settings) {
  values |>
    dplyr::group_by(trait_name, cell_x, cell_y, woodiness_group) |>
    dplyr::summarise(.groups = "drop",
                     n_species = dplyr::n(),
                     mean_log = mean(log_value),
                     var_log = stats::var(log_value)) |>
    dplyr::mutate(
      enough = n_species >= settings$min_species,
      community_mean = dplyr::if_else(enough, exp(mean_log), NA_real_),
      community_variance = dplyr::if_else(enough, var_log, NA_real_)
    ) |>
    dplyr::select(-mean_log, -var_log, -enough)
}

#' Observed community layers from the trait records collected in each cell.
#'
#' Each species' records in a cell are averaged on the log scale (a geometric
#' mean), then the species are averaged. Uses the field records of the
#' modelled traits (`records_all`: raw, mean and median values, and the AVH
#' height bounds), with infraspecific taxa rolled up to species.
#'
#' @param records `records_all`.
#' @param groups woody lookup (`build_woody_lookup()`).
#' @return long tibble: trait, cell, group, `n_species` (species sampled),
#'   `n_records`, `community_mean`, `community_variance`.
community_layers_local <- function(records, groups, settings = community_layer_settings()) {
  per_species <- records |>
    dplyr::select(trait_name, taxon_name, log_value, decimal_longitude, decimal_latitude) |>
    dplyr::bind_cols(lonlat_to_cell(records$decimal_longitude, records$decimal_latitude, settings)) |>
    dplyr::group_by(trait_name, cell_x, cell_y, taxon_name) |>
    dplyr::summarise(.groups = "drop", n_records = dplyr::n(), log_value = mean(log_value))
  values <- .with_groups(per_species, groups)
  n_rec <- values |>
    dplyr::group_by(trait_name, cell_x, cell_y, woodiness_group) |>
    dplyr::summarise(n_records = sum(n_records), .groups = "drop")
  .community_statistics(values, settings) |>
    dplyr::left_join(n_rec, by = c("trait_name", "cell_x", "cell_y", "woodiness_group")) |>
    dplyr::relocate(n_records, .after = n_species) |>
    add_cell_centres(settings) |>
    dplyr::arrange(trait_name, woodiness_group, cell_y, cell_x)
}

#' Observed community layers from GBIF presence and species means.
#'
#' @param presence output of `gbif_presence_cells()`.
#' @param name_map output of `gbif_to_species()`.
#' @param groups woody lookup.
#' @param means species means: taxon_name, trait_name, log_value (the
#'   geometric observed means of the modelled traits).
#' @return long tibble: trait, cell, group, `n_species_present` (GBIF species
#'   in the cell, those without an AusTraits name included), `n_species`
#'   (present with a mean), `coverage`, `community_mean`, `community_variance`.
community_layers_presence <- function(presence, name_map, groups, means,
                                      settings = community_layer_settings()) {
  present <- presence |>
    dplyr::left_join(name_map, by = "gbif_taxon_name") |>
    dplyr::mutate(taxon_name = dplyr::coalesce(taxon_name, gbif_taxon_name)) |>
    dplyr::distinct(taxon_name, cell_x, cell_y) |>
    .with_groups(groups)
  n_present <- dplyr::count(present, cell_x, cell_y, woodiness_group,
                            name = "n_species_present")
  stats <- present |>
    dplyr::inner_join(dplyr::select(means, taxon_name, trait_name, log_value),
                      by = "taxon_name", relationship = "many-to-many") |>
    .community_statistics(settings)

  tidyr::crossing(n_present, trait_name = sort(unique(means$trait_name))) |>
    dplyr::left_join(stats, by = c("trait_name", "cell_x", "cell_y", "woodiness_group")) |>
    dplyr::mutate(n_species = dplyr::coalesce(n_species, 0L),
                  coverage = n_species / n_species_present) |>
    dplyr::relocate(n_species_present, .before = n_species) |>
    dplyr::relocate(coverage, .after = n_species) |>
    add_cell_centres(settings) |>
    dplyr::arrange(trait_name, woodiness_group, cell_y, cell_x)
}

#' Geometric observed means of the modelled traits, as the presence layers'
#' species values.
observed_species_values <- function(observed_means, traits) {
  observed_means |>
    dplyr::filter(type == "numerical", trait_name %in% traits) |>
    dplyr::select(taxon_name, trait_name, unit, log_value)
}

#' Write a set of layers: a long parquet, one GeoTIFF per trait x group x
#' statistic, and a metadata file.
#' @param statistics columns written as GeoTIFFs.
#' @return the file paths.
write_community_layers <- function(layers, dir, metadata, statistics,
                                   settings = community_layer_settings()) {
  dir.create(file.path(dir, "geotiff"), recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, c("community_layers.parquet", "community_layers_metadata.json"))
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

#' A map panel of one statistic on the layer grid, with state outlines.
#' Cells with `grey` set are drawn grey (species present, statistic empty).
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
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                   strip.text = ggplot2::element_text(hjust = 0, size = 9),
                   legend.key.height = grid::unit(1.2, "cm"))
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(facet, ncol = 1)
  p
}

layer_colours <- list(
  mean = c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"),
  variance = c("#fde3d6", "#f6a57f", "#eb6834", "#b8461d", "#7a2c10")
)

#' Quick-look map of one trait from an observed product: community mean and
#' variance, woody and non-woody. Grey: species present but fewer than the
#' minimum with a value.
plot_community_layers <- function(layers, trait, unit = NULL, product = "") {
  d <- dplyr::filter(layers, trait_name == trait,
                     woodiness_group %in% c("woody", "non_woody")) |>
    dplyr::mutate(woodiness_group = factor(woodiness_group, c("woody", "non_woody"),
                                           c("Woody", "Non-woody")))
  unit_label <- if (is.null(unit) || is.na(unit)) "" else paste0(" (", unit, ")")
  patchwork::wrap_plots(
    map_cells(d, "community_mean", layer_colours$mean,
              paste0("Community mean", unit_label), "log10", facet = ~woodiness_group),
    map_cells(d, "community_variance", layer_colours$variance,
              "Community variance (ln units)", facet = ~woodiness_group),
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = paste0(trait, ": ", product),
      subtitle = "10 km equal-area grid; grey: fewer than 10 species with a value",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                             plot.subtitle = ggplot2::element_text(size = 9, colour = "#6b6a66")))
}

#' Save one quick-look map per trait, named `<trait>_<key>`; returns the file
#' paths.
write_community_maps <- function(layers, units, dir, product, key) {
  unlist(lapply(sort(unique(layers$trait_name)), function(tr) {
    p <- plot_community_layers(layers, tr, units$unit[units$trait_name == tr], product)
    save_figure(p, file.path(dir, "figures", paste0(tr, "_", key)), width = 9, height = 8)
  }), use.names = FALSE)
}

#' Units of the modelled traits.
trait_units <- function(observed_means, traits) {
  observed_means |>
    dplyr::filter(type == "numerical", trait_name %in% traits) |>
    dplyr::distinct(trait_name, unit)
}

#' Metadata common to every community product.
community_metadata <- function(data_basis, product, description, extra = list(),
                               groups = paste("all species, and woody / non-woody by the rule of #19",
                                              "(woody: woodiness 'woody' or 'woody_like_stem';",
                                              "unassigned species only in 'all')"),
                               settings = community_layer_settings()) {
  c(list(
    product = paste("AusTraits community layers:", product),
    description = description,
    weighting = "presence only: every species counts equally (#18)",
    groups = groups,
    spatial_unit = paste0("10 x 10 km cells of an equal-area grid (", settings$crs,
                          "); x and y give cell centres in metres, longitude and latitude in degrees"),
    min_species = paste0("statistics are reported where at least ", settings$min_species,
                         " species in the group have a value; counts are always given"),
    definitions = list(
      community_mean = "exp(mean of the species' log values); a geometric mean, in trait units",
      community_variance = "variance of the species' log values; between species only, natural-log units")
  ), extra, data_basis, repo_provenance())
}

community_local_metadata <- function(data_basis) {
  community_metadata(
    data_basis, "observed, local records",
    paste("For each cell, each species' trait records collected in the cell are averaged",
          "on the log scale, then the species are averaged. No model and no occurrence",
          "data: only cells where the trait was measured have values."),
    extra = list(
      records = paste("field records of the modelled traits: value types raw, mean and",
                      "median, plus the AVH height bounds; infraspecific taxa rolled up to species"),
      columns = list(n_species = "species with a record of the trait in the cell and group",
                     n_records = "records behind those species' values")))
}

community_presence_metadata <- function(data_basis) {
  community_metadata(
    data_basis, "observed, GBIF presence x species means",
    paste("For each cell, the species with GBIF occurrences there, each given its",
          "geometric observed mean (one value per species, wherever it was measured)."),
    extra = list(
      occurrences = "cleaned GBIF occurrences (CC0 and CC BY), names aligned to APC",
      species_means = "geometric observed means of the development build (estimate_type observed_mean, method geometric)",
      columns = list(n_species_present = "GBIF species in the cell and group",
                     n_species = "of those, species with an observed mean of the trait",
                     coverage = "n_species / n_species_present")))
}
