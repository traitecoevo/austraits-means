# Community mean and community variance layers, first draft (#23).
#
# Provisional choices, listed in docs/community-layers-design.md and to be
# settled by #18, #19 and #20:
#   - presence only, every species present counting equally (#18), so the
#     statistics are a community mean and variance, never CWM / CWV;
#   - woody / non-woody from the provisional lookup (#19, #22);
#   - a regular 0.25° grid in EPSG:4326 (#20);
#   - a statistic only where at least 10 species in the group have a mean (#20);
#   - presence from the cleaned GBIF occurrences (#20);
#   - species means: variants are built for review (#20): `combined` (the
#     modelled mean where there is one, otherwise the geometric observed
#     mean) and `modelled_centre` (modelled means only); the model surface
#     (model_surface.R) is a further, model-based variant.
#
# Definitions (terminology.md): community mean = exp(mean of the log species
# means); community variance = variance of the log species means, between
# species only, in natural-log units.

#' Settings of the first-draft layers (see the header for their status).
community_layer_settings <- function() {
  list(resolution = 0.25, min_species = 10, crs = "EPSG:4326",
       extent = c(xmin = 112, xmax = 154, ymin = -44, ymax = -9))
}

#' Combined species means: modelled where available, otherwise observed.
#'
#' @param fits list of `trait_climate_fit` (the `tmb_fits` target).
#' @param observed_means geometric observed means.
#' @return one row per taxon × trait: value, log_value, unit, estimate_type,
#'   method.
combine_species_means <- function(fits, observed_means) {
  fits <- purrr::compact(fits)
  modelled <- purrr::map_dfr(fits, \(f) f$centres |>
    dplyr::transmute(taxon_name, trait_name = f$trait, value = estimate,
                     log_value = log_estimate, estimate_type = "modelled_mean",
                     method = f$backend))
  traits <- unique(modelled$trait_name)
  observed <- observed_means |>
    dplyr::filter(type == "numerical", trait_name %in% traits) |>
    dplyr::transmute(taxon_name, trait_name, value, log_value,
                     estimate_type = "observed_mean", method = "geometric")
  units <- observed_means |>
    dplyr::filter(type == "numerical", trait_name %in% traits) |>
    dplyr::distinct(trait_name, unit)
  stopifnot(!anyDuplicated(units$trait_name))

  dplyr::bind_rows(
    modelled,
    dplyr::anti_join(observed, modelled, by = c("taxon_name", "trait_name"))
  ) |>
    dplyr::left_join(units, by = "trait_name") |>
    dplyr::relocate(unit, estimate_type, method, .after = trait_name) |>
    dplyr::arrange(trait_name, taxon_name)
}

#' The species means a layer variant uses.
#'
#' @param means output of `combine_species_means()`.
#' @param input "combined" (all rows) or "modelled_centre" (modelled means,
#'   at each species' climate centre, only).
select_species_means <- function(means, input = c("combined", "modelled_centre")) {
  input <- match.arg(input)
  if (input == "modelled_centre") {
    dplyr::filter(means, estimate_type == "modelled_mean")
  } else {
    means
  }
}

species_means_input_description <- c(
  combined = "combined: the modelled mean (TMB) where there is one, otherwise the geometric observed mean",
  modelled_centre = "modelled means (TMB, at each species' climate centre) only; species without a modelled mean are counted as present but have no mean"
)

#' Species present per grid cell, from the GBIF occurrences.
#'
#' @return tibble(gbif_taxon_name, cell_x, cell_y); cell (x, y) spans
#'   [x, x + 1) × resolution in longitude and likewise in latitude.
gbif_presence_cells <- function(gbif_path, resolution = 0.25) {
  arrow::open_dataset(gbif_path) |>
    dplyr::select(taxon_name, decimal_longitude, decimal_latitude) |>
    dplyr::filter(!is.na(taxon_name), !is.na(decimal_longitude),
                  !is.na(decimal_latitude)) |>
    dplyr::mutate(cell_x = floor(decimal_longitude / resolution),
                  cell_y = floor(decimal_latitude / resolution)) |>
    dplyr::distinct(gbif_taxon_name = taxon_name, cell_x, cell_y) |>
    dplyr::collect()
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

#' Community statistics per trait × cell × group.
#'
#' @param presence output of `gbif_presence_cells()`.
#' @param name_map output of `gbif_to_species()`.
#' @param groups woody lookup (`build_woody_lookup()`); unassigned taxa are
#'   left out of the woody and non-woody groups but counted in `all`.
#' @param means output of `combine_species_means()`.
#' @param settings `community_layer_settings()`.
#' @return long tibble: trait_name, cell_x, cell_y, longitude, latitude (cell
#'   centre), woodiness_group, n_species_present, n_species, coverage,
#'   community_mean, community_variance.
compute_community_layers <- function(presence, name_map, groups, means,
                                     settings = community_layer_settings()) {
  res <- settings$resolution
  mapped <- presence |>
    dplyr::inner_join(name_map, by = "gbif_taxon_name") |>
    dplyr::distinct(taxon_name, cell_x, cell_y)
  present <- dplyr::bind_rows(
    mapped |>
      dplyr::inner_join(
        groups |>
          dplyr::filter(woodiness_group != "unassigned") |>
          dplyr::select(taxon_name, woodiness_group),
        by = "taxon_name"),
    dplyr::mutate(mapped, woodiness_group = "all")
  )

  n_present <- dplyr::count(present, cell_x, cell_y, woodiness_group,
                            name = "n_species_present")

  stats <- present |>
    dplyr::inner_join(dplyr::select(means, taxon_name, trait_name, log_value),
                      by = "taxon_name", relationship = "many-to-many") |>
    dplyr::group_by(trait_name, cell_x, cell_y, woodiness_group) |>
    dplyr::summarise(.groups = "drop",
                     n_species = dplyr::n(),
                     mean_log = mean(log_value),
                     var_log = stats::var(log_value))

  tidyr::crossing(n_present, trait_name = sort(unique(means$trait_name))) |>
    dplyr::left_join(stats, by = c("trait_name", "cell_x", "cell_y",
                                   "woodiness_group")) |>
    dplyr::mutate(
      n_species = dplyr::coalesce(n_species, 0L),
      coverage = n_species / n_species_present,
      enough = n_species >= settings$min_species,
      community_mean = dplyr::if_else(enough, exp(mean_log), NA_real_),
      community_variance = dplyr::if_else(enough, var_log, NA_real_),
      longitude = (cell_x + 0.5) * res,
      latitude = (cell_y + 0.5) * res
    ) |>
    dplyr::select(trait_name, cell_x, cell_y, longitude, latitude,
                  woodiness_group, n_species_present, n_species, coverage,
                  community_mean, community_variance) |>
    dplyr::arrange(trait_name, woodiness_group, cell_y, cell_x)
}

#' Write the layers: a long parquet, one GeoTIFF per trait × group ×
#' statistic, and a metadata file.
#' @return the file paths.
write_community_layers <- function(layers, dir, metadata,
                                   settings = community_layer_settings()) {
  dir.create(file.path(dir, "geotiff"), recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, c("community_layers.parquet",
                            "community_layers_metadata.json"))
  arrow::write_parquet(layers, paths[1])
  jsonlite::write_json(c(metadata, list(n_rows = nrow(layers))), paths[2],
                       auto_unbox = TRUE, pretty = TRUE)

  e <- settings$extent
  template <- terra::rast(xmin = e[["xmin"]], xmax = e[["xmax"]],
                          ymin = e[["ymin"]], ymax = e[["ymax"]],
                          resolution = settings$resolution, crs = settings$crs)
  statistics <- c("community_mean", "community_variance", "n_species", "coverage")
  combos <- dplyr::distinct(layers, trait_name, woodiness_group)
  tifs <- purrr::pmap_chr(combos, function(trait_name, woodiness_group) {
    d <- layers[layers$trait_name == trait_name &
                  layers$woodiness_group == woodiness_group, ]
    cells <- terra::cellFromXY(template, cbind(d$longitude, d$latitude))
    vapply(statistics, function(st) {
      r <- template
      v <- rep(NA_real_, terra::ncell(r))
      v[cells] <- as.numeric(d[[st]])
      terra::values(r) <- v
      names(r) <- st
      f <- file.path(dir, "geotiff", paste0(trait_name, "_", woodiness_group,
                                            "_", st, ".tif"))
      terra::writeRaster(r, f, overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
      f
    }, character(1)) |> paste(collapse = "\n")
  })
  c(paths, unlist(strsplit(tifs, "\n", fixed = TRUE)))
}

#' Quick-look map of one trait: community mean and variance, woody and
#' non-woody. Grey cells have species present but fewer than the minimum with
#' a mean.
plot_community_layers <- function(layers, trait, unit = NULL, input = "combined") {
  d <- dplyr::filter(layers, trait_name == trait,
                     woodiness_group %in% c("woody", "non_woody")) |>
    dplyr::mutate(woodiness_group = factor(woodiness_group,
                                           c("woody", "non_woody"),
                                           c("Woody", "Non-woody")))
  oz <- ozmaps::ozmap_states
  panel <- function(stat, colours, title, trans = "identity") {
    ggplot2::ggplot(d, ggplot2::aes(longitude, latitude)) +
      ggplot2::geom_tile(data = dplyr::filter(d, is.na(.data[[stat]])),
                         fill = "#e3e2de", width = 0.25, height = 0.25) +
      ggplot2::geom_tile(data = dplyr::filter(d, !is.na(.data[[stat]])),
                         ggplot2::aes(fill = .data[[stat]]),
                         width = 0.25, height = 0.25) +
      ggplot2::geom_sf(data = oz, inherit.aes = FALSE, fill = NA,
                       colour = "#6b6a66", linewidth = 0.2) +
      ggplot2::scale_fill_gradientn(colours = colours, transform = trans,
                                    name = NULL) +
      ggplot2::facet_wrap(~woodiness_group, ncol = 1) +
      ggplot2::coord_sf(xlim = c(112, 154), ylim = c(-44, -9.5), expand = FALSE) +
      ggplot2::labs(title = title, x = NULL, y = NULL) +
      ggplot2::theme_void(base_size = 10) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                     strip.text = ggplot2::element_text(hjust = 0, size = 9),
                     legend.key.height = grid::unit(1.2, "cm"))
  }
  blue <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
  orange <- c("#fde3d6", "#f6a57f", "#eb6834", "#b8461d", "#7a2c10")
  unit_label <- if (is.null(unit) || is.na(unit)) "" else paste0(" (", unit, ")")
  patchwork::wrap_plots(
    panel("community_mean", blue, paste0("Community mean", unit_label), "log10"),
    panel("community_variance", orange, "Community variance (ln units)"),
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = trait,
      subtitle = paste0("First draft, provisional choices; species means: ", input,
                        "; grey: fewer than 10 species with a mean"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                             plot.subtitle = ggplot2::element_text(size = 9, colour = "#6b6a66")))
}

#' Save one quick-look map per trait; returns the file paths.
write_community_maps <- function(layers, means, dir, input = "combined") {
  units <- dplyr::distinct(means, trait_name, unit)
  unlist(lapply(sort(unique(layers$trait_name)), function(tr) {
    p <- plot_community_layers(layers, tr, units$unit[units$trait_name == tr],
                               input = input)
    save_figure(p, file.path(dir, "figures", paste0(tr, "_community")),
                width = 9, height = 8)
  }), use.names = FALSE)
}

#' Metadata written with the layers: what they are and their provisional basis.
community_layers_metadata <- function(data_basis, input = "combined",
                                      settings = community_layer_settings()) {
  c(list(
    product = "AusTraits community layers, observed (first draft)",
    status = "provisional: built on choices still open in issues #18, #19 and #20",
    statistics = list(
      community_mean = "exp(mean of the log species means of the species present); a geometric mean, in trait units",
      community_variance = "variance of the log species means of the species present; between species only, natural-log units",
      n_species = "species present in the cell and group that have a species mean",
      n_species_present = "species present in the cell and group",
      coverage = "n_species / n_species_present"),
    weighting = "presence only: every species present counts equally (provisional, #18)",
    groups = "woody / non-woody from the provisional lookup (#19, #22); unassigned taxa left out",
    spatial_unit = paste0(settings$resolution, " degree grid, ", settings$crs,
                          " (provisional, #20); longitude and latitude give cell centres"),
    min_species = paste0(settings$min_species, " species with a mean per cell and group; statistics are empty below that (provisional, #20)"),
    occurrence_source = "cleaned GBIF occurrences, points (provisional, #20)",
    species_means_input = paste0(species_means_input_description[[input]],
                                 " (one of two variants built; the choice is #20)")
  ), data_basis, repo_provenance())
}
