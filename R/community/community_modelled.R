# Modelled community layer (#28): the fitted trait-climate model's net
# average trait value for a cell, averaging over between- and within-species
# effects. No species means are used.
#
# A record at climate x from a species centred at m has expected log trait
#   b0 + b_between * m + b_within * (x - m)
# (species random effects average to zero over species). Averaging over the
# species present in a cell gives
#   b0 + b_within * x + (b_between - b_within) * mean(m),
# so the cell value needs only the cell's climate x and the mean climate
# centre of the species found there (GBIF presence and the species medoids),
# for each climate variable. Species found in the cell need no trait data.

#' Mean climate of each grid cell (WorldClim averaged onto the layer grid,
#' cells touching Australian land only).
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

#' Mean climate centre of the species present in each cell and group.
#'
#' @param presence output of `gbif_presence_cells()` (GBIF names x cells).
#' @param medoids_path medoid parquet, keyed by GBIF name.
#' @param name_map,groups to place each GBIF species in woody / non-woody;
#'   species without a group count only in "all".
#' @return tibble(cell_x, cell_y, woodiness_group, n_species,
#'   centre_log_annual_precipitation, centre_average_temperature).
cell_species_centres <- function(presence, medoids_path, name_map, groups) {
  centres <- arrow::read_parquet(medoids_path) |>
    dplyr::filter(!is.na(annual_precipitation), annual_precipitation > 0) |>
    dplyr::transmute(gbif_taxon_name = taxon_name,
                     log_p = log(annual_precipitation), t = average_temperature)
  grp <- name_map |>
    dplyr::inner_join(dplyr::filter(groups, woodiness_group != "unassigned") |>
                        dplyr::select(taxon_name, woodiness_group), by = "taxon_name") |>
    dplyr::distinct(gbif_taxon_name, woodiness_group)
  present <- dplyr::inner_join(presence, centres, by = "gbif_taxon_name")
  dplyr::bind_rows(
    dplyr::inner_join(present, grp, by = "gbif_taxon_name"),
    dplyr::mutate(present, woodiness_group = "all")
  ) |>
    dplyr::group_by(cell_x, cell_y, woodiness_group) |>
    dplyr::summarise(.groups = "drop",
                     n_species = dplyr::n(),
                     centre_log_annual_precipitation = mean(log_p),
                     centre_average_temperature = mean(t))
}

#' The model's net average per trait x cell x group.
#'
#' @param fits list of TMB `trait_climate_fit` (`tmb_fits`).
#' @param climate output of `cell_climate()`.
#' @param centres output of `cell_species_centres()`.
#' @param groups which groups to predict for (the per-trait models are fitted
#'   on all species together, so "all" only until #19's woody split is fitted).
#' @param level interval mass (delta method on the log scale, fixed-effect
#'   uncertainty only).
#' @return long tibble: trait, cell, group, `n_species`, cell climate, mean
#'   species centre, `community_mean` with `lower`/`upper`, `extrapolated`.
community_layers_modelled <- function(fits, climate, centres, groups = "all",
                                      level = 0.94,
                                      settings = community_layer_settings()) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  cells <- centres |>
    dplyr::filter(woodiness_group %in% groups) |>
    dplyr::inner_join(climate, by = c("cell_x", "cell_y"))

  purrr::map_dfr(purrr::compact(fits), function(f) {
    stopifnot(f$backend == "tmb")
    b <- glmmTMB::fixef(f$fit)$cond
    V <- stats::vcov(f$fit)$cond[names(b), names(b)]
    X <- matrix(0, nrow(cells), length(b), dimnames = list(NULL, names(b)))
    X[, "(Intercept)"] <- 1
    X[, "log_annual_precipitation_mn"] <- cells$centre_log_annual_precipitation
    X[, "average_temperature_mn"] <- cells$centre_average_temperature
    X[, "log_annual_precipitation_c"] <- cells$log_annual_precipitation - cells$centre_log_annual_precipitation
    X[, "average_temperature_c"] <- cells$average_temperature - cells$centre_average_temperature
    eta <- as.numeric(X %*% b)
    se <- sqrt(rowSums((X %*% V) * X))

    # outside the fit's range: mean centre beyond the species centres fitted,
    # or cell-minus-centre beyond the within-species deviations fitted
    rng <- function(v) range(f$data[[v]])
    outside <- function(v, r) v < r[1] | v > r[2]
    extrapolated <-
      outside(cells$centre_log_annual_precipitation, rng("log_annual_precipitation_mn")) |
      outside(cells$centre_average_temperature, rng("average_temperature_mn")) |
      outside(X[, "log_annual_precipitation_c"], rng("log_annual_precipitation_c")) |
      outside(X[, "average_temperature_c"], rng("average_temperature_c"))

    enough <- cells$n_species >= settings$min_species
    cells |>
      dplyr::mutate(
        trait_name = f$trait, .before = 1,
        community_mean = dplyr::if_else(enough, exp(eta), NA_real_),
        lower = dplyr::if_else(enough, exp(eta - z * se), NA_real_),
        upper = dplyr::if_else(enough, exp(eta + z * se), NA_real_),
        extrapolated = extrapolated
      )
  }) |>
    add_cell_centres(settings) |>
    dplyr::arrange(trait_name, woodiness_group, cell_y, cell_x)
}

community_modelled_metadata <- function(data_basis) {
  community_metadata(
    data_basis, "modelled",
    paste("For each cell, the fitted trait-climate model's net average trait value over",
          "the species present: exp(b0 + b_within * x + (b_between - b_within) * m),",
          "where x is the cell's climate and m the mean climate centre (medoid) of the",
          "GBIF species recorded there, for log annual precipitation and mean annual",
          "temperature. Averages over between- and within-species effects; uses no",
          "species means."),
    groups = "all species (the per-trait models are fitted on all species together)",
    extra = list(
      model = "log_value ~ P_mn + T_mn + P_c + T_c + (1 | taxon) + (0 + P_c | taxon) + (0 + T_c | taxon), glmmTMB",
      climate = "WorldClim 2.1 (bio1, bio12) averaged onto the grid",
      interval = "lower, upper: 94% interval, delta method on the log scale, fixed-effect uncertainty only",
      columns = list(n_species = "GBIF species in the cell with a climate centre",
                     extrapolated = paste("TRUE where the mean centre or the cell-minus-centre",
                                          "difference lies outside the range the model was fitted on"))))
}

#' Map of the modelled layer beside the observed presence layer (all species)
#' on one colour scale. Extrapolated cells are shown grey.
plot_community_modelled <- function(modelled, observed, trait, unit = NULL) {
  m <- dplyr::filter(modelled, trait_name == trait, woodiness_group == "all")
  o <- dplyr::filter(observed, trait_name == trait, woodiness_group == "all")
  limits <- range(c(m$community_mean[!m$extrapolated], o$community_mean), na.rm = TRUE)
  unit_label <- if (is.null(unit) || is.na(unit)) "" else paste0(" (", unit, ")")
  patchwork::wrap_plots(
    map_cells(dplyr::mutate(m, community_mean = dplyr::if_else(extrapolated, NA_real_, community_mean)),
              "community_mean", layer_colours$mean, "Modelled: net average from the fitted model",
              "log10", limits = limits),
    map_cells(o, "community_mean", layer_colours$mean,
              "Observed: species present (GBIF) x observed species means",
              "log10", limits = limits),
    ncol = 2, guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = paste0(trait, ", all species, community mean", unit_label),
      subtitle = "10 km equal-area grid; grey: fewer than 10 species, or (modelled) outside the range the model was fitted on",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                             plot.subtitle = ggplot2::element_text(size = 9, colour = "#6b6a66")))
}

#' Save one modelled-vs-observed map per trait; returns the file paths.
write_modelled_maps <- function(modelled, observed, units, dir) {
  unlist(lapply(sort(unique(modelled$trait_name)), function(tr) {
    p <- plot_community_modelled(modelled, observed, tr, units$unit[units$trait_name == tr])
    save_figure(p, file.path(dir, "figures", paste0(tr, "_modelled")), width = 10, height = 5)
  }), use.names = FALSE)
}
