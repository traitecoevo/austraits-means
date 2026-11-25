# All-species community statistics, for comparison only (roadmap.md, 25
# November 2026): published community-trait maps, such as Andrew et al. 2025
# (Global Ecology and Biogeography), pool woody and non-woody species. Our
# statistics are per growth form; the all-species version mixes them in each
# cell in proportion to the two forms' species counts there (GBIF presence).
# The 24 November decision not to publish an all-species statistic stands:
# it needs each location's woody share, which comes from occurrence data.
#
# With weights w_g = n_g / n, the mixture of the two forms' (log-normal)
# distributions of species means has
#   mean     = sum_g w_g mu_g
#   variance = sum_g w_g (sd_g^2 + mu_g^2) - mean^2
# i.e. the forms' own spreads plus the spread between their means.

#' @param stats output of community_statistics_climate() (per cell and form).
#' @param species output of cell_species() (GBIF species per cell, with form).
#' @param level mass of the community range.
#' @return tibble per trait x cell for traits fitted for both forms:
#'   n_woody, n_non_woody, share_woody, community_mean, community_sd,
#'   range_lower, range_upper, and the cell coordinates.
community_all_species <- function(stats, species, level = 0.90,
                                  settings = community_layer_settings()) {
  counts <- species |>
    dplyr::count(cell_x, cell_y, woodiness_group) |>
    tidyr::pivot_wider(names_from = woodiness_group, values_from = n,
                       names_prefix = "n_", values_fill = 0L)
  both <- stats |>
    dplyr::group_by(trait_name) |>
    dplyr::filter(dplyr::n_distinct(woodiness_group) == 2) |>
    dplyr::ungroup()
  both |>
    dplyr::transmute(trait_name, woodiness_group, cell_x, cell_y,
                     mu = log(community_mean), v = community_sd^2) |>
    tidyr::pivot_wider(names_from = woodiness_group, values_from = c(mu, v)) |>
    dplyr::inner_join(counts, by = c("cell_x", "cell_y")) |>
    dplyr::filter(n_woody + n_non_woody > 0) |>
    dplyr::mutate(
      share_woody = n_woody / (n_woody + n_non_woody),
      mean_log = share_woody * mu_woody + (1 - share_woody) * mu_non_woody,
      var_log = share_woody * (v_woody + mu_woody^2) +
        (1 - share_woody) * (v_non_woody + mu_non_woody^2) - mean_log^2) |>
    (\(d) dplyr::bind_cols(
      dplyr::select(d, trait_name, cell_x, cell_y, n_woody, n_non_woody, share_woody),
      community_range(d$mean_log, d$var_log, level)))() |>
    add_cell_centres(settings) |>
    dplyr::arrange(trait_name, cell_y, cell_x)
}

#' Map the all-species statistics: community mean and fold range per trait.
plot_community_all_species <- function(all, units, traits = unique(all$trait_name)) {
  panels <- unlist(lapply(traits, function(tr) {
    d <- dplyr::filter(all, trait_name == tr) |>
      dplyr::mutate(range_fold = range_upper / range_lower)
    u <- units$unit[units$trait_name == tr]
    list(map_cells(d, "community_mean", layer_colours$mean,
                   paste0(tr, ": community mean (", u, ")"), "log10"),
         map_cells(d, "range_fold", layer_colours$variance,
                   paste0(tr, ": 90% range, fold span"), "log10"))
  }), recursive = FALSE)
  patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = "All species (woody and non-woody mixed by their GBIF species counts): for comparison only",
      subtitle = paste0("Location model, ", location_cell_km(),
                        " km cells; not a product (no all-species statistic is published)"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11),
                             plot.subtitle = ggplot2::element_text(size = 9, colour = "#6b6a66")))
}

#' Write the all-species table and map to `dir`; returns the file paths.
write_community_all_species <- function(all, units, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  p <- file.path(dir, "community_all_species.parquet")
  arrow::write_parquet(all, p)
  n <- dplyr::n_distinct(all$trait_name)
  c(p, save_figure(plot_community_all_species(all, units), file.path(dir, "community_all_species"),
                   width = 10, height = 4.5 * n))
}
