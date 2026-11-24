# Diagnostic plot for the location model prototype (#41): for a set of
# cells, the trait records and species means there, against the modelled
# community mean and range, from the location model and from the published
# model. One panel per cell, woody and non-woody side by side.

#' Cells to show: up to `n_cells`, evenly spread along the precipitation
#' gradient, preferring cells with both growth forms.
select_diagnostic_cells <- function(cells, n_cells = 12) {
  by_cell <- cells |>
    dplyr::group_by(location) |>
    dplyr::summarise(.groups = "drop", n_groups = dplyr::n(),
                     n_species = sum(n_species),
                     annual_precipitation = exp(mean(log(annual_precipitation))),
                     average_temperature = mean(average_temperature))
  pool <- if (sum(by_cell$n_groups == 2) >= n_cells) dplyr::filter(by_cell, n_groups == 2) else by_cell
  pool <- dplyr::arrange(pool, annual_precipitation)
  idx <- unique(round(seq(1, nrow(pool), length.out = min(n_cells, nrow(pool)))))
  pool[idx, ]
}

#' The diagnostic plot for one trait's prototype result.
#'
#' @param res one element of the prototype result (`run_location_prototype()`),
#'   with `diagnostics`.
#' @param cells optional output of `select_diagnostic_cells()`.
plot_community_range_cells <- function(res, n_cells = 12, cells = NULL, seed = 1) {
  dl <- res$diagnostics$location
  dp <- res$diagnostics$published
  if (is.null(cells)) cells <- select_diagnostic_cells(dl$cells, n_cells)
  lab <- cells |>
    dplyr::mutate(panel = sprintf("%.0f mm, %.1f °C\n%d species",
                                  annual_precipitation, average_temperature, n_species),
                  panel = factor(panel, panel)) |>
    dplyr::select(location, panel)
  xpos <- c(woody = 1, non_woody = 2)
  set.seed(seed)
  prep <- function(d) d |>
    dplyr::inner_join(lab, by = "location") |>
    dplyr::mutate(x = xpos[woodiness_group])
  rec <- prep(dl$records) |> dplyr::mutate(x = x - 0.12 + stats::runif(dplyr::n(), -0.1, 0.1))
  sp <- prep(dl$species) |> dplyr::mutate(x = x - 0.12 + stats::runif(dplyr::n(), -0.1, 0.1))
  rng <- dplyr::bind_rows(
    `location model` = prep(dl$cells) |> dplyr::mutate(x = x + 0.17),
    `published model` = prep(dp$cells) |> dplyr::mutate(x = x + 0.32),
    .id = "model")
  cols <- c(woody = "#1b7837", non_woody = "#b2182b")

  ggplot2::ggplot() +
    ggplot2::geom_point(data = rec, ggplot2::aes(x, exp(log_value), colour = woodiness_group),
                        size = 0.4, alpha = 0.2, shape = 16) +
    ggplot2::geom_point(data = sp, ggplot2::aes(x, exp(y), colour = woodiness_group, size = n),
                        alpha = 0.75, shape = 16) +
    ggplot2::geom_linerange(data = rng, ggplot2::aes(x, ymin = exp(range_lower_log),
                                                     ymax = exp(range_upper_log), linetype = model),
                            linewidth = 0.7) +
    ggplot2::geom_point(data = rng, ggplot2::aes(x, exp(community_mean_log), shape = model), size = 1.8) +
    ggplot2::facet_wrap(~panel, nrow = 2) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(breaks = xpos, labels = c("woody", "non-woody"),
                                limits = c(0.6, 2.5)) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    # log scale with a floor: one species with hundreds of records in a cell
    # would otherwise shrink every other species mean to nothing
    ggplot2::scale_size(range = c(0.9, 3.2), trans = "log10", name = "records of\nthe species") +
    ggplot2::scale_linetype_manual(values = c(`location model` = "solid", `published model` = "22"),
                                   name = "modelled\ncommunity\nmean and\n90% range") +
    ggplot2::scale_shape_manual(values = c(`location model` = 16, `published model` = 1),
                                name = "modelled\ncommunity\nmean and\n90% range") +
    ggplot2::labs(x = NULL, y = paste0(res$trait, " (", res$unit, ")"),
                  title = paste0(res$trait, ": records, species means and modelled community range"),
                  subtitle = paste0(res$cell_km, " km cells along the rainfall gradient. ",
                                    "Faint dots: records; solid dots: species means there; ",
                                    "bars: modelled range of species means (solid: location model, dashed: published)")) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 7))
}
