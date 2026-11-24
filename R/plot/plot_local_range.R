# Diagnostic plot for the location model prototype (#41): for a set of
# cells, the trait records and species means there, against the modelled
# community mean and range, from the location model and from the pooled
# model. One panel per cell, woody and non-woody side by side.

#' Community SD against climate, observed and modelled, per growth form.
#'
#' Points: each cell's observed SD among the species recorded there (less
#' the sampling term), coloured by how many species were recorded. Lines
#' (loess over cells): the SD each model implies for those same species:
#' pooled, location model with a constant species variance, and location
#' model with the species variance by climate class along that axis; and, in
#' red, the observed SD in cells with at least 10 species.
#'
#' @param res one element of the prototype result, with `diagnostics` and
#'   `spread`.
plot_spread_by_climate <- function(res) {
  sigma2 <- .vc(res$varcomp_new, "residual", "sigma")
  obs <- res$diagnostics$location$species |>
    dplyr::group_by(location, woodiness_group) |>
    dplyr::summarise(.groups = "drop", n_species = dplyr::n(),
                     sd = sqrt(max(stats::var(y) - mean(sigma2 / n), 0))) |>
    dplyr::inner_join(dplyr::select(res$diagnostics$location$cells, location, woodiness_group,
                                    annual_precipitation, average_temperature),
                      by = c("location", "woodiness_group")) |>
    dplyr::mutate(species_recorded = cut(n_species, c(4, 9, 19, Inf), labels = c("5-9", "10-19", "20+")))
  model_cells <- function(cells, label) dplyr::transmute(
    cells, location, woodiness_group, annual_precipitation, average_temperature,
    sd = community_sd, model = label)
  base_lines <- dplyr::bind_rows(
    model_cells(res$diagnostics$pooled$cells, "pooled"),
    model_cells(res$diagnostics$location$cells, "location, constant spread"))
  form_lab <- ggplot2::labeller(woodiness_group = c(woody = "woody", non_woody = "non-woody"))

  one_axis <- function(axis, col, xlab, log_x) {
    lines <- dplyr::bind_rows(
      base_lines,
      model_cells(res$spread[[axis]]$cells, "location, spread by climate"),
      if (!is.null(res$spread$continuous))
        model_cells(res$spread$continuous$cells, "location, continuous spread"))
    p <- ggplot2::ggplot(obs, ggplot2::aes(.data[[col]], sd)) +
      ggplot2::geom_point(ggplot2::aes(colour = species_recorded), size = 0.9, alpha = 0.6) +
      ggplot2::geom_smooth(data = dplyr::filter(obs, n_species >= 10), method = "loess",
                           formula = y ~ x, se = FALSE, span = 0.9,
                           colour = "#d73027", linewidth = 1.1) +
      ggplot2::geom_smooth(data = lines, ggplot2::aes(linetype = model), colour = "grey15",
                           method = "loess", formula = y ~ x, se = FALSE, span = 0.9,
                           linewidth = 0.8) +
      ggplot2::facet_grid(woodiness_group ~ ., labeller = form_lab) +
      ggplot2::scale_colour_manual(values = c(`5-9` = "#fdae61", `10-19` = "#74add1",
                                              `20+` = "#313695"),
                                   name = "species recorded\nin the cell") +
      ggplot2::scale_linetype_manual(values = c(pooled = "dotted",
                                                `location, constant spread` = "solid",
                                                `location, spread by climate` = "longdash",
                                                `location, continuous spread` = "twodash"),
                                     name = "modelled\ncommunity SD") +
      ggplot2::labs(x = xlab, y = "Community SD among species (natural-log units)") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    if (log_x) p <- p + ggplot2::scale_x_log10()
    p
  }
  patchwork::wrap_plots(
    one_axis("precipitation", "annual_precipitation", "Annual precipitation of the cell (mm)", TRUE),
    one_axis("temperature", "average_temperature", "Mean annual temperature of the cell (°C)", FALSE),
    guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(res$trait, ": spread among species within a cell, along climate"),
      subtitle = paste("Points: observed SD in each cell (less sampling noise). Red: trend (loess) of",
                       "the observed SD in cells with 10+ species. Black lines: trend of the SD each",
                       "model implies for the same cells"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11),
                             plot.subtitle = ggplot2::element_text(size = 7)))
}

#' The species SD (between-species spread, location effects removed) as a
#' function of climate, per growth form: the continuous RTMB estimate with
#' its 90% interval (the other climate variable at the growth form's median
#' species centre), the binned glmmTMB estimates per climate class (points
#' at each class's median centre) and the constant-spread estimate (dashed).
#'
#' @param res one element of the prototype result, with `spread$continuous`.
plot_spread_curves <- function(res) {
  cont <- res$spread$continuous$curves
  classes <- purrr::imap(res$spread[c("precipitation", "temperature")], function(s, ax)
    dplyr::mutate(s$classes, axis = ax)) |> purrr::list_rbind()
  constant <- tibble::tibble(woodiness_group = res$diagnostics$location$cells$woodiness_group |> unique()) |>
    dplyr::mutate(species_sd = sqrt(vapply(woodiness_group, function(g)
      .vc_species(res$varcomp_new, g), numeric(1))))
  form_lab <- ggplot2::labeller(woodiness_group = c(woody = "woody", non_woody = "non-woody"))
  one_axis <- function(ax, xlab, log_x) {
    p <- ggplot2::ggplot(dplyr::filter(cont, axis == ax), ggplot2::aes(climate, species_sd)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey80") +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_hline(data = constant, ggplot2::aes(yintercept = species_sd), linetype = "22") +
      ggplot2::geom_point(data = dplyr::filter(classes, axis == ax), colour = "#d73027", size = 2) +
      ggplot2::facet_grid(woodiness_group ~ ., labeller = form_lab) +
      ggplot2::expand_limits(y = 0) +
      ggplot2::labs(x = xlab, y = "Species SD (natural-log units)") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    if (log_x) p <- p + ggplot2::scale_x_log10()
    p
  }
  patchwork::wrap_plots(
    one_axis("precipitation", "Annual precipitation of the species' climate centre (mm)", TRUE),
    one_axis("temperature", "Mean annual temperature of the species' climate centre (°C)", FALSE)) +
    patchwork::plot_annotation(
      title = paste0(res$trait, ": spread among species as a function of climate"),
      subtitle = paste("Line and band: continuous model (RTMB), 90% interval, other climate variable",
                       "at the median species centre. Red points: one variance per climate quartile",
                       "of species (glmmTMB). Dashed: constant spread"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11),
                             plot.subtitle = ggplot2::element_text(size = 7)))
}

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
  dp <- res$diagnostics$pooled
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
    `pooled model` = prep(dp$cells) |> dplyr::mutate(x = x + 0.32),
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
    ggplot2::scale_linetype_manual(values = c(`location model` = "solid", `pooled model` = "22"),
                                   name = "modelled\ncommunity\nmean and\n90% range") +
    ggplot2::scale_shape_manual(values = c(`location model` = 16, `pooled model` = 1),
                                name = "modelled\ncommunity\nmean and\n90% range") +
    ggplot2::labs(x = NULL, y = paste0(res$trait, " (", res$unit, ")"),
                  title = paste0(res$trait, ": records, species means and modelled community range"),
                  subtitle = paste0(res$cell_km, " km cells along the rainfall gradient. ",
                                    "Faint dots: records; solid dots: species means there; ",
                                    "bars: modelled range of species means (solid: location model, dashed: pooled)")) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 7))
}
