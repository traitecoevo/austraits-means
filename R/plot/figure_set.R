# Figure set for presenting the species means and community products (#25),
# written to outputs/figures/ as PNG and PDF via save_figure().
#
# Colours follow one fixed assignment across figures: the portal-method
# observed mean is blue, the geometric observed mean orange, the modelled
# mean aqua (the first three slots of the reference categorical palette, which
# stay distinguishable for colour-blind readers); other two-way contrasts use
# dark and light ink so they never reuse a product colour; magnitudes use a
# single-hue blue ramp.

product_colours <- c(
  "Observed, portal method (7.0.0)" = "#2a78d6",
  "Observed, geometric (development build)" = "#eb6834",
  "Modelled (TMB, at climate centre)" = "#1baf7a"
)

#' Species with a mean per trait and product (the modelled traits).
coverage_table <- function(portal, geometric, modelled) {
  traits <- unique(modelled$trait_name)
  count_species <- function(d, label) {
    d |>
      dplyr::filter(trait_name %in% traits) |>
      dplyr::distinct(trait_name, taxon_name) |>
      dplyr::count(trait_name, name = "n_species") |>
      dplyr::mutate(product = label)
  }
  dplyr::bind_rows(
    count_species(dplyr::filter(portal, type == "numerical"), names(product_colours)[1]),
    count_species(dplyr::filter(geometric, type == "numerical"), names(product_colours)[2]),
    count_species(modelled, names(product_colours)[3])
  ) |>
    dplyr::mutate(product = factor(product, names(product_colours)))
}

#' Coverage: species per trait for each product (dot plot, log axis).
plot_coverage <- function(coverage) {
  order <- coverage |>
    dplyr::filter(product == names(product_colours)[2]) |>
    dplyr::arrange(n_species) |>
    dplyr::pull(trait_name)
  coverage |>
    dplyr::mutate(trait_name = factor(trait_name, order)) |>
    ggplot2::ggplot(ggplot2::aes(n_species, trait_name, colour = product)) +
    ggplot2::geom_line(ggplot2::aes(group = trait_name), colour = "#d6d5d1",
                       linewidth = 0.6) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(values = product_colours, name = NULL) +
    ggplot2::scale_x_log10(labels = scales::label_comma()) +
    ggplot2::labs(x = "Species with a mean (log scale)", y = NULL,
                  title = "Coverage of the species means",
                  caption = "Leaf N and P per area are not portal core traits, so have no portal-method mean.") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
}

#' Between- and within-species climate slopes per trait, with 95% intervals.
plot_climate_responses <- function(fits) {
  fx <- purrr::map_dfr(purrr::compact(fits), \(f) dplyr::mutate(f$fixed, trait = f$trait)) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      variable = dplyr::if_else(grepl("precipitation", term),
                                "log annual precipitation (per log mm)",
                                "mean annual temperature (per degC)"),
      level = dplyr::if_else(grepl("_mn$", term), "Between species", "Within species")
    )
  order <- fx |>
    dplyr::filter(level == "Between species", grepl("precipitation", variable)) |>
    dplyr::arrange(estimate) |>
    dplyr::pull(trait)
  fx |>
    dplyr::mutate(trait = factor(trait, order)) |>
    ggplot2::ggplot(ggplot2::aes(estimate, trait, colour = level)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#6b6a66", linewidth = 0.3) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = conf.low, xmax = conf.high),
                            position = ggplot2::position_dodge(width = 0.5),
                            linewidth = 0.6) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), size = 1.8) +
    ggplot2::facet_wrap(~variable, scales = "free_x") +
    # dark and light ink, so these two series do not reuse the product colours
    ggplot2::scale_colour_manual(values = c("Between species" = "#2b2b28",
                                            "Within species" = "#9c9b96"),
                                 name = NULL) +
    ggplot2::labs(x = "Slope of log trait value (95% interval)", y = NULL,
                  title = "Trait responses to climate, between and within species") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
}

#' Community mean maps for a few traits, woody and non-woody (observed, GBIF
#' presence x species means).
plot_community_overview <- function(layers, traits = c("plant_height", "leaf_mass_per_area",
                                                       "leaf_area", "seed_dry_mass")) {
  panels <- lapply(traits, function(tr) {
    d <- dplyr::filter(layers, trait_name == tr, woodiness_group %in% c("woody", "non_woody")) |>
      dplyr::mutate(woodiness_group = factor(woodiness_group, c("woody", "non_woody"),
                                             c("Woody", "Non-woody")))
    map_cells(d, "community_mean", layer_colours$mean, tr, "log10") +
      ggplot2::facet_wrap(~woodiness_group, ncol = 2)
  })
  patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = "Community mean, observed: GBIF presence x observed species means (10 km grid)",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11)))
}

#' One-page overview of the products: coverage, observed vs modelled for
#' height, climate responses, and the height community mean.
plot_product_overview <- function(coverage, obs_mod_pairs, fits, layers) {
  height_pairs <- dplyr::filter(obs_mod_pairs, trait_name == "plant_height")
  p_compare <- height_pairs |>
    ggplot2::ggplot(ggplot2::aes(exp(log_observed), exp(log_modelled))) +
    ggplot2::geom_abline(linetype = "dashed", colour = "#6b6a66") +
    ggplot2::geom_point(alpha = 0.2, stroke = 0, size = 0.7, colour = "#1baf7a") +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Observed mean (m)", y = "Modelled mean (m)",
                  title = "Plant height: modelled vs observed") +
    ggplot2::theme_classic(base_size = 9)
  p_map <- plot_community_overview(layers, "plant_height")[[1]] +
    ggplot2::labs(title = "Plant height: community mean (m)")
  patchwork::wrap_plots(
    plot_coverage(coverage) + ggplot2::theme(text = ggplot2::element_text(size = 8)),
    p_compare,
    plot_climate_responses(fits) + ggplot2::theme(text = ggplot2::element_text(size = 8)),
    p_map,
    ncol = 2, heights = c(1, 1)
  ) +
    patchwork::plot_annotation(
      title = "AusTraits species means and community layers: first draft",
      tag_levels = "a",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 13)))
}

#' Write the figure set to `dir`; returns the file paths. The #24 comparison
#' figures and per-trait community maps are written by their own targets and
#' copied here so the set is in one place.
write_figure_set <- function(portal, geometric, modelled, fits, obs_mod, layers,
                             comparison_files, community_map_files,
                             dir = "outputs/figures") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  coverage <- coverage_table(portal, geometric, modelled)
  readr::write_csv(coverage, file.path(dir, "coverage_by_product.csv"))
  made <- c(
    save_figure(plot_coverage(coverage), file.path(dir, "coverage_by_product"),
                width = 7, height = 5),
    save_figure(plot_climate_responses(fits), file.path(dir, "climate_responses"),
                width = 9, height = 5),
    save_figure(plot_community_overview(layers), file.path(dir, "community_mean_overview"),
                width = 10, height = 9),
    save_figure(plot_product_overview(coverage, obs_mod$pairs, fits, layers),
                file.path(dir, "product_overview"), width = 12, height = 10)
  )
  copied <- c(comparison_files, community_map_files)
  copied <- copied[grepl("\\.(png|pdf)$", copied)]
  dest <- file.path(dir, "parts", basename(copied))
  dir.create(file.path(dir, "parts"), showWarnings = FALSE)
  file.copy(copied, dest, overwrite = TRUE)
  c(file.path(dir, "coverage_by_product.csv"), unname(made), dest)
}
