# Observed-vs-modelled comparison (#24): coverage, agreement and where the
# two species estimates differ, per trait; and the same comparison between the
# portal-method and the geometric observed means.
#
# Differences are on the log scale: diff_log = log(modelled) - log(observed),
# so 0.1 is about a 10% difference.

#' Pair modelled and geometric observed means, with drivers of their difference.
#'
#' @param fits list of `trait_climate_fit` (the `tmb_fits` target).
#' @param observed_means geometric observed means.
#' @param unified the record table the fits used, for the climate offset.
#' @return list: `pairs` (one row per taxon × trait with both), `summary`
#'   (per trait), `drivers` (per trait × driver bin).
compare_observed_modelled <- function(fits, observed_means, unified) {
  fits <- purrr::compact(fits)
  modelled <- purrr::map_dfr(fits, \(f) f$centres |>
    dplyr::transmute(taxon_name, trait_name = f$trait,
                     log_modelled = log_estimate, n_modelled = n))
  traits <- unique(modelled$trait_name)
  observed <- observed_means |>
    dplyr::filter(type == "numerical", trait_name %in% traits) |>
    dplyr::transmute(taxon_name, trait_name, log_observed = log_value,
                     n_observed = n_records)

  # How far a species' records sit from its climate centre: the mean offset
  # of its records (x_c) in each climate variable, scaled by the across-
  # species spread of the centres, combined as a distance.
  offset <- unified |>
    dplyr::group_by(trait_name) |>
    dplyr::mutate(sd_p = stats::sd(log_annual_precipitation_mn),
                  sd_t = stats::sd(average_temperature_mn)) |>
    dplyr::group_by(trait_name, taxon_name) |>
    dplyr::summarise(.groups = "drop",
                     climate_offset = sqrt(
                       (mean(log_annual_precipitation_c) / dplyr::first(sd_p))^2 +
                         (mean(average_temperature_c) / dplyr::first(sd_t))^2))

  pairs <- dplyr::inner_join(modelled, observed, by = c("taxon_name", "trait_name")) |>
    dplyr::left_join(offset, by = c("trait_name", "taxon_name")) |>
    dplyr::mutate(diff_log = log_modelled - log_observed)

  coverage <- dplyr::full_join(
    dplyr::count(observed, trait_name, name = "n_species_observed"),
    dplyr::count(modelled, trait_name, name = "n_species_modelled"),
    by = "trait_name")

  summary <- pairs |>
    dplyr::group_by(trait_name) |>
    dplyr::summarise(
      n_species_both = dplyr::n(),
      cor_log = stats::cor(log_modelled, log_observed),
      median_diff_log = stats::median(diff_log),
      q05_diff_log = stats::quantile(diff_log, 0.05, names = FALSE),
      q95_diff_log = stats::quantile(diff_log, 0.95, names = FALSE),
      median_abs_diff_log = stats::median(abs(diff_log)),
      pct_abs_diff_over_0.2 = 100 * mean(abs(diff_log) > 0.2),
      rho_absdiff_n = stats::cor(abs(diff_log), n_modelled, method = "spearman"),
      rho_absdiff_offset = stats::cor(abs(diff_log), climate_offset,
                                      method = "spearman", use = "complete.obs")
    ) |>
    dplyr::right_join(coverage, by = "trait_name") |>
    dplyr::relocate(n_species_observed, n_species_modelled, .after = trait_name)

  bin <- function(x, breaks, labels) cut(x, breaks, labels, right = FALSE)
  drivers <- dplyr::bind_rows(
    pairs |>
      dplyr::mutate(driver = "records in the model",
                    bin = bin(n_modelled, c(6, 10, 20, 50, Inf),
                              c("6-9", "10-19", "20-49", "50+"))),
    pairs |>
      dplyr::mutate(driver = "climate offset of records",
                    bin = bin(climate_offset, c(0, 0.1, 0.25, 0.5, Inf),
                              c("<0.1", "0.1-0.25", "0.25-0.5", "0.5+")))
  ) |>
    dplyr::filter(!is.na(bin)) |>
    dplyr::group_by(trait_name, driver, bin) |>
    dplyr::summarise(.groups = "drop", n_species = dplyr::n(),
                     median_abs_diff_log = stats::median(abs(diff_log)))

  list(pairs = pairs, summary = summary, drivers = drivers)
}

#' Pair the portal-method (7.0.0) and geometric (development build) observed
#' means for numeric traits.
#'
#' Two contrasts: `diff_log` = log(geometric value) - log(portal value_mean)
#' combines the change of method and of data; `diff_log_data` =
#' log(geometric arithmetic_mean) - log(portal value_mean) isolates the data
#' (and roll-up), since both are arithmetic.
compare_portal_geometric <- function(portal, geometric) {
  p <- portal |>
    dplyr::filter(type == "numerical", value_mean > 0) |>
    dplyr::transmute(taxon_name, trait_name, unit, portal_mean = value_mean)
  g <- geometric |>
    dplyr::filter(type == "numerical") |>
    dplyr::transmute(taxon_name, trait_name, unit, geometric = value,
                     geometric_arithmetic = arithmetic_mean)
  pairs <- dplyr::inner_join(p, g, by = c("taxon_name", "trait_name", "unit")) |>
    dplyr::mutate(diff_log = log(geometric) - log(portal_mean),
                  diff_log_data = log(geometric_arithmetic) - log(portal_mean))

  summary <- pairs |>
    dplyr::group_by(trait_name) |>
    dplyr::summarise(
      n_species_both = dplyr::n(),
      cor_log = stats::cor(log(geometric), log(portal_mean)),
      median_diff_log = stats::median(diff_log),
      q05_diff_log = stats::quantile(diff_log, 0.05, names = FALSE),
      q95_diff_log = stats::quantile(diff_log, 0.95, names = FALSE),
      median_diff_log_data = stats::median(diff_log_data),
      pct_data_unchanged = 100 * mean(abs(diff_log_data) < 1e-9)
    ) |>
    dplyr::right_join(
      dplyr::full_join(
        dplyr::count(p, trait_name, name = "n_species_portal"),
        dplyr::count(g, trait_name, name = "n_species_geometric"),
        by = "trait_name"),
      by = "trait_name") |>
    dplyr::relocate(n_species_portal, n_species_geometric, .after = trait_name) |>
    dplyr::arrange(trait_name)

  list(pairs = pairs, summary = summary)
}

#' Modelled against observed means, one panel per trait (log scales).
plot_observed_modelled <- function(pairs) {
  pairs |>
    ggplot2::ggplot(ggplot2::aes(exp(log_observed), exp(log_modelled))) +
    ggplot2::geom_abline(linetype = "dashed", colour = "#6b6a66") +
    ggplot2::geom_point(alpha = 0.3, stroke = 0, size = 0.9, colour = "#2a78d6") +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~trait_name, scales = "free") +
    ggplot2::labs(x = "Observed mean (geometric method)",
                  y = "Modelled mean (TMB, at the climate centre)") +
    ggplot2::theme_classic(base_size = 10)
}

#' Median absolute difference by driver bin, per trait.
plot_difference_drivers <- function(drivers) {
  drivers |>
    ggplot2::ggplot(ggplot2::aes(bin, median_abs_diff_log, group = trait_name)) +
    ggplot2::geom_line(colour = "#86b6ef", linewidth = 0.5) +
    ggplot2::geom_point(colour = "#2a78d6", size = 1.5) +
    ggplot2::facet_wrap(~driver, scales = "free_x") +
    ggplot2::labs(x = NULL, y = "Median |log(modelled) - log(observed)|",
                  subtitle = "One line per trait") +
    ggplot2::theme_classic(base_size = 10)
}

#' Portal against geometric observed means: median log ratio and 5-95% range
#' per trait.
plot_portal_geometric <- function(summary) {
  summary |>
    dplyr::filter(!is.na(median_diff_log)) |>
    dplyr::mutate(trait_name = stats::reorder(trait_name, median_diff_log)) |>
    ggplot2::ggplot(ggplot2::aes(median_diff_log, trait_name)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#6b6a66", linewidth = 0.3) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = q05_diff_log, xmax = q95_diff_log),
                            colour = "#86b6ef", linewidth = 0.6) +
    ggplot2::geom_point(colour = "#2a78d6", size = 1.6) +
    ggplot2::labs(x = "log(geometric observed mean) - log(portal arithmetic mean): median and 5-95%",
                  y = NULL) +
    ggplot2::theme_classic(base_size = 9)
}

#' Write the comparison tables and figures to `dir`; returns the paths.
write_means_comparison <- function(obs_mod, portal_geo, dir = "outputs/comparisons") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tables <- file.path(dir, c("observed_vs_modelled_summary.csv",
                             "observed_vs_modelled_drivers.csv",
                             "portal_vs_geometric_summary.csv"))
  readr::write_csv(obs_mod$summary, tables[1])
  readr::write_csv(obs_mod$drivers, tables[2])
  readr::write_csv(portal_geo$summary, tables[3])
  figs <- c(
    save_figure(plot_observed_modelled(obs_mod$pairs),
                file.path(dir, "observed_vs_modelled"), width = 10, height = 9),
    save_figure(plot_difference_drivers(obs_mod$drivers),
                file.path(dir, "observed_vs_modelled_drivers"), width = 8, height = 4),
    save_figure(plot_portal_geometric(portal_geo$summary),
                file.path(dir, "portal_vs_geometric"), width = 8, height = 9)
  )
  c(tables, unname(figs))
}
