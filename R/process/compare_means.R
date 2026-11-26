# Observed-vs-modelled comparison (#24): coverage, agreement and where the
# two species estimates differ, per trait.
#
# Differences are on the log scale: diff_log = log(modelled) - log(observed),
# so 0.1 is about a 10% difference.

#' Pair modelled and geometric observed means, with drivers of their difference.
#'
#' @param fits list of location-model fits (the `location_fits` target).
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

#' Modelled against observed means, one panel per trait (log scales).
plot_observed_modelled <- function(pairs) {
  pairs |>
    ggplot2::ggplot(ggplot2::aes(exp(log_observed), exp(log_modelled))) +
    ggplot2::geom_abline(linetype = "dashed", colour = "#6b6a66") +
    ggplot2::geom_point(alpha = 0.3, stroke = 0, size = 0.9, colour = "#2a78d6") +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~trait_name, scales = "free") +
    ggplot2::labs(x = "Observed mean (geometric method)",
                  y = "Modelled mean (location model, at the climate centre)") +
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

#' Write the comparison tables and figures to `dir`; returns the paths.
write_means_comparison <- function(obs_mod, dir = "outputs/comparisons") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tables <- file.path(dir, c("observed_vs_modelled_summary.csv",
                             "observed_vs_modelled_drivers.csv"))
  readr::write_csv(obs_mod$summary, tables[1])
  readr::write_csv(obs_mod$drivers, tables[2])
  figs <- c(
    save_figure(plot_observed_modelled(obs_mod$pairs),
                file.path(dir, "observed_vs_modelled"), width = 10, height = 9),
    save_figure(plot_difference_drivers(obs_mod$drivers),
                file.path(dir, "observed_vs_modelled_drivers"), width = 8, height = 4)
  )
  c(tables, unname(figs))
}
