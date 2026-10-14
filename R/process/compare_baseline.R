# Compare a refit of the modelled means with a saved baseline
# (`docs/baselines/<label>/`, written before the data basis changed): fixed
# effects, record and species counts, and the species centres themselves.

#' Compare TMB fits with a baseline snapshot.
#'
#' @param fits list of `trait_climate_fit` (the `tmb_fits` target).
#' @param unified the unified record table the fits came from.
#' @param baseline_dir directory holding `fixed_effects.csv`,
#'   `trait_counts.csv` and `species_centres.csv.gz`.
#' @return list of tibbles: `fixed` (per trait × term), `counts` (per trait),
#'   `centres` (per-trait summary of species-level shifts) and `pairs` (the
#'   species present in both, for plotting).
compare_with_baseline <- function(fits, unified,
                                  baseline_dir = "docs/baselines/2026-06") {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  names(fits) <- vapply(fits, \(f) f$trait, character(1))
  read_base <- function(f) readr::read_csv(file.path(baseline_dir, f),
                                           show_col_types = FALSE)

  fixed <- dplyr::full_join(
    read_base("fixed_effects.csv"),
    purrr::imap_dfr(fits, \(f, tr) dplyr::mutate(f$fixed, trait = tr, .before = 1)),
    by = c("trait", "term"), suffix = c("_baseline", "_new")
  ) |>
    dplyr::transmute(
      trait, term,
      estimate_baseline, std.error_baseline, estimate_new, std.error_new,
      shift = estimate_new - estimate_baseline,
      shift_in_se = shift / std.error_baseline
    )

  counts_new <- unified |>
    dplyr::group_by(trait_name) |>
    dplyr::summarise(n_records_unified = dplyr::n(),
                     n_species_unified = dplyr::n_distinct(taxon_name)) |>
    dplyr::left_join(
      purrr::imap_dfr(fits, \(f, tr) tibble::tibble(
        trait_name = tr, n_records_fitted = nrow(f$data),
        n_species_fitted = dplyr::n_distinct(f$data$taxon_name))),
      by = "trait_name")
  counts <- dplyr::full_join(
    dplyr::select(read_base("trait_counts.csv"), -dplyr::any_of("n_species_observed_mean")),
    counts_new, by = "trait_name", suffix = c("_baseline", "_new")
  ) |>
    dplyr::select(trait_name, dplyr::starts_with("n_records_fitted"),
                  dplyr::starts_with("n_species_fitted"),
                  dplyr::starts_with("n_records_unified"),
                  dplyr::starts_with("n_species_unified"))

  centres_new <- purrr::imap_dfr(fits, \(f, tr) f$centres |>
    dplyr::transmute(trait_name = tr, taxon_name, n, log_estimate,
                     log_conf.low, log_conf.high))
  centres_base <- read_base("species_centres.csv.gz") |>
    dplyr::select(trait_name, taxon_name, n, log_estimate, log_conf.low, log_conf.high)
  pairs <- dplyr::inner_join(centres_base, centres_new,
                             by = c("trait_name", "taxon_name"),
                             suffix = c("_baseline", "_new")) |>
    dplyr::mutate(
      diff_log = log_estimate_new - log_estimate_baseline,
      # shift beyond the baseline interval half-width on its own side
      outside_baseline_ci = log_estimate_new < log_conf.low_baseline |
        log_estimate_new > log_conf.high_baseline
    )

  centres <- dplyr::full_join(
    dplyr::count(centres_base, trait_name, name = "n_species_baseline"),
    dplyr::count(centres_new, trait_name, name = "n_species_new"),
    by = "trait_name"
  ) |>
    dplyr::left_join(
      pairs |>
        dplyr::group_by(trait_name) |>
        dplyr::summarise(
          n_species_both = dplyr::n(),
          cor_log = stats::cor(log_estimate_baseline, log_estimate_new),
          median_diff_log = stats::median(diff_log),
          median_abs_diff_log = stats::median(abs(diff_log)),
          q95_abs_diff_log = stats::quantile(abs(diff_log), 0.95, names = FALSE),
          pct_outside_baseline_ci = 100 * mean(outside_baseline_ci)
        ),
      by = "trait_name")

  list(fixed = fixed, counts = counts, centres = centres, pairs = pairs)
}

#' Write the comparison tables to `dir` and a figure of species centres.
#' @return the file paths.
write_baseline_comparison <- function(cmp, dir = "outputs/comparisons",
                                      label = "baseline_2026-06") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(dir, paste0(label, "_", c("fixed_effects", "counts",
                                               "centres"), ".csv"))
  readr::write_csv(cmp$fixed, paths[1])
  readr::write_csv(cmp$counts, paths[2])
  readr::write_csv(cmp$centres, paths[3])
  figs <- save_figure(plot_baseline_centres(cmp$pairs),
                      file.path(dir, paste0(label, "_centres")),
                      width = 10, height = 9)
  c(paths, unname(figs))
}

#' Species centres, baseline against refit, one panel per trait (log scale).
plot_baseline_centres <- function(pairs) {
  pairs |>
    ggplot2::ggplot(ggplot2::aes(log_estimate_baseline, log_estimate_new)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(alpha = 0.25, stroke = 0, size = 0.8, colour = "#2C6E8E") +
    ggplot2::facet_wrap(~trait_name, scales = "free") +
    ggplot2::labs(x = "Modelled mean, June 2026 baseline (log)",
                  y = "Modelled mean, development build (log)") +
    ggplot2::theme_classic()
}
