# Medoid reasonableness checks (data-quality gate; refactor_plan.md §3.2a).
#
# Two distinct questions, two functions:
#
# 1. check_medoid(): EXTRAPOLATION RISK relative to the trait records we centre
#    on the medoid. Because x_c = obs - medoid, if the medoid lies outside the
#    span of a species' trait-record climate then every x_c shares a sign and the
#    centre prediction is an extrapolation. Robust (range-based); the per-species
#    trait sample is often tiny, so SD-based z-scores are unstable and only
#    reported when n is adequate.
#
# 2. check_medoid_vs_gbif(): is the medoid a SENSIBLE CENTRE of the species'
#    actual distribution? Compares the medoid against the species' GBIF
#    occurrence cloud (large, well-sampled) — a small z here confirms the medoid
#    sits in the dense core, not at an edge.

#' Extrapolation check: medoid vs the trait records' climate span.
#'
#' @param data record-level data from `add_medoid_centring()`.
#' @param tol fraction of the trait-data range tolerated beyond min/max.
#' @param min_n_z minimum records before an SD-based z-score is reported.
#' @return tibble, one row per species, with `medoid_ok = in_span`.
check_medoid <- function(data, tol = 0.10, min_n_z = 5) {
  span_ok <- function(centre, x, tol) {
    lo <- min(x); hi <- max(x); rng <- hi - lo
    centre >= lo - tol * rng & centre <= hi + tol * rng
  }
  z_dist <- function(centre, x, n) {
    s <- stats::sd(x)
    if (n < min_n_z || is.na(s) || s == 0) return(NA_real_)
    abs(centre - mean(x)) / s
  }

  data |>
    dplyr::group_by(taxon_name) |>
    dplyr::summarise(
      .groups = "drop",
      n = dplyr::n(),
      P_centre = dplyr::first(log_annual_precipitation_mn),
      P_med = stats::median(log_annual_precipitation),
      P_in_span = span_ok(dplyr::first(log_annual_precipitation_mn), log_annual_precipitation, tol),
      P_z = z_dist(dplyr::first(log_annual_precipitation_mn), log_annual_precipitation, dplyr::n()),
      T_centre = dplyr::first(average_temperature_mn),
      T_med = stats::median(average_temperature),
      T_in_span = span_ok(dplyr::first(average_temperature_mn), average_temperature, tol),
      T_z = z_dist(dplyr::first(average_temperature_mn), average_temperature, dplyr::n())
    ) |>
    dplyr::mutate(
      in_span = P_in_span & T_in_span,
      medoid_ok = in_span
    )
}

#' Short text summary of a `check_medoid()` result.
summarise_medoid_check <- function(check) {
  n <- nrow(check)
  message(sprintf(
    "Medoid extrapolation check: %d species | within trait span %d (%.1f%%) | outside span %d",
    n, sum(check$in_span), 100 * mean(check$in_span), sum(!check$in_span)
  ))
  invisible(check)
}

#' Reasonableness vs the GBIF occurrence cloud (the medoid's own distribution).
#'
#' For a sample of species, loads GBIF occurrences and measures how central the
#' medoid is within that cloud: distance from the cloud mean in SDs (`P_z`, `T_z`)
#' and whether it lies within the interquartile box. A reasonable medoid has
#' small z and `in_iqr = TRUE`.
#'
#' @param medoid_lookup output of `build_medoid_lookup()` (centre per taxon).
#' @param gbif_path GBIF parquet.
#' @param n_sample number of species to sample (cloud loading is the cost).
#' @param min_n minimum GBIF occurrences for a species to be checkable.
#' @param seed sampling seed.
check_medoid_vs_gbif <- function(medoid_lookup = build_medoid_lookup(),
                                 gbif_path = "data-processed/data_gbif_cleaned.parquet",
                                 n_sample = 300, min_n = 30, seed = 1) {
  set.seed(seed)
  taxa <- sample(unique(medoid_lookup$taxon_name),
                 min(n_sample, length(unique(medoid_lookup$taxon_name))))

  gbif <- arrow::open_dataset(gbif_path) |>
    dplyr::filter(taxon_name %in% taxa) |>
    dplyr::select(taxon_name, annual_precipitation, average_temperature) |>
    dplyr::filter(!is.na(annual_precipitation), annual_precipitation > 0,
                  !is.na(average_temperature)) |>
    dplyr::collect() |>
    dplyr::mutate(log_annual_precipitation = log(annual_precipitation))

  centres <- medoid_lookup |> dplyr::filter(taxon_name %in% taxa)

  gbif |>
    dplyr::group_by(taxon_name) |>
    dplyr::summarise(
      .groups = "drop",
      n_gbif = dplyr::n(),
      P_mean = mean(log_annual_precipitation), P_sd = stats::sd(log_annual_precipitation),
      P_q1 = stats::quantile(log_annual_precipitation, .25),
      P_q3 = stats::quantile(log_annual_precipitation, .75),
      T_mean = mean(average_temperature), T_sd = stats::sd(average_temperature),
      T_q1 = stats::quantile(average_temperature, .25),
      T_q3 = stats::quantile(average_temperature, .75)
    ) |>
    dplyr::filter(n_gbif >= min_n) |>
    dplyr::left_join(centres, by = "taxon_name") |>
    dplyr::mutate(
      P_z = abs(log_annual_precipitation_centre - P_mean) / P_sd,
      T_z = abs(average_temperature_centre - T_mean) / T_sd,
      in_iqr = log_annual_precipitation_centre >= P_q1 & log_annual_precipitation_centre <= P_q3 &
               average_temperature_centre >= T_q1 & average_temperature_centre <= T_q3,
      central = P_z <= 1 & T_z <= 1
    )
}

#' Summary of `check_medoid_vs_gbif()`.
summarise_medoid_vs_gbif <- function(check) {
  message(sprintf(
    "Medoid vs GBIF cloud: %d species | median P_z %.2f, T_z %.2f | within IQR box %.1f%% | central (z<=1 both) %.1f%%",
    nrow(check), stats::median(check$P_z), stats::median(check$T_z),
    100 * mean(check$in_iqr), 100 * mean(check$central)
  ))
  invisible(check)
}

#' Diagnostic scatter of medoid-to-cloud distances (z), saved png+pdf if `file` set.
plot_medoid_check <- function(check, file = NULL, threshold = 1) {
  p <- check |>
    ggplot2::ggplot(ggplot2::aes(P_z, T_z)) +
    ggplot2::geom_point(alpha = 0.4, stroke = 0, colour = "#2C6E8E") +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed") +
    ggplot2::labs(
      x = "Medoid distance from cloud centre — precip (SDs)",
      y = "Medoid distance from cloud centre — temperature (SDs)",
      title = "Medoid centrality within the occurrence cloud"
    ) +
    ggplot2::theme_classic()

  if (!is.null(file)) save_figure(p, file, width = 7, height = 6)
  p
}
