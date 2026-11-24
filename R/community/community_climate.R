# Community statistics from climate alone (roadmap.md, Later work, next
# step 1): the product. community_statistics() (community_statistics.R)
# needs the species GBIF records in each cell; here the location model needs
# only the climate x, because what it uses of the species present is the
# distribution of their climate centres, and that is learned as a smooth
# function of climate from the GBIF cells.
#
# For growth form g at climate x, with the co-occurring species' centres m
# having mean mbar(x) and covariance S(x) (2 x 2: log precipitation,
# temperature), the location model (community_statistics.R, header) gives
#   community mean (log) = b0 + b_within . x + d . mbar,   d = b_between - b_within
#   community variance   = d' S d                                     between-species turnover
#                        + E_m[sd_u(m)^2]                             species spread
#                        + sd_P^2 E[(x_P - m_P)^2] + sd_T^2 E[(x_T - m_T)^2]   species slopes
#                        + sd_w^2                                     species x location
# with E[(x_k - m_k)^2] = (x_k - mbar_k)^2 + S_kk, and E_m[sd_u(m)^2] taken
# over m ~ N(mbar, S) by Gauss-Hermite quadrature, each node's SD kept within
# the range of the fitted species' SDs as in community_statistics().
#
# mbar(x) and S(x) are smooth surfaces (mgcv, a tensor product of log
# precipitation and temperature) fitted per growth form to the cells with at
# least `min_species` species: the mean centres, the log SDs of the centres
# and the Fisher-transformed correlation.

#' Per cell and growth form: the climate, the number of species and the mean,
#' SDs and correlation of their climate centres.
#'
#' @param species output of `cell_species()`.
#' @param climate output of `cell_climate()`.
cell_centre_summaries <- function(species, climate) {
  species |>
    dplyr::group_by(cell_x, cell_y, woodiness_group) |>
    dplyr::summarise(
      .groups = "drop", n_species = dplyr::n(),
      mean_p = mean(centre_log_annual_precipitation),
      mean_t = mean(centre_average_temperature),
      sd_p = stats::sd(centre_log_annual_precipitation),
      sd_t = stats::sd(centre_average_temperature),
      cor_pt = suppressWarnings(stats::cor(centre_log_annual_precipitation,
                                           centre_average_temperature))) |>
    dplyr::inner_join(climate, by = c("cell_x", "cell_y"))
}

#' Smooth surfaces of the centres' distribution against climate, per growth
#' form.
#'
#' @param summaries output of `cell_centre_summaries()`.
#' @param min_species cells used: at least this many species of the form.
#' @return list per growth form: the five mgcv fits (`mean_p`, `mean_t`,
#'   `log_sd_p`, `log_sd_t`, `z_cor`) and the climate range they were fitted
#'   on.
fit_centre_surfaces <- function(summaries, min_species = 20, k = 8) {
  responses <- c(mean_p = "mean_p", mean_t = "mean_t", log_sd_p = "log(sd_p)",
                 log_sd_t = "log(sd_t)", z_cor = "atanh(pmin(pmax(cor_pt, -0.99), 0.99))")
  groups <- sort(unique(summaries$woodiness_group))
  stats::setNames(lapply(groups, function(g) {
    d <- dplyr::filter(summaries, woodiness_group == g, n_species >= min_species)
    fits <- lapply(responses, function(r) mgcv::gam(
      stats::as.formula(paste(r, "~ te(log_annual_precipitation, average_temperature, k =", k, ")")),
      data = d, method = "REML"))
    list(fits = fits, n_cells = nrow(d),
         range_p = range(d$log_annual_precipitation), range_t = range(d$average_temperature))
  }), groups)
}

#' The centres' predicted mean and covariance at each cell's climate.
#' @return tibble(cell_x, cell_y, woodiness_group, cell climate, mean_p,
#'   mean_t, sd_p, sd_t, cor_pt, centres_extrapolated).
predict_centre_surfaces <- function(surfaces, climate) {
  purrr::imap(surfaces, function(s, g) {
    pr <- lapply(s$fits, function(f) as.numeric(stats::predict(f, newdata = climate)))
    climate |>
      dplyr::mutate(woodiness_group = g, .after = cell_y,
                    mean_p = pr$mean_p, mean_t = pr$mean_t,
                    sd_p = exp(pr$log_sd_p), sd_t = exp(pr$log_sd_t), cor_pt = tanh(pr$z_cor),
                    centres_extrapolated =
                      log_annual_precipitation < s$range_p[1] | log_annual_precipitation > s$range_p[2] |
                      average_temperature < s$range_t[1] | average_temperature > s$range_t[2])
  }) |> purrr::list_rbind()
}

# Nodes and weights for E[f(Z)], Z ~ N(0, 1) (Gauss-Hermite, Golub-Welsch).
.gauss_hermite_normal <- function(n = 7) {
  J <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1) / 2)
  J[cbind(1:(n - 1), 2:n)] <- off
  J[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(z = sqrt(2) * e$values, w = e$vectors[1, ]^2)
}

#' Community mean and range per trait x cell x growth form, from the
#' location model and the cell's climate alone (see the file header).
#'
#' @param fits list of fit_location_target() results.
#' @param centres output of `predict_centre_surfaces()`.
#' @param level mass of the community range.
#' @return long tibble in the schema of community_statistics(), with
#'   `n_species` replaced by the surfaces' `centres_extrapolated` flag; every
#'   cell gets statistics.
community_statistics_climate <- function(fits, centres, level = 0.90,
                                         settings = community_layer_settings()) {
  gh <- .gauss_hermite_normal(7)
  nodes <- expand.grid(i = seq_along(gh$z), j = seq_along(gh$z))
  z1 <- gh$z[nodes$i]
  z2 <- gh$z[nodes$j]
  w <- gh$w[nodes$i] * gh$w[nodes$j]
  purrr::map(purrr::compact(fits), function(m) {
    vc <- m$vc
    sP <- vc$sd[vc$grouping == "taxon_name" & vc$term == "log_annual_precipitation_c"]
    sT <- vc$sd[vc$grouping == "taxon_name" & vc$term == "average_temperature_c"]
    sw <- vc$sd[vc$grouping == "taxon_name:location"]
    s <- m$centre_scaling
    th <- m$spread$theta
    fitted_sp <- dplyr::distinct(m$data, taxon_name, woodiness_group)
    purrr::map(m$groups, function(g) {
      b <- m$fixed[.group_coef_names(m, g)]
      names(b) <- c("b0", "bP", "bT", "wP", "wT")
      d <- dplyr::filter(centres, woodiness_group == g)
      if (nrow(d) == 0) return(NULL)
      xP <- d$log_annual_precipitation
      xT <- d$average_temperature
      dP <- b[["bP"]] - b[["wP"]]
      dT <- b[["bT"]] - b[["wT"]]
      vPP <- d$sd_p^2
      vTT <- d$sd_t^2
      vPT <- d$cor_pt * d$sd_p * d$sd_t
      mean_log <- b[["b0"]] + b[["wP"]] * xP + b[["wT"]] * xT + dP * d$mean_p + dT * d$mean_t
      var_between <- dP^2 * vPP + dT^2 * vTT + 2 * dP * dT * vPT
      var_slopes <- sP^2 * ((xP - d$mean_p)^2 + vPP) + sT^2 * ((xT - d$mean_t)^2 + vTT)
      # E[sd_u(m)^2] over m ~ N(mean, S): quadrature on the Cholesky factor,
      # each node's log SD clamped to the fitted species' range
      fitted_log_sd <- 0.5 * log(m$species_variance$species_variance[
        m$species_variance$taxon_name %in% fitted_sp$taxon_name[fitted_sp$woodiness_group == g]])
      l11 <- d$sd_p
      l21 <- vPT / pmax(l11, 1e-12)
      l22 <- sqrt(pmax(vTT - l21^2, 0))
      var_u <- Reduce(`+`, lapply(seq_along(w), function(k) {
        mp <- d$mean_p + l11 * z1[k]
        mt <- d$mean_t + l21 * z1[k] + l22 * z2[k]
        p <- (mp - s$p_mean) / s$p_sd
        t <- (mt - s$t_mean) / s$t_sd
        log_sd <- th[[paste0(g, ":(Intercept)")]] +
          Reduce(`+`, lapply(m$spread$terms, function(term)
            th[[paste0(g, ":", term)]] * .spread_term_cols(p, t, term)[, 1]), 0)
        w[k] * exp(2 * pmin(pmax(log_sd, min(fitted_log_sd)), max(fitted_log_sd)))
      }))
      rec <- dplyr::filter(m$data, woodiness_group == g)
      outside <- function(v, r) v < min(r) | v > max(r)
      d |>
        dplyr::transmute(
          trait_name = m$trait, woodiness_group = g, cell_x, cell_y,
          annual_precipitation, log_annual_precipitation, average_temperature,
          extrapolated = outside(log_annual_precipitation, rec$log_annual_precipitation) |
            outside(average_temperature, rec$average_temperature),
          centres_extrapolated) |>
        dplyr::bind_cols(community_range(mean_log, var_between + var_u + var_slopes + sw^2,
                                         level))
    }) |> purrr::list_rbind()
  }) |>
    purrr::list_rbind() |>
    add_cell_centres(settings) |>
    dplyr::arrange(trait_name, woodiness_group, cell_y, cell_x)
}

#' The climate-alone statistics against those from each cell's own GBIF
#' species (community_statistics()), in cells with at least `min_species`
#' species: how much is lost by predicting the species present from climate.
#'
#' @return list(summary (per trait x growth form: cells, correlation and
#'   median absolute difference of the log community mean and of the
#'   community SD), pairs).
compare_community_climate_gbif <- function(climate_stats, gbif_stats, min_species = 20) {
  pairs <- dplyr::inner_join(
    dplyr::select(climate_stats, trait_name, woodiness_group, cell_x, cell_y,
                  mean_climate = community_mean, sd_climate = community_sd),
    dplyr::filter(gbif_stats, n_species >= min_species) |>
      dplyr::select(trait_name, woodiness_group, cell_x, cell_y, n_species,
                    mean_gbif = community_mean, sd_gbif = community_sd),
    by = c("trait_name", "woodiness_group", "cell_x", "cell_y"))
  summary <- pairs |>
    dplyr::group_by(trait_name, woodiness_group) |>
    dplyr::summarise(
      .groups = "drop", cells = dplyr::n(),
      cor_log_mean = stats::cor(log(mean_climate), log(mean_gbif)),
      median_abs_diff_log_mean = stats::median(abs(log(mean_climate / mean_gbif))),
      cor_sd = stats::cor(sd_climate, sd_gbif),
      median_abs_diff_sd = stats::median(abs(sd_climate - sd_gbif)),
      median_sd_gbif = stats::median(sd_gbif))
  list(summary = summary, pairs = pairs)
}

#' Write the comparison: the summary table and a figure of the cell pairs.
#' @return the file paths.
write_community_climate_check <- function(check, dir) {
  dir.create(file.path(dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  p <- file.path(dir, "climate_vs_gbif_species_summary.csv")
  readr::write_csv(check$summary, p)
  d <- check$pairs |>
    dplyr::mutate(woodiness_group = factor(woodiness_group, c("woody", "non_woody"),
                                           c("woody", "non-woody")))
  fig <- patchwork::wrap_plots(
    ggplot2::ggplot(d, ggplot2::aes(mean_gbif, mean_climate, colour = woodiness_group)) +
      ggplot2::geom_abline(linetype = "dashed", colour = "grey50") +
      ggplot2::geom_point(size = 0.3, alpha = 0.3) +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggplot2::facet_wrap(~trait_name, scales = "free", ncol = 3) +
      ggplot2::labs(x = "From the cell's GBIF species", y = "From climate alone",
                    colour = NULL, title = "Community mean") +
      ggplot2::theme_bw(base_size = 8),
    ggplot2::ggplot(d, ggplot2::aes(sd_gbif, sd_climate, colour = woodiness_group)) +
      ggplot2::geom_abline(linetype = "dashed", colour = "grey50") +
      ggplot2::geom_point(size = 0.3, alpha = 0.3) +
      ggplot2::facet_wrap(~trait_name, scales = "free", ncol = 3) +
      ggplot2::labs(x = "From the cell's GBIF species", y = "From climate alone",
                    colour = NULL, title = "Community SD (natural-log units)") +
      ggplot2::theme_bw(base_size = 8),
    ncol = 2, guides = "collect")
  c(p, save_figure(fig, file.path(dir, "figures", "climate_vs_gbif_species"),
                   width = 14, height = 7))
}
