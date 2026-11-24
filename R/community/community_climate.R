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
# the range of the fitted species' SDs as in community_statistics(). The
# species-spread term is scaled by a sorting factor k <= 1
# (community_sorting.R). Intervals come from draws of the model's fixed
# parameters.
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
    # mgcv::predict.gam explicitly: a fit read back from the store dispatches
    # to predict.glm when mgcv is not loaded
    pr <- lapply(s$fits, function(f) as.numeric(mgcv::predict.gam(f, newdata = climate)))
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

# The parameters the community statistics use, from a vector of the fit's
# fixed parameters (its estimates `par_fixed`, or a draw around them).
.community_params <- function(m, v) {
  terms <- m$spread$terms
  a_names <- paste0(m$groups, ":(Intercept)")
  b_names <- if (length(terms)) unlist(lapply(m$groups, function(g) paste0(g, ":", terms)))
  list(fixed = stats::setNames(v[names(v) == "beta"], names(m$fixed)),
       theta = c(stats::setNames(v[names(v) == "a"], a_names),
                 if (length(terms)) stats::setNames(v[names(v) == "b"], b_names)),
       sd_p = exp(v[["log_sd_P"]]), sd_t = exp(v[["log_sd_T"]]), sd_w = exp(v[["log_sd_w"]]))
}

# log SD of the species spread at standardised centres (p, t), growth form g.
.log_sd_u <- function(m, theta, g, p, t) {
  theta[[paste0(g, ":(Intercept)")]] +
    Reduce(`+`, lapply(m$spread$terms, function(k)
      theta[[paste0(g, ":", k)]] * .spread_term_cols(p, t, k)[, 1]), 0)
}

# Community mean (log) and variance of growth form g at the cells of `d`
# (predicted centres), for parameters `pr`; k scales the species spread.
.community_climate_moments <- function(m, pr, g, d, fitted_z, k, nodes) {
  b <- pr$fixed[.group_coef_names(m, g)]
  names(b) <- c("b0", "bP", "bT", "wP", "wT")
  s <- m$centre_scaling
  xP <- d$log_annual_precipitation
  xT <- d$average_temperature
  dP <- b[["bP"]] - b[["wP"]]
  dT <- b[["bT"]] - b[["wT"]]
  vPP <- d$sd_p^2
  vTT <- d$sd_t^2
  vPT <- d$cor_pt * d$sd_p * d$sd_t
  mean_log <- b[["b0"]] + b[["wP"]] * xP + b[["wT"]] * xT + dP * d$mean_p + dT * d$mean_t
  var_between <- dP^2 * vPP + dT^2 * vTT + 2 * dP * dT * vPT
  var_slopes <- pr$sd_p^2 * ((xP - d$mean_p)^2 + vPP) + pr$sd_t^2 * ((xT - d$mean_t)^2 + vTT)
  # E[sd_u(m)^2] over m ~ N(mean, S): quadrature on the Cholesky factor,
  # each node's log SD clamped to the range of the fitted species'
  sd_range <- range(.log_sd_u(m, pr$theta, g, fitted_z$p, fitted_z$t))
  l11 <- d$sd_p
  l21 <- vPT / pmax(l11, 1e-12)
  l22 <- sqrt(pmax(vTT - l21^2, 0))
  var_u <- Reduce(`+`, lapply(seq_along(nodes$w), function(i) {
    mp <- d$mean_p + l11 * nodes$z1[i]
    mt <- d$mean_t + l21 * nodes$z1[i] + l22 * nodes$z2[i]
    log_sd <- .log_sd_u(m, pr$theta, g, (mp - s$p_mean) / s$p_sd, (mt - s$t_mean) / s$t_sd)
    nodes$w[i] * exp(2 * pmin(pmax(log_sd, sd_range[1]), sd_range[2]))
  }))
  list(mean_log = mean_log, variance = var_between + k * var_u + var_slopes + pr$sd_w^2)
}

#' Community mean and range per trait x cell x growth form, from the
#' location model and the cell's climate alone (see the file header), with
#' intervals for parameter uncertainty.
#'
#' @param fits list of fit_location_target() results.
#' @param centres output of `predict_centre_surfaces()`.
#' @param sorting optional `estimate_sorting_factors()` result: the species
#'   spread is scaled by its `k` (community_sorting.R).
#' @param level mass of the community range.
#' @param n_draws draws of the fit's fixed parameters from their estimated
#'   normal distribution (the inverse Hessian), for the intervals; 0 for
#'   none.
#' @param ci_level mass of those intervals.
#' @param spatial `centres` are grid cells (add their coordinates); FALSE for
#'   points of a climate grid (community_surface()).
#' @return long tibble: trait_name, woodiness_group, cell, cell climate,
#'   `extrapolated` (outside the trait's record climate for the growth form),
#'   `centres_extrapolated` (outside the cells the centre surfaces were fitted
#'   to), sorting factor `k`, community_mean, community_sd, range_lower,
#'   range_upper, and the `ci_level` intervals community_mean_lower / _upper
#'   and community_sd_lower / _upper. The intervals carry the uncertainty of
#'   the model's parameters (those at the boundary held fixed), not of the
#'   centre surfaces or of k.
community_statistics_climate <- function(fits, centres, sorting = NULL, level = 0.90,
                                         n_draws = 200, ci_level = 0.94, spatial = TRUE,
                                         settings = community_layer_settings()) {
  gh <- .gauss_hermite_normal(7)
  grid <- expand.grid(i = seq_along(gh$z), j = seq_along(gh$z))
  nodes <- list(z1 = gh$z[grid$i], z2 = gh$z[grid$j], w = gh$w[grid$i] * gh$w[grid$j])
  qs <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  purrr::map(purrr::compact(fits), function(m) {
    s <- m$centre_scaling
    draws <- NULL
    if (n_draws > 0 && !is.null(m$cov_fixed)) {
      # variance parameters at the boundary (zero rows; fit_spread_rtmb()) are
      # held at their estimate
      V <- m$cov_fixed
      free <- diag(V) > 0
      R <- tryCatch(chol(V[free, free, drop = FALSE]), error = function(e) NULL)
      if (!is.null(R)) {
        set.seed(1)
        z <- matrix(stats::rnorm(n_draws * sum(free)), n_draws)
        draws <- matrix(m$par_fixed, n_draws, length(m$par_fixed), byrow = TRUE)
        draws[, free] <- sweep(z %*% R, 2, m$par_fixed[free], `+`)
        colnames(draws) <- names(m$par_fixed)
      }
    }
    fitted <- dplyr::distinct(m$data, taxon_name, woodiness_group,
                              log_annual_precipitation_mn, average_temperature_mn)
    purrr::map(m$groups, function(g) {
      d <- dplyr::filter(centres, woodiness_group == g)
      if (nrow(d) == 0) return(NULL)
      fg <- dplyr::filter(fitted, woodiness_group == g)
      fitted_z <- list(p = (fg$log_annual_precipitation_mn - s$p_mean) / s$p_sd,
                       t = (fg$average_temperature_mn - s$t_mean) / s$t_sd)
      k <- .sorting_k(sorting, m$trait, g)
      est <- .community_climate_moments(m, .community_params(m, m$par_fixed), g, d,
                                        fitted_z, k, nodes)
      ci <- if (is.null(draws)) NULL else {
        sims <- lapply(seq_len(nrow(draws)), function(i) .community_climate_moments(
          m, .community_params(m, draws[i, ]), g, d, fitted_z, k, nodes))
        mean_q <- apply(do.call(cbind, lapply(sims, `[[`, "mean_log")), 1, stats::quantile, qs)
        sd_q <- apply(sqrt(do.call(cbind, lapply(sims, `[[`, "variance"))), 1, stats::quantile, qs)
        tibble::tibble(community_mean_lower = exp(mean_q[1, ]), community_mean_upper = exp(mean_q[2, ]),
                       community_sd_lower = sd_q[1, ], community_sd_upper = sd_q[2, ])
      }
      rec <- dplyr::filter(m$data, woodiness_group == g)
      outside <- function(v, r) v < min(r) | v > max(r)
      d |>
        dplyr::transmute(
          trait_name = m$trait, woodiness_group = g, cell_x, cell_y,
          annual_precipitation, log_annual_precipitation, average_temperature,
          extrapolated = outside(log_annual_precipitation, rec$log_annual_precipitation) |
            outside(average_temperature, rec$average_temperature),
          centres_extrapolated, k = k) |>
        dplyr::bind_cols(community_range(est$mean_log, est$variance, level)) |>
        dplyr::bind_cols(if (is.null(ci)) tibble::tibble(
          community_mean_lower = NA_real_, community_mean_upper = NA_real_,
          community_sd_lower = NA_real_, community_sd_upper = NA_real_) else ci)
    }) |> purrr::list_rbind()
  }) |>
    purrr::list_rbind() |>
    (\(d) if (spatial) add_cell_centres(d, settings) else d)() |>
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
