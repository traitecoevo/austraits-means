# Forecasting from the location model's exported parameters (roadmap.md, 25
# November 2026): a trait value for a named species at any climate, for a
# species without trait data from its growth form and climate centre
# (forecast_unmeasured()), and the community mean and range of a growth form
# at any climate, each with uncertainty, from published tables and no model
# machinery.
#
# Species. In the location model a species' expected log trait value at a
# site of climate x, at an average location, is
#   mu_s(x) = c_s + beta_P,s (x_P - m_P,s) + beta_T,s (x_T - m_T,s)
# with c_s its modelled mean at its climate centre m_s and beta_.,s its total
# within-species slopes (the growth form's plus its own). The table holds
# c_s, beta_P,s, beta_T,s and their 3 x 3 covariance Sigma_s (including the
# uncertainty of the model's fixed parameters), so with a = (1, dP, dT)
#   var(mu_s(x)) = a' Sigma_s a.
# For the value at an actual new site add the location and species x
# location variances (sd_location^2 + sd_species_location^2); for a single
# new record also the residual variance (sigma^2): the site-variances table.
#
# Community. The community statistics (community_climate.R) depend on the
# location only through its climate, so they are tabulated on a regular grid
# of log annual precipitation x mean annual temperature, per trait and growth
# form, and a forecast is a bilinear interpolation on that grid.

#' The species-responses table: per trait x species, the three numbers that
#' give its expected value at any climate, their covariance, and the climate
#' range of its records.
#' @param fits list of fit_location_target() results.
species_response_table <- function(fits) {
  purrr::map(purrr::compact(fits), function(m) {
    rec <- m$data |>
      dplyr::group_by(taxon_name) |>
      dplyr::summarise(n_records = dplyr::n(),
                       record_precipitation_min = exp(min(log_annual_precipitation)),
                       record_precipitation_max = exp(max(log_annual_precipitation)),
                       record_temperature_min = min(average_temperature),
                       record_temperature_max = max(average_temperature))
    m$centres |>
      dplyr::transmute(
        taxon_name, trait_name = m$trait, unit = m$data$unit[1],
        woodiness_group = as.character(woodiness_group),
        centre_precipitation = exp(log_annual_precipitation_mn),
        centre_temperature = average_temperature_mn,
        log_value_centre = log_estimate,
        slope_precipitation, slope_temperature,
        var_centre = log_se^2, var_slope_precipitation, var_slope_temperature,
        cov_centre_slope_precipitation, cov_centre_slope_temperature,
        cov_slope_precipitation_temperature) |>
      dplyr::left_join(rec, by = "taxon_name")
  }) |> purrr::list_rbind()
}

# The five fixed coefficients of a growth form, in this order, with short
# names for the export.
.coef_short <- c("intercept", "precipitation_between", "temperature_between",
                 "precipitation_within", "temperature_within")

#' The growth-form table: per trait x growth form, what forecasting a species
#' without trait data needs (see forecast_unmeasured()): the five fixed
#' coefficients and their covariance, the species spread as a function of the
#' climate centre (and the range it is clamped to), the SDs of the species'
#' climate slopes, and the climate ranges of the fitted species' centres and
#' of the records.
growth_form_response_table <- function(fits) {
  purrr::map(purrr::compact(fits), function(m) {
    beta_i <- which(names(m$par_fixed) == "beta")
    Vb <- m$cov_fixed[beta_i, beta_i]
    dimnames(Vb) <- list(names(m$fixed), names(m$fixed))
    vc <- m$vc
    s <- m$centre_scaling
    th <- m$spread$theta
    fitted <- dplyr::distinct(m$data, taxon_name, woodiness_group,
                              log_annual_precipitation_mn, average_temperature_mn)
    purrr::map(m$groups, function(g) {
      nm <- .group_coef_names(m, g)
      b <- stats::setNames(unname(m$fixed[nm]), paste0("coef_", .coef_short))
      V <- Vb[nm, nm]
      ij <- which(upper.tri(V, diag = TRUE), arr.ind = TRUE)
      cv <- stats::setNames(V[ij], paste0("cov_", .coef_short[ij[, 1]], "__", .coef_short[ij[, 2]]))
      fg <- dplyr::filter(fitted, woodiness_group == g)
      sp <- .log_sd_u(m, th, g, (fg$log_annual_precipitation_mn - s$p_mean) / s$p_sd,
                      (fg$average_temperature_mn - s$t_mean) / s$t_sd)
      rec <- dplyr::filter(m$data, woodiness_group == g)
      tib <- function(x) tibble::as_tibble(as.list(x))
      dplyr::bind_cols(
        tibble::tibble(trait_name = m$trait, unit = m$data$unit[1], woodiness_group = g,
                       n_species = nrow(fg)),
        tib(b), tib(cv),
        tibble::tibble(
          spread_intercept = th[[paste0(g, ":(Intercept)")]],
          spread_precipitation = if (paste0(g, ":p") %in% names(th)) th[[paste0(g, ":p")]] else 0,
          spread_temperature = if (paste0(g, ":t") %in% names(th)) th[[paste0(g, ":t")]] else 0,
          spread_scale_precipitation_mean = s$p_mean, spread_scale_precipitation_sd = s$p_sd,
          spread_scale_temperature_mean = s$t_mean, spread_scale_temperature_sd = s$t_sd,
          spread_log_sd_min = min(sp), spread_log_sd_max = max(sp),
          sd_slope_precipitation = vc$sd[vc$grouping == "taxon_name" &
                                           vc$term == "log_annual_precipitation_c"],
          sd_slope_temperature = vc$sd[vc$grouping == "taxon_name" &
                                         vc$term == "average_temperature_c"],
          centre_precipitation_min = exp(min(fg$log_annual_precipitation_mn)),
          centre_precipitation_max = exp(max(fg$log_annual_precipitation_mn)),
          centre_temperature_min = min(fg$average_temperature_mn),
          centre_temperature_max = max(fg$average_temperature_mn),
          record_precipitation_min = exp(min(rec$log_annual_precipitation)),
          record_precipitation_max = exp(max(rec$log_annual_precipitation)),
          record_temperature_min = min(rec$average_temperature),
          record_temperature_max = max(rec$average_temperature)))
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind()
}

#' The species that can be forecast without trait data: every species with a
#' climate centre (the medoid of its GBIF occurrences) and a woody or
#' non-woody class.
species_centre_table <- function(medoid_lookup, woody_lookup) {
  medoid_lookup |>
    dplyr::inner_join(dplyr::filter(woody_lookup, woodiness_group %in% c("woody", "non_woody")) |>
                        dplyr::select(taxon_name, woodiness_group), by = "taxon_name") |>
    dplyr::transmute(taxon_name, woodiness_group,
                     centre_precipitation = annual_precipitation_centre,
                     centre_temperature = average_temperature_centre) |>
    dplyr::distinct(taxon_name, .keep_all = TRUE) |>
    dplyr::arrange(taxon_name)
}

#' Forecast the trait value of a species without trait data, from its growth
#' form and climate centre, at given climates, from the growth-form and
#' site-variances tables only.
#'
#' The species' expected log value at climate x, at an average location, is
#' the growth form's line b0 + b_between . m + b_within . (x - m), with m its
#' centre. Its variance adds to the line's (from the coefficients'
#' covariance) the spread of species about it: the species SD at its centre,
#' sd_u(m)^2 (clamped to the range of the fitted species'), and the species'
#' unknown climate slopes, sd_P^2 dP^2 + sd_T^2 dT^2. "site" and "record" add
#' the site variances as in forecast_species().
#'
#' @param groups,variances the exported tables (growth_form_response_table(),
#'   site_variance_table()).
#' @param centre_precipitation,centre_temperature the species' climate centre
#'   (mm, deg C); from species_centre_table() or the user's own.
#' @return tibble as forecast_species(), with `extrapolated` (the climate is
#'   outside the growth form's records) and `centre_extrapolated` (the centre
#'   is outside those of the fitted species).
forecast_unmeasured <- function(groups, variances, trait, woodiness_group, centre_precipitation,
                                centre_temperature, annual_precipitation, average_temperature,
                                type = c("mean", "site", "record"), level = 0.94) {
  type <- match.arg(type)
  g <- groups[groups$trait_name == trait & groups$woodiness_group == woodiness_group, ]
  if (nrow(g) != 1) stop("no growth-form row for ", trait, ", ", woodiness_group)
  mP <- log(centre_precipitation)
  mT <- centre_temperature
  dp <- log(annual_precipitation) - mP
  dt <- average_temperature - mT
  X <- cbind(1, mP, mT, dp, dt)
  b <- unlist(g[paste0("coef_", .coef_short)])
  V <- matrix(0, 5, 5)
  for (i in 1:5) for (j in i:5)
    V[i, j] <- V[j, i] <- g[[paste0("cov_", .coef_short[i], "__", .coef_short[j])]]
  est <- as.vector(X %*% b)
  log_sd <- g$spread_intercept +
    g$spread_precipitation * (mP - g$spread_scale_precipitation_mean) / g$spread_scale_precipitation_sd +
    g$spread_temperature * (mT - g$spread_scale_temperature_mean) / g$spread_scale_temperature_sd
  log_sd <- pmin(pmax(log_sd, g$spread_log_sd_min), g$spread_log_sd_max)
  v <- rowSums((X %*% V) * X) + exp(2 * log_sd) +
    g$sd_slope_precipitation^2 * dp^2 + g$sd_slope_temperature^2 * dt^2
  if (type != "mean") {
    sv <- variances[variances$trait_name == trait & variances$woodiness_group == woodiness_group, ]
    v <- v + sv$sd_location^2 + sv$sd_species_location^2 + if (type == "record") sv$sd_residual^2 else 0
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  se <- sqrt(v)
  tibble::tibble(annual_precipitation, average_temperature,
                 estimate = exp(est), lower = exp(est - z * se), upper = exp(est + z * se),
                 log_estimate = est, log_se = se,
                 extrapolated = annual_precipitation < g$record_precipitation_min |
                   annual_precipitation > g$record_precipitation_max |
                   average_temperature < g$record_temperature_min |
                   average_temperature > g$record_temperature_max,
                 centre_extrapolated = centre_precipitation < g$centre_precipitation_min |
                   centre_precipitation > g$centre_precipitation_max |
                   centre_temperature < g$centre_temperature_min |
                   centre_temperature > g$centre_temperature_max)
}

#' The site-variances table: per trait x growth form, the SDs to add to a
#' species forecast for a new site or a new record (natural-log units).
site_variance_table <- function(fits) {
  purrr::map(purrr::compact(fits), function(m) {
    vc <- m$vc
    tibble::tibble(
      trait_name = m$trait, woodiness_group = m$groups,
      sd_location = vc$sd[vc$grouping == "location"],
      sd_species_location = vc$sd[vc$grouping == "taxon_name:location"],
      sd_residual = vc$sd[vc$grouping == "residual"])
  }) |> purrr::list_rbind()
}

#' A regular climate grid for the community surface: log-spaced annual
#' precipitation and evenly spaced temperature over the range of the
#' continent's cells.
climate_grid <- function(climate, n = 60) {
  lp <- seq(min(climate$log_annual_precipitation), max(climate$log_annual_precipitation),
            length.out = n)
  tt <- seq(min(climate$average_temperature), max(climate$average_temperature), length.out = n)
  tidyr::expand_grid(log_annual_precipitation = lp, average_temperature = tt) |>
    dplyr::mutate(annual_precipitation = exp(log_annual_precipitation),
                  cell_x = dplyr::row_number(), cell_y = 0L)
}

#' The community surface: the community statistics (from climate alone) on a
#' climate grid, per trait x growth form, with `inside_australia`: whether the
#' grid point is within reach of a real cell's climate (some cell within one
#' grid step on both axes); forecasts outside it are extrapolations of the
#' co-occurrence surfaces too.
community_surface <- function(fits, centre_surfaces, climate, sorting = NULL, n = 60) {
  grid <- climate_grid(climate, n)
  step_p <- diff(range(grid$log_annual_precipitation)) / (n - 1)
  step_t <- diff(range(grid$average_temperature)) / (n - 1)
  inside <- vapply(seq_len(nrow(grid)), function(i)
    any(abs(climate$log_annual_precipitation - grid$log_annual_precipitation[i]) <= step_p &
          abs(climate$average_temperature - grid$average_temperature[i]) <= step_t), logical(1))
  community_statistics_climate(fits, predict_centre_surfaces(centre_surfaces, grid),
                               sorting = sorting, spatial = FALSE) |>
    dplyr::left_join(dplyr::transmute(grid, cell_x, inside_australia = inside), by = "cell_x") |>
    dplyr::select(-cell_x, -cell_y)
}

#' Forecast a named species' trait value at given climates, from the
#' species-responses and site-variances tables only.
#'
#' @param responses,variances the exported tables (species_response_table(),
#'   site_variance_table()).
#' @param taxon,trait the species and trait.
#' @param annual_precipitation,average_temperature the climates (mm, deg C),
#'   vectors of equal length.
#' @param type "mean": the species' expected value at that climate, at an
#'   average location; "site": its mean at a new site there (adds the location
#'   and species x location variances); "record": a single new record there
#'   (adds the residual too).
#' @param level interval mass.
#' @return tibble: the climate, `estimate` with `lower` / `upper` (trait units),
#'   `log_estimate`, `log_se`, and `extrapolated` (outside the climate range
#'   of the species' records).
forecast_species <- function(responses, variances, taxon, trait, annual_precipitation,
                             average_temperature, type = c("mean", "site", "record"),
                             level = 0.94) {
  type <- match.arg(type)
  r <- responses[responses$taxon_name == taxon & responses$trait_name == trait, ]
  if (nrow(r) != 1) stop("no single row for ", taxon, ", ", trait)
  dp <- log(annual_precipitation) - log(r$centre_precipitation)
  dt <- average_temperature - r$centre_temperature
  est <- r$log_value_centre + r$slope_precipitation * dp + r$slope_temperature * dt
  v <- r$var_centre + dp^2 * r$var_slope_precipitation + dt^2 * r$var_slope_temperature +
    2 * dp * r$cov_centre_slope_precipitation + 2 * dt * r$cov_centre_slope_temperature +
    2 * dp * dt * r$cov_slope_precipitation_temperature
  if (type != "mean") {
    sv <- variances[variances$trait_name == trait & variances$woodiness_group == r$woodiness_group, ]
    v <- v + sv$sd_location^2 + sv$sd_species_location^2 + if (type == "record") sv$sd_residual^2 else 0
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  se <- sqrt(v)
  tibble::tibble(annual_precipitation, average_temperature,
                 estimate = exp(est), lower = exp(est - z * se), upper = exp(est + z * se),
                 log_estimate = est, log_se = se,
                 extrapolated = annual_precipitation < r$record_precipitation_min |
                   annual_precipitation > r$record_precipitation_max |
                   average_temperature < r$record_temperature_min |
                   average_temperature > r$record_temperature_max)
}

#' Forecast the community statistics of a growth form at given climates,
#' from the community-surface table only: bilinear interpolation on its grid
#' (log precipitation x temperature), of the log community mean and range and
#' of the community SD and the interval ends.
#' @return tibble: the climate and the interpolated community_mean,
#'   community_sd, range_lower, range_upper, their interval ends, and the
#'   flags (TRUE if any surrounding grid point carries them).
forecast_community <- function(surface, trait, woodiness_group, annual_precipitation,
                               average_temperature) {
  s <- surface[surface$trait_name == trait & surface$woodiness_group == woodiness_group, ]
  if (nrow(s) == 0) stop("no surface for ", trait, ", ", woodiness_group)
  lp_axis <- sort(unique(s$log_annual_precipitation))
  t_axis <- sort(unique(s$average_temperature))
  s <- s[order(s$log_annual_precipitation, s$average_temperature), ]
  lp <- log(annual_precipitation)
  i <- pmin(pmax(findInterval(lp, lp_axis), 1), length(lp_axis) - 1)
  j <- pmin(pmax(findInterval(average_temperature, t_axis), 1), length(t_axis) - 1)
  fx <- pmin(pmax((lp - lp_axis[i]) / (lp_axis[i + 1] - lp_axis[i]), 0), 1)
  fy <- pmin(pmax((average_temperature - t_axis[j]) / (t_axis[j + 1] - t_axis[j]), 0), 1)
  at <- function(ii, jj) (ii - 1) * length(t_axis) + jj
  corners <- list(at(i, j), at(i + 1, j), at(i, j + 1), at(i + 1, j + 1))
  w <- list((1 - fx) * (1 - fy), fx * (1 - fy), (1 - fx) * fy, fx * fy)
  interp <- function(v) Reduce(`+`, Map(function(k, wk) wk * v[k], corners, w))
  any_at <- function(v) Reduce(`|`, lapply(corners, function(k) v[k]))
  logged <- c("community_mean", "range_lower", "range_upper",
              "community_mean_lower", "community_mean_upper")
  out <- tibble::tibble(annual_precipitation, average_temperature)
  for (col in logged) out[[col]] <- exp(interp(log(s[[col]])))
  for (col in c("community_sd", "community_sd_lower", "community_sd_upper"))
    out[[col]] <- interp(s[[col]])
  out$extrapolated <- any_at(s$extrapolated)
  out$outside_australia <- any_at(!s$inside_australia) |
    annual_precipitation < exp(min(lp_axis)) | annual_precipitation > exp(max(lp_axis)) |
    average_temperature < min(t_axis) | average_temperature > max(t_axis)
  out
}

#' Write the three forecasting tables (parquet, csv.gz, dictionary) to `dir`.
#' @return the file paths.
write_forecast_tables <- function(responses, variances, surface, groups, centres, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(species_climate_responses = responses, site_variances = variances,
                 community_climate_surface = surface, growth_form_responses = groups,
                 species_climate_centres = centres)
  dict <- forecast_dictionary()
  unlist(purrr::imap(tables, function(t, stem) {
    p <- file.path(dir, paste0(stem, c(".parquet", ".csv.gz", "_dictionary.csv")))
    arrow::write_parquet(t, p[1])
    readr::write_csv(t, p[2])
    readr::write_csv(dplyr::filter(dict, column %in% names(t)), p[3])
    p
  }), use.names = FALSE)
}

#' Column dictionary of the forecasting tables.
forecast_dictionary <- function() {
  tibble::tribble(
    ~column, ~definition,
    "taxon_name", "Species name, as in the modelled means.",
    "trait_name", "AusTraits trait name.",
    "unit", "Unit of the trait (estimates are back-transformed to it).",
    "woodiness_group", "woody or non_woody: the growth form whose coefficients apply.",
    "centre_precipitation", "Annual precipitation at the species' climate centre (mm).",
    "centre_temperature", "Mean annual temperature at the species' climate centre (deg C).",
    "log_value_centre", "c_s: the species' modelled mean at its climate centre and an average location, natural log (equals log_value of the modelled means).",
    "slope_precipitation", "beta_P,s: the species' total within-species slope on log annual precipitation (growth form's + species' own), log trait per log mm.",
    "slope_temperature", "beta_T,s: likewise on mean annual temperature, log trait per deg C.",
    "var_centre", "Variance of log_value_centre.",
    "var_slope_precipitation", "Variance of slope_precipitation.",
    "var_slope_temperature", "Variance of slope_temperature.",
    "cov_centre_slope_precipitation", "Covariance of log_value_centre and slope_precipitation.",
    "cov_centre_slope_temperature", "Covariance of log_value_centre and slope_temperature.",
    "cov_slope_precipitation_temperature", "Covariance of the two slopes.",
    "n_records", "Records of the species in the fit.",
    "record_precipitation_min", "Lowest annual precipitation among the species' records (mm); below it a forecast extrapolates.",
    "record_precipitation_max", "Highest annual precipitation among the species' records (mm).",
    "record_temperature_min", "Lowest mean annual temperature among the species' records (deg C).",
    "record_temperature_max", "Highest mean annual temperature among the species' records (deg C).",
    "sd_location", "SD of the location effect (natural log): how sites of the same climate differ; add its square for a new site.",
    "sd_species_location", "SD of the species x location effect (natural log); add its square for a new site.",
    "sd_residual", "Residual SD (natural log) among records of a species at a site; add its square for a single new record.",
    "annual_precipitation", "Annual precipitation of the grid point (mm).",
    "log_annual_precipitation", "Its natural log (the grid is regular in it).",
    "average_temperature", "Mean annual temperature of the grid point (deg C).",
    "extrapolated", "The climate lies outside that of the trait's records for the growth form.",
    "centres_extrapolated", "The climate lies outside the cells the co-occurring species' centre surfaces were fitted to.",
    "inside_australia", "Some 50 km cell of Australia has a climate within one grid step on both axes.",
    "k", "Sorting factor applied to the species spread (community_sorting.R).",
    "community_mean", "Community mean at that climate: geometric mean of the species' means, at an average location (trait units).",
    "community_sd", "Community SD: SD of the species' means within a location (natural log).",
    "range_lower", "Lower end of the community range, the central 90% of species' means (trait units).",
    "range_upper", "Upper end of the community range (trait units).",
    "community_mean_lower", "Lower end of the 94% interval of community_mean (parameter uncertainty).",
    "community_mean_upper", "Upper end of the 94% interval of community_mean.",
    "community_sd_lower", "Lower end of the 94% interval of community_sd.",
    "community_sd_upper", "Upper end of the 94% interval of community_sd.",
    "n_species", "Species of the growth form in the fit.",
    "coef_intercept", "Growth form's intercept (natural log of the trait).",
    "coef_precipitation_between", "Between-species slope on the log annual precipitation of the species' centre.",
    "coef_temperature_between", "Between-species slope on the mean annual temperature of the species' centre.",
    "coef_precipitation_within", "Within-species slope on log annual precipitation (climate minus the species' centre).",
    "coef_temperature_within", "Within-species slope on mean annual temperature (climate minus the species' centre).",
    "cov_*__*", "Covariance of two of the five coefficients (upper triangle, diagonal included).",
    "spread_intercept", "log SD of species about the growth form's line, at the mean standardised centre.",
    "spread_precipitation", "Change in that log SD per SD of the species' log-precipitation centre (spread_scale_precipitation_*).",
    "spread_temperature", "Change in that log SD per SD of the species' temperature centre (spread_scale_temperature_*).",
    "spread_scale_precipitation_mean", "Mean of the fitted species' log-precipitation centres, for standardising.",
    "spread_scale_precipitation_sd", "SD of the fitted species' log-precipitation centres.",
    "spread_scale_temperature_mean", "Mean of the fitted species' temperature centres.",
    "spread_scale_temperature_sd", "SD of the fitted species' temperature centres.",
    "spread_log_sd_min", "Lowest log SD among the fitted species: the spread is clamped to this range.",
    "spread_log_sd_max", "Highest log SD among the fitted species.",
    "sd_slope_precipitation", "SD of species' own within-species slopes on log precipitation, about the growth form's.",
    "sd_slope_temperature", "SD of species' own within-species slopes on temperature.",
    "centre_precipitation_min", "Lowest annual precipitation among the fitted species' centres (mm).",
    "centre_precipitation_max", "Highest annual precipitation among the fitted species' centres (mm).",
    "centre_temperature_min", "Lowest mean annual temperature among the fitted species' centres (deg C).",
    "centre_temperature_max", "Highest mean annual temperature among the fitted species' centres (deg C)."
  )
}
