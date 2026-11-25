# Location model prototype (#41, R/model/local_range.R): the community
# variance is the spread among species of one growth form at one location,
# and must not be built from a model that pools all records of a climate.

test_that("add_location_key re-keys data that already has a location", {
  d <- tibble::tibble(decimal_longitude = c(150.0, 150.3, 151.5),
                      decimal_latitude = c(-33.0, -33.1, -34.0),
                      dataset_id = "a", location_id = c("x", "y", "z"))
  k10 <- add_location_key(d, "cell", cell_km = 10)
  k100 <- add_location_key(k10, "cell", cell_km = 100)
  expect_equal(dplyr::n_distinct(k10$location), 3)
  expect_lt(dplyr::n_distinct(k100$location), 3)
  expect_equal(add_location_key(k10, "site")$location, c("a|x", "a|y", "a|z"))
})

test_that("the location model recovers known variances; the pooled model overstates the community SD", {
  sim <- simulate_location_data(seed = 2)
  old <- suppressMessages(fit_woodiness_model(sim$data, "sim", sim$groups, n_min = 0))
  new <- suppressMessages(fit_woodiness_model(sim$data, "sim", sim$groups, n_min = 0,
                                              location = "location"))
  vc <- tidy_varcomp(new)
  # every component is found under the name the benchmark looks up
  expect_setequal(unique(vc$grouping),
                  c("taxon_name", "location", "taxon_name:location", "residual"))
  expect_equal(sqrt(.vc(vc, "taxon_name:location", "(Intercept)")), 0.15, tolerance = 0.3)
  expect_equal(sqrt(.vc_species(vc, "woody")), 0.4, tolerance = 0.25)
  expect_equal(sqrt(.vc(vc, "location", "woodiness_groupwoody")), 0.4, tolerance = 0.25)
  expect_equal(sqrt(.vc(vc, "residual", "sigma")), 0.2, tolerance = 0.1)

  b_new <- summarise_community_benchmark(benchmark_community_variance(new, noise_vc = vc))
  b_old <- summarise_community_benchmark(benchmark_community_variance(old, noise_vc = vc))
  expect_true(all(abs(b_new$ratio - 1) < 0.2))
  expect_true(all(b_old$ratio > 1.25))
  expect_true(all(b_old$model_community_sd > b_new$model_community_sd))
})

test_that("community_diagnostics gives one row per location and form, with the location effect", {
  sim <- simulate_location_data(seed = 3)
  new <- suppressMessages(fit_woodiness_model(sim$data, "sim", sim$groups, n_min = 0,
                                              location = "location"))
  dg <- community_diagnostics(new)
  expect_equal(nrow(dg$cells), nrow(dplyr::distinct(dg$species, location, woodiness_group)))
  expect_true(all(dg$cells$n_species >= 5))
  expect_true(any(dg$cells$location_effect != 0))
  expect_equal(dg$cells$community_mean_log, dg$cells$fixed_mean + dg$cells$location_effect)
  expect_true(all(dg$cells$range_lower_log < dg$cells$community_mean_log))
  # the cell means track the observed species means better with the location effect
  obs <- dplyr::summarise(dplyr::group_by(dg$species, location, woodiness_group),
                          y = mean(y), .groups = "drop")
  m <- dplyr::inner_join(dg$cells, obs, by = c("location", "woodiness_group"))
  expect_lt(mean((m$y - m$community_mean_log)^2), mean((m$y - m$fixed_mean)^2))
  p <- plot_community_range_cells(list(trait = "sim", unit = "u", cell_km = 50,
                                       diagnostics = list(location = dg, pooled = dg)))
  expect_s3_class(p, "ggplot")
})

test_that("species variance by climate class: balanced classes, one variance each", {
  sim <- simulate_location_data(seed = 4)
  d <- dplyr::left_join(sim$data, sim$groups, by = "taxon_name")
  cl <- add_spread_class(d, "temperature", n_bins = 3)
  expect_equal(nrow(cl$bins), 6)
  expect_true(all(abs(cl$bins$n_species - 50) <= 1))
  expect_true(all(cl$bins$lower <= cl$bins$centre & cl$bins$centre <= cl$bins$upper))
  expect_false(anyNA(cl$data$spread_class))

  m <- suppressMessages(fit_woodiness_model(sim$data, "sim", sim$groups, n_min = 0,
                                            location = "location", spread_by = "temperature",
                                            n_bins = 3))
  s <- spread_by_climate(m)
  expect_equal(nrow(s), 6)
  expect_true(all(s$species_sd > 0))
  # the simulation's species SD is 0.4 in every class
  expect_equal(mean(s$species_sd), 0.4, tolerance = 0.3)
  sv <- .species_variances(m, tidy_varcomp(m))
  expect_equal(nrow(sv), dplyr::n_distinct(sim$data$taxon_name))
  b <- benchmark_community_variance(m)
  expect_true(all(is.finite(b$model_community_variance)))
})

test_that("RTMB spread model reproduces glmmTMB with constant spread", {
  skip_if_not_installed("RTMB")
  sim <- simulate_location_data(seed = 5)
  g <- suppressMessages(fit_woodiness_model(sim$data, "sim", sim$groups, n_min = 0,
                                            location = "location"))
  r <- suppressMessages(fit_spread_rtmb(sim$data, "sim", sim$groups, n_min = 0,
                                        spread_terms = character(0)))
  expect_true(r$converged)
  expect_equal(r$logLik, as.numeric(stats::logLik(g$fit)), tolerance = 1e-4)
  expect_equal(r$df, attr(stats::logLik(g$fit), "df"))
  expect_equal(unname(exp(r$spread$theta)),
               sqrt(c(.vc_species(tidy_varcomp(g), "woody"),
                      .vc_species(tidy_varcomp(g), "non_woody"))), tolerance = 1e-3)
  b_g <- benchmark_community_variance(g)
  b_r <- benchmark_community_variance(r)
  expect_equal(b_r$model_community_variance, b_g$model_community_variance, tolerance = 1e-3)
})

test_that("RTMB spread model recovers spread that changes with climate", {
  skip_if_not_installed("RTMB")
  sim <- simulate_location_data(seed = 6, n_species = 600, spread_slope = 0.4)
  r <- suppressMessages(fit_spread_rtmb(sim$data, "sim", sim$groups, n_min = 0,
                                        spread_terms = c("p", "t")))
  expect_true(r$converged)
  th <- r$spread$theta
  se <- sqrt(diag(r$spread$cov))
  # true: 0.4 per SD of temperature centre, 0 for precipitation
  expect_true(all(abs(th[c("woody:t", "non_woody:t")] - 0.4) < 3 * se[c("woody:t", "non_woody:t")]))
  expect_true(all(abs(th[c("woody:p", "non_woody:p")]) < 3 * se[c("woody:p", "non_woody:p")]))
  cv <- spread_curve(r, "temperature")
  expect_true(all(diff(cv$species_sd[cv$woodiness_group == "woody"]) > 0))
})

test_that("shrinking the spread terms keeps a real slope and drops the rest", {
  skip_if_not_installed("RTMB")
  # true: log species SD linear in temperature (0.4 per SD), no precipitation
  # or quadratic effect. The quadratic terms are not the default; they are
  # fitted here to check that shrinkage drops a degree the data don't support.
  all_terms <- c("p", "t", "p2", "t2")
  sim <- simulate_location_data(seed = 6, n_species = 600, spread_slope = 0.4)
  free <- suppressMessages(fit_spread_rtmb(sim$data, "sim", sim$groups, n_min = 0,
                                           spread_terms = all_terms))
  r <- suppressMessages(fit_spread_shrunk(sim$data, "sim", sim$groups, free, n_min = 0))
  expect_true(r$converged)
  sh <- r$spread$shrinkage
  expect_gt(sh$tau[["degree 1"]], 0.1)
  expect_lt(sh$tau[["degree 2"]], sh$tau[["degree 1"]] / 3)
  th <- r$spread$theta
  se <- sqrt(diag(r$spread$cov))
  expect_true(all(abs(th[c("woody:t", "non_woody:t")] - 0.4) < 3 * se[c("woody:t", "non_woody:t")]))
  # every null term is pulled toward 0, relative to the free fit
  null_terms <- names(th)[grepl(":(p|p2|t2)$", names(th))]
  expect_true(all(abs(th[null_terms]) < abs(free$spread$theta[null_terms])))
  # the refit's mode agrees with the Gaussian approximation that chose tau
  expect_equal(unname(th[names(sh$b_post)]), unname(sh$b_post), tolerance = 0.05)
  expect_gt(sh$edf, 2)
  expect_lt(sh$edf, 8)
  expect_lt(r$aic, free$aic)
  expect_gt(sh$log_bf, 10)

  # no climate effect at all: every degree is dropped
  sim0 <- simulate_location_data(seed = 7, n_species = 300, spread_slope = 0)
  free0 <- suppressMessages(fit_spread_rtmb(sim0$data, "sim", sim0$groups, n_min = 0,
                                            spread_terms = all_terms))
  sh0 <- estimate_spread_shrinkage(free0$spread$theta, free0$spread$cov)
  expect_true(sh0$at_lower[["degree 1"]])
  expect_true(all(abs(sh0$b_post) < 0.03))
  expect_lt(sh0$log_bf, 1)
})

test_that("species means from the RTMB fit match glmmTMB's, at an average location", {
  skip_if_not_installed("RTMB")
  sim <- simulate_location_data(seed = 5)
  g <- suppressMessages(fit_woodiness_model(sim$data, "sim", sim$groups, n_min = 0,
                                            location = "location"))
  r <- suppressMessages(fit_spread_rtmb(sim$data, "sim", sim$groups, n_min = 0,
                                        spread_terms = character(0), centres = TRUE))
  cn <- r$centres
  expect_equal(nrow(cn), dplyr::n_distinct(sim$data$taxon_name))
  nd <- dplyr::mutate(cn, log_annual_precipitation_c = 0, average_temperature_c = 0,
                      location = "__new__")
  p <- stats::predict(g$fit, newdata = nd, se.fit = TRUE, allow.new.levels = TRUE)
  expect_equal(cn$log_estimate, unname(p$fit), tolerance = 1e-4)
  # glmmTMB's SE is for a new location: it adds the location and species x
  # location variances, which an average location leaves out
  v_loc <- stats::setNames(r$vc$sd[r$vc$grouping == "location"]^2, r$groups)
  v_w <- r$vc$sd[r$vc$grouping == "taxon_name:location"]^2
  expect_equal(unname(sqrt(cn$log_se^2 + v_loc[as.character(cn$woodiness_group)] + v_w)),
               unname(p$se.fit), tolerance = 1e-4)
  expect_true(all(cn$log_conf.low < cn$log_estimate & cn$log_estimate < cn$log_conf.high))
  # the species' own effect, as glmmTMB estimates it
  re <- glmmTMB::ranef(g$fit)$cond$taxon_name
  own <- re[cbind(match(cn$taxon_name, rownames(re)),
                  match(paste0("woodiness_group", cn$woodiness_group), colnames(re)))]
  expect_equal(cn$g_0, own, tolerance = 1e-4)
  # with the climate terms shrunk, the centres come with the fit
  free <- suppressMessages(fit_spread_rtmb(sim$data, "sim", sim$groups, n_min = 0))
  s <- suppressMessages(fit_spread_shrunk(sim$data, "sim", sim$groups, free, n_min = 0,
                                          centres = TRUE))
  expect_equal(nrow(s$centres), nrow(cn))
  expect_true(all(s$centres$log_se > 0))
})

test_that("community statistics from the centres' distribution equal those from the species", {
  skip_if_not_installed("RTMB")
  gh <- .gauss_hermite_normal(7)
  expect_equal(c(sum(gh$w), sum(gh$w * gh$z^2), sum(gh$w * gh$z^4)), c(1, 1, 3))

  sim <- simulate_location_data(seed = 5)
  sim$data <- dplyr::mutate(
    sim$data,
    log_annual_precipitation = log_annual_precipitation_mn + log_annual_precipitation_c,
    average_temperature = average_temperature_mn + average_temperature_c)
  m <- suppressMessages(fit_spread_rtmb(sim$data, "sim", sim$groups, n_min = 0,
                                        spread_terms = character(0)))
  # species present in 30 cells, drawn from the fitted species' centres
  set.seed(1)
  sp <- dplyr::distinct(sim$data, taxon_name, log_annual_precipitation_mn, average_temperature_mn) |>
    dplyr::inner_join(sim$groups, by = "taxon_name")
  climate <- tibble::tibble(cell_x = 1:30, cell_y = 0L,
                            log_annual_precipitation = stats::runif(30, 6, 7.5),
                            average_temperature = stats::runif(30, 12, 22))
  climate$annual_precipitation <- exp(climate$log_annual_precipitation)
  species <- purrr::map(1:30, function(i) dplyr::slice_sample(sp, n = 40) |>
    dplyr::transmute(cell_x = i, cell_y = 0L, gbif_taxon_name = taxon_name, woodiness_group,
                     centre_log_annual_precipitation = log_annual_precipitation_mn,
                     centre_average_temperature = average_temperature_mn)) |>
    purrr::list_rbind()
  by_species <- community_statistics(list(m), climate, species)
  # the cells' own centre distributions, as population SDs
  centres <- cell_centre_summaries(species, climate) |>
    dplyr::mutate(sd_p = sd_p * sqrt((n_species - 1) / n_species),
                  sd_t = sd_t * sqrt((n_species - 1) / n_species),
                  centres_extrapolated = FALSE)
  by_climate <- community_statistics_climate(list(m), centres, n_draws = 0)
  both <- dplyr::inner_join(by_species, by_climate,
                            by = c("woodiness_group", "cell_x", "cell_y"), suffix = c("_s", "_c"))
  expect_equal(nrow(both), 60)
  expect_equal(both$community_mean_c, both$community_mean_s, tolerance = 1e-10)
  expect_equal(both$community_sd_c, both$community_sd_s, tolerance = 1e-10)

  # intervals from parameter draws hold the estimate
  ci <- community_statistics_climate(list(m), centres, n_draws = 100)
  expect_gt(mean(ci$community_mean_lower <= ci$community_mean &
                   ci$community_mean <= ci$community_mean_upper), 0.95)
  expect_true(all(ci$community_sd_lower < ci$community_sd_upper))
  # a sorting factor below 1 narrows the range, in both versions alike
  srt <- tibble::tibble(trait = "sim", woodiness_group = c("woody", "non_woody"), k = 0.5)
  half_c <- community_statistics_climate(list(m), centres, sorting = srt, n_draws = 0)
  half_s <- community_statistics(list(m), climate, species, sorting = srt)
  expect_true(all(half_c$community_sd < by_climate$community_sd))
  expect_equal(dplyr::arrange(half_c, woodiness_group, cell_x)$community_sd,
               dplyr::arrange(half_s, woodiness_group, cell_x)$community_sd, tolerance = 1e-10)
})

test_that("forecasts from the exported tables reproduce the model", {
  skip_if_not_installed("RTMB")
  sim <- simulate_location_data(seed = 5)
  sim$data <- dplyr::mutate(
    sim$data, unit = "u",
    log_annual_precipitation = log_annual_precipitation_mn + log_annual_precipitation_c,
    average_temperature = average_temperature_mn + average_temperature_c)
  g <- suppressMessages(fit_woodiness_model(sim$data, "sim", sim$groups, n_min = 0,
                                            location = "location"))
  m <- suppressMessages(fit_spread_rtmb(sim$data, "sim", sim$groups, n_min = 0,
                                        spread_terms = character(0), centres = TRUE))
  m$centres <- location_species_centres(m)
  responses <- species_response_table(list(m))
  variances <- site_variance_table(list(m))
  # a new site away from each species' centre, against glmmTMB's prediction
  # for a new location (which carries the location and species x location
  # variances, as type = "site" does)
  sp <- responses[1:20, ]
  f <- purrr::pmap(sp, function(taxon_name, centre_precipitation, centre_temperature, ...)
    forecast_species(responses, variances, taxon_name, "sim",
                     centre_precipitation * exp(0.3), centre_temperature - 2, type = "site")) |>
    dplyr::bind_rows()
  nd <- tibble::tibble(taxon_name = sp$taxon_name, woodiness_group = sp$woodiness_group,
                       log_annual_precipitation_mn = log(sp$centre_precipitation),
                       average_temperature_mn = sp$centre_temperature,
                       log_annual_precipitation_c = 0.3, average_temperature_c = -2,
                       location = "__new__")
  p <- stats::predict(g$fit, newdata = nd, se.fit = TRUE, allow.new.levels = TRUE)
  expect_equal(f$log_estimate, unname(p$fit), tolerance = 1e-4)
  expect_equal(unname(f$log_se), unname(p$se.fit), tolerance = 1e-4)
  expect_true(all(forecast_species(responses, variances, sp$taxon_name[1], "sim",
                                   1e5, 50)$extrapolated))
})

test_that("forecasts for species without trait data are the growth form's line plus the species spread", {
  skip_if_not_installed("RTMB")
  sim <- simulate_location_data(seed = 5)
  sim$data <- dplyr::mutate(
    sim$data, unit = "u",
    log_annual_precipitation = log_annual_precipitation_mn + log_annual_precipitation_c,
    average_temperature = average_temperature_mn + average_temperature_c)
  g <- suppressMessages(fit_woodiness_model(sim$data, "sim", sim$groups, n_min = 0,
                                            location = "location"))
  m <- suppressMessages(fit_spread_rtmb(sim$data, "sim", sim$groups, n_min = 0,
                                        spread_terms = character(0), centres = TRUE))
  groups <- growth_form_response_table(list(m))
  variances <- site_variance_table(list(m))
  expect_equal(nrow(groups), 2)
  # a species centred at (log 700 mm, 16 deg C), forecast at 900 mm and 18 deg C
  f <- forecast_unmeasured(groups, variances, "sim", "woody", 700, 16, 900, 18)
  nd <- tibble::tibble(woodiness_group = factor("woody", m$groups),
                       log_annual_precipitation_mn = log(700), average_temperature_mn = 16,
                       log_annual_precipitation_c = log(900 / 700), average_temperature_c = 2,
                       taxon_name = NA, location = NA)
  p <- stats::predict(g$fit, newdata = nd, re.form = NA, se.fit = TRUE)
  expect_equal(f$log_estimate, unname(p$fit), tolerance = 1e-4)
  # the species terms: its effect (constant spread here) and its slopes
  vc <- m$vc
  sd_u2 <- exp(2 * m$spread$theta[["woody:(Intercept)"]])
  slopes <- vc$sd[vc$term == "log_annual_precipitation_c"]^2 * log(900 / 700)^2 +
    vc$sd[vc$term == "average_temperature_c"]^2 * 4
  expect_equal(f$log_se^2 - sd_u2 - slopes, unname(p$se.fit)^2, tolerance = 1e-3)
  # at a fitted species' centre, the line is its modelled mean less its own effect
  ce <- m$centres[1, ]
  f0 <- forecast_unmeasured(groups, variances, "sim", as.character(ce$woodiness_group),
                            exp(ce$log_annual_precipitation_mn), ce$average_temperature_mn,
                            exp(ce$log_annual_precipitation_mn), ce$average_temperature_mn)
  expect_equal(f0$log_estimate, ce$log_estimate - ce$g_0, tolerance = 1e-10)
  expect_true(forecast_unmeasured(groups, variances, "sim", "woody", 5, 50, 900, 18)$centre_extrapolated)
})

test_that("forecast_community interpolates the surface exactly at its nodes", {
  surface <- tidyr::expand_grid(log_annual_precipitation = log(c(200, 600, 1800)),
                                average_temperature = c(10, 20, 30)) |>
    dplyr::mutate(trait_name = "sim", woodiness_group = "woody",
                  annual_precipitation = exp(log_annual_precipitation),
                  community_mean = exp(log_annual_precipitation / 3 + average_temperature / 10),
                  range_lower = community_mean / 4, range_upper = community_mean * 4,
                  community_mean_lower = community_mean / 1.1,
                  community_mean_upper = community_mean * 1.1,
                  community_sd = 0.5 + average_temperature / 100,
                  community_sd_lower = community_sd - 0.05, community_sd_upper = community_sd + 0.05,
                  extrapolated = average_temperature > 25, inside_australia = TRUE)
  f <- forecast_community(surface, "sim", "woody", surface$annual_precipitation,
                          surface$average_temperature)
  expect_equal(f$community_mean, surface$community_mean, tolerance = 1e-10)
  expect_equal(f$community_sd, surface$community_sd, tolerance = 1e-10)
  # between nodes: linear in log precipitation and temperature
  mid <- forecast_community(surface, "sim", "woody", sqrt(200 * 600), 15)
  expect_equal(log(mid$community_mean), log(sqrt(200 * 600)) / 3 + 1.5, tolerance = 1e-10)
  expect_false(mid$extrapolated)
  expect_true(forecast_community(surface, "sim", "woody", 100, 15)$outside_australia)
})

test_that("community_range back-transforms a central interval", {
  r <- community_range(log(10), 0.5^2, level = 0.90)
  expect_equal(r$community_mean, 10)
  expect_equal(r$range_lower * r$range_upper, 100)
  expect_equal(log(r$range_upper / 10), stats::qnorm(0.95) * 0.5)
})
