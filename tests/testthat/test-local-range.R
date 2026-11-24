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

test_that("community_range back-transforms a central interval", {
  r <- community_range(log(10), 0.5^2, level = 0.90)
  expect_equal(r$community_mean, 10)
  expect_equal(r$range_lower * r$range_upper, 100)
  expect_equal(log(r$range_upper / 10), stats::qnorm(0.95) * 0.5)
})
