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
  expect_equal(sqrt(.vc_species(vc, "woody")), 0.4, tolerance = 0.25)
  expect_equal(sqrt(.vc(vc, "location", "woodiness_groupwoody")), 0.4, tolerance = 0.25)
  expect_equal(sqrt(.vc(vc, "residual", "sigma")), 0.2, tolerance = 0.1)

  b_new <- summarise_community_benchmark(benchmark_community_variance(new, noise_vc = vc))
  b_old <- summarise_community_benchmark(benchmark_community_variance(old, noise_vc = vc))
  expect_true(all(abs(b_new$ratio - 1) < 0.2))
  expect_true(all(b_old$ratio > 1.25))
  expect_true(all(b_old$model_community_sd > b_new$model_community_sd))
})

test_that("community_range back-transforms a central interval", {
  r <- community_range(log(10), 0.5^2, level = 0.90)
  expect_equal(r$community_mean, 10)
  expect_equal(r$range_lower * r$range_upper, 100)
  expect_equal(log(r$range_upper / 10), stats::qnorm(0.95) * 0.5)
})
