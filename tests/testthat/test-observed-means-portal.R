# The portal method's behaviour, quirks included, on a small constructed
# dataset; then the reproduction check against the deployed portal file.

toy_records <- function() {
  tibble::tribble(
    ~taxon_name, ~trait_name, ~unit, ~value_type, ~basis_of_record, ~value,
    ~dataset_id, ~location_id, ~observation_id, ~original_name,
    # two locations of field records: location means 2 and 6
    "Acacia a", "plant_height", "m", "raw", "field", "1", "D1", "L1", "01", "A a",
    "Acacia a", "plant_height", "m", "raw", "field", "3", "D1", "L1", "02", "A a",
    "Acacia a", "plant_height", "m", "mean", "field", "6", "D1", "L2", "03", "A a",
    # a flora observation with both bounds: midpoint 10
    "Acacia a", "plant_height", "m", "minimum", "preserved_specimen", "8", "F1", NA, "01", "A a",
    "Acacia a", "plant_height", "m", "maximum", "preserved_specimen", "12", "F1", NA, "01", "A a",
    # a lone flora bound, used as its own midpoint
    "Acacia b", "plant_height", "m", "maximum", "literature", "4", "F1", NA, "02", "A b",
    # a raw record from a glasshouse: not filtered out by basis_of_record
    "Acacia b", "plant_height", "m", "raw", "experiment", "2", "D2", "L9", "01", "A b",
    # a bound from the field: not a flora record, so not used
    "Acacia b", "plant_height", "m", "maximum", "field", "100", "D3", "L8", "01", "A b",
    # categorical, two space-separated values in one record
    "Acacia a", "woodiness", NA, "mode", "literature", "woody", "D1", NA, "04", "A a",
    "Acacia a", "woodiness", NA, "mode", "literature", "woody semi_woody", "F1", NA, "03", "A a"
  ) |>
    dplyr::mutate(
      `latitude (deg)` = NA_character_, `longitude (deg)` = NA_character_,
      location_name = NA_character_, source_primary_key = dataset_id,
      source_primary_citation = paste("Citation", dataset_id),
      taxon_rank = "species", taxon_name_alternatives = NA_character_,
      family = "Fabaceae", scientific_name_id = paste0("id:", taxon_name)
    )
}

test_that("species value is the arithmetic mean of location means and flora midpoints", {
  out <- estimate_observed_means_portal(toy_records())
  a <- dplyr::filter(out, taxon_name == "Acacia a", trait_name == "plant_height")
  expect_equal(nrow(a), 1)
  expect_equal(a$value_mean, mean(c(2, 6, 10)))
  expect_equal(a$location_replicates, 2)
  expect_equal(a$flora_replicates, 1)
  expect_equal(a$all_replicates, 3 + 1)
  expect_equal(a$value_min, 1)
  expect_equal(a$value_max, 12)
})

test_that("value_geom_mean equals value_mean (kept quirk)", {
  out <- estimate_observed_means_portal(toy_records())
  num <- dplyr::filter(out, type == "numerical")
  expect_equal(num$value_geom_mean, num$value_mean)
})

test_that("lone flora bounds are used and raw records are not filtered by basis_of_record", {
  out <- estimate_observed_means_portal(toy_records())
  b <- dplyr::filter(out, taxon_name == "Acacia b", trait_name == "plant_height")
  expect_equal(b$value_mean, mean(c(4, 2)))
  expect_equal(b$dataset_id, "D2; F1")
})

test_that("categorical values are split on spaces and counted", {
  out <- estimate_observed_means_portal(toy_records())
  w <- dplyr::filter(out, trait_name == "woodiness")
  expect_true(all(w$value_count == "semi_woody (1); woody (2)"))
  expect_equal(unique(w$type), "categorical")
})

test_that("the pipeline reproduces the deployed portal species averages for 7.0.0", {
  store <- file.path(project_root, "_targets")
  skip_if_not(dir.exists(store), "targets store not built")
  check <- tryCatch(targets::tar_read(portal_reproduction, store = store),
                    error = function(e) NULL)
  skip_if(is.null(check), "portal_reproduction target not built")
  expect_true(check$pass)
  expect_equal(sum(check$columns$n_mismatch), 0)
})
