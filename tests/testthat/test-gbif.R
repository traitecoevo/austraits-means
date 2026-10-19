# The record-level GBIF filters (filter_gbif_records()) on constructed records:
# each row breaks exactly one rule, except the first, which passes them all.

toy_gbif <- function() {
  base <- tibble::tibble(
    license = "CC_BY_4_0", countryCode = "AU", occurrenceStatus = "PRESENT",
    decimalLongitude = 150.1, decimalLatitude = -33.5,
    basisOfRecord = "PRESERVED_SPECIMEN", coordinatePrecision = NA_real_,
    year = 2001L, coordinateUncertaintyInMeters = 50, issue = "",
    taxonRank = "SPECIES", species = "Acacia linifolia", institutionCode = "NSW"
  )
  breaks <- list(
    keep = list(),
    licence_nc = list(license = "CC_BY_NC_4_0"),
    country = list(countryCode = "NZ"),
    absent = list(occurrenceStatus = "ABSENT"),
    no_coordinates = list(decimalLatitude = NA_real_),
    living = list(basisOfRecord = "LIVING_SPECIMEN"),
    imprecise = list(coordinatePrecision = 0.1),
    old = list(year = 1850L),
    uncertain = list(coordinateUncertaintyInMeters = 20000),
    zero = list(decimalLongitude = 0),
    mismatch = list(issue = "COUNTRY_COORDINATE_MISMATCH;TAXON_MATCH_FUZZY"),
    genus_rank = list(taxonRank = "GENUS"),
    no_species = list(species = ""),
    nhmuk = list(institutionCode = "NHMUK")
  )
  purrr::imap_dfr(breaks, function(change, id) {
    row <- base
    for (nm in names(change)) row[[nm]] <- change[[nm]]
    dplyr::mutate(row, id = id)
  })
}

test_that("each GBIF record filter removes exactly the records that break it", {
  kept <- filter_gbif_records(toy_gbif())
  expect_equal(kept$id, "keep")
})

test_that("CC0 and CC BY records are kept, CC BY-NC records dropped", {
  x <- dplyr::filter(toy_gbif(), id == "keep") |>
    tidyr::uncount(3) |>
    dplyr::mutate(license = c("CC0_1_0", "CC_BY_4_0", "CC_BY_NC_4_0"))
  expect_equal(filter_gbif_records(x)$license, c("CC0_1_0", "CC_BY_4_0"))
})

test_that("missing coordinate precision and uncertainty are allowed", {
  x <- dplyr::filter(toy_gbif(), id == "keep") |>
    dplyr::mutate(coordinatePrecision = NA_real_, coordinateUncertaintyInMeters = NA_real_)
  expect_equal(nrow(filter_gbif_records(x)), 1)
})

test_that("a download's citation uses its DOI and date", {
  expect_equal(gbif_citation("0027765-260928105237408"),
               "GBIF.org (15 October 2026) GBIF Occurrence Download https://doi.org/10.15468/dl.766erw")
})

test_that("names align from the provider name first, GBIF's name as fallback, rolled up to species", {
  skip_if_not(dir.exists(tools::R_user_dir("APCalign")), "APCalign resources not cached")
  cleaned <- tibble::tribble(
    ~verbatimScientificName,                  ~species,
    "Corymbia maculata",                      "Blakella maculata",   # GBIF genus APC does not use
    "",                                       "Acacia aneura",       # no provider name
    "Eucalyptus leucoxylon subsp. megalocarpa", "Eucalyptus leucoxylon"
  )
  a <- align_gbif_names(cleaned, apc_version = "2026-07-30")
  a <- a[match(cleaned$species, a$gbif_name), ]
  expect_equal(a$name_source, c("provider", "gbif", "provider"))
  expect_equal(a$taxon_name, c("Corymbia maculata", "Acacia aneura", "Eucalyptus leucoxylon"))
  expect_equal(a$taxon_name_infraspecific[3], "Eucalyptus leucoxylon subsp. megalocarpa")
})

test_that("fuzzy matches that change the genus are rejected, spelling variants kept", {
  cleaned <- tibble::tribble(
    ~verbatimScientificName,   ~species,
    "Acanthus spinosus",       "Acanthus spinosus",
    "Leschenaultia linioides", "Lechenaultia linarioides"
  )
  lookup <- tibble::tribble(
    ~original_name,             ~aligned_name,              ~accepted_name,             ~taxon_rank, ~genus,         ~aligned_reason,
    "Acanthus spinosus",        "Amaranthus spinosus",      "Amaranthus spinosus",      "species",   "Amaranthus",   "Fuzzy match of taxon name",
    "Leschenaultia linioides",  "Lechenaultia linarioides", "Lechenaultia linarioides", "species",   "Lechenaultia", "Fuzzy match of taxon name",
    "Lechenaultia linarioides", "Lechenaultia linarioides", "Lechenaultia linarioides", "species",   "Lechenaultia", "Exact match"
  )
  a <- align_gbif_names(cleaned, apc_version = "2026-07-30", lookup = lookup)
  a <- a[match(cleaned$species, a$gbif_name), ]
  expect_equal(a$name_source, c("none", "provider"))
  expect_equal(a$taxon_name[2], "Lechenaultia linarioides")
})
