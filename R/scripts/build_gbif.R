# Build the cleaned, climate-tagged GBIF extract from a download (#32).
#
# Run from the project root, deliberately (it reads a ~17 GB table and takes
# tens of minutes; it is not part of tar_make()):
#   Rscript R/scripts/build_gbif.R <download key> <folder holding the download>
# e.g.
#   Rscript R/scripts/build_gbif.R 0027765-260928105237408 data/GBIF/2026-10-15
#
# Writes, in data-processed/:
#   gbif_records_cleaned_<key>.parquet  the cleaned records before name
#                                     alignment (reused if present, so a new
#                                     APC release or APCalign version only
#                                     reruns the alignment and climate steps)
#   data_gbif_cleaned_<key>.parquet   the cleaned extract with climate
#   gbif_names_aligned_<key>.parquet  the APC alignment of its species names
#   data_gbif_cleaned_<key>.json      provenance: download, licences, APC
#                                     release, record count after each step
# It never overwrites data_gbif_cleaned.parquet, the extract the pipeline
# currently reads; switching to a new extract is a separate step (#33).

source("R/init.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
key <- args[[1]]
dir <- args[[2]]
stopifnot(key %in% gbif_downloads$key)

# The APC release pinned for the crosswalk in _targets.R.
apc_version <- local({
  line <- grep('^apc_version <- ', readLines("_targets.R"), value = TRUE)
  sub('^apc_version <- "([^"]+)".*', "\\1", line)
})

t0 <- Sys.time()
stage <- file.path("data-processed", paste0("gbif_records_cleaned_", key, ".parquet"))
if (file.exists(stage)) {
  message("Reusing ", stage)
  cleaned <- arrow::read_parquet(stage)
  steps <- unlist(jsonlite::read_json(sub("\\.parquet$", "_steps.json", stage)))
} else {
  message("Cleaning ", key)
  cleaned <- clean_gbif(download_gbif(key, dir))
  steps <- attr(cleaned, "steps")
  arrow::write_parquet(cleaned, stage)
  jsonlite::write_json(as.list(steps), sub("\\.parquet$", "_steps.json", stage),
                       auto_unbox = TRUE)
}
message("Aligning ", dplyr::n_distinct(cleaned$species), " names to APC ", apc_version)
aligned <- align_gbif_names(cleaned$species, apc_version)
message("Attaching climate")
gbif <- build_gbif_climate(cleaned, aligned)
steps["with_climate"] <- nrow(gbif)

out <- file.path("data-processed", paste0(c("data_gbif_cleaned_", "gbif_names_aligned_"),
                                          key, ".parquet"))
arrow::write_parquet(gbif, out[1])
arrow::write_parquet(aligned, out[2])
jsonlite::write_json(list(
  gbif_download = key,
  doi = gbif_downloads$doi[gbif_downloads$key == key],
  citation = gbif_citation(key),
  licences_kept = gbif_open_licences,
  apc_version = apc_version,
  apcalign_package = as.character(utils::packageVersion("APCalign")),
  records = as.list(steps),
  species_names = dplyr::n_distinct(gbif$original_name),
  taxa_aligned = dplyr::n_distinct(stats::na.omit(gbif$taxon_name)),
  built = format(Sys.time(), "%Y-%m-%d"),
  minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
), sub("\\.parquet$", ".json", out[1]), auto_unbox = TRUE, pretty = TRUE)
message("Done: ", paste(out, collapse = ", "))
print(steps)
