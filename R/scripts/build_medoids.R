# Build the climate medoid of every taxon in a cleaned GBIF extract (#33).
#
# Run from the project root, deliberately (not part of tar_make()):
#   Rscript R/scripts/build_medoids.R <download key> [workers]
# e.g.
#   Rscript R/scripts/build_medoids.R 0027765-260928105237408 10
#
# Reads data-processed/data_gbif_cleaned_<key>.parquet and writes
# data-processed/climate_centres_<key>.parquet plus a provenance JSON. Never
# overwrites climate_centres.parquet; the pipeline points at a medoid file
# through `medoids_file` in _targets.R.

source("R/init.R")

args <- commandArgs(trailingOnly = TRUE)
key <- args[[1]]
workers <- if (length(args) >= 2) as.integer(args[[2]]) else 10L
gbif_path <- file.path("data-processed", paste0("data_gbif_cleaned_", key, ".parquet"))
out <- file.path("data-processed", paste0("climate_centres_", key, ".parquet"))
stopifnot(file.exists(gbif_path))

t0 <- Sys.time()
medoids <- build_climate_medoids(gbif_path, workers = workers)
arrow::write_parquet(medoids, out)
jsonlite::write_json(list(
  gbif_extract = basename(gbif_path),
  gbif_download = key,
  method = "locate_medoid_2d() on distinct points in standardised log precipitation and temperature; taxa with more than 20,000 points subsampled (seed 123); 1% most outlying points dropped (at most 30)",
  taxa = nrow(medoids),
  workers = workers,
  built = format(Sys.time(), "%Y-%m-%d"),
  minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
), sub("\\.parquet$", ".json", out), auto_unbox = TRUE, pretty = TRUE)
message("Done: ", nrow(medoids), " taxa in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min -> ", out)
