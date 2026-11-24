# Observed means, geometric method (`estimate_type = observed_mean`,
# `method = geometric`; #15), as AusTraits publishes them.
#
# This repository no longer computes them. austraits.build computes them from
# each build, in `R/observed_means.R`, and writes them at export with its
# `post_export()` (traitecoevo/austraits.build#899) as
# `austraits-<version>-observed-means.parquet`. The method lives there; change
# it there. The pipeline reads that file, for the same build as the records
# (`austraits_dev_rds` and `austraits_dev_observed_means` in `_targets.R`).

#' Columns of the published observed means this pipeline reads.
observed_means_columns <- c(
  "taxon_name", "trait_name", "unit", "estimate_type", "method", "type",
  "value", "log_value", "arithmetic_mean", "n_records", "value_mode", "value_count"
)

#' Read the observed means an austraits.build export published.
#'
#' Stops if the file is not the geometric observed means, lacks a column the
#' pipeline reads, or lacks one of `traits`.
#'
#' @param path the export's `austraits-<version>-observed-means.parquet`.
#' @param traits traits the pipeline needs means for.
#' @return the table, as published.
read_observed_means <- function(path, traits) {
  means <- arrow::read_parquet(path)

  missing <- setdiff(observed_means_columns, names(means))
  if (length(missing) > 0) {
    stop("`", path, "` lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!identical(unique(means$estimate_type), "observed_mean") ||
      !all(means$method %in% c("geometric", "arithmetic", "mode", "majority"))) {
    stop("`", path, "` is not the observed means", call. = FALSE)
  }
  # The modelled traits are on the log scale: their observed means must be
  # geometric (a signed trait's are arithmetic, terminology.md)
  not_geometric <- unique(means$trait_name[means$trait_name %in% traits &
                                             means$type == "numerical" &
                                             means$method != "geometric"])
  if (length(not_geometric) > 0) {
    stop("`", path, "` has arithmetic, not geometric, means for: ",
         paste(not_geometric, collapse = ", "), call. = FALSE)
  }
  absent <- setdiff(traits, means$trait_name)
  if (length(absent) > 0) {
    stop("`", path, "` has no observed means for: ", paste(absent, collapse = ", "),
         call. = FALSE)
  }

  means
}
