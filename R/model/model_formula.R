# Build the within-group-centred mixed-model formula.
#
# For climate variables v1, v2, ... (each present as `<v>_mn` = species medoid and
# `<v>_c` = obs - medoid), the model is the Bafumi-Gelman / Mundlak form with
# uncorrelated random intercept and per-variable random slopes (no climate
# interaction — matches the convergence-driven choices in the height analysis):
#
#   response ~ v1_mn + v2_mn + v1_c + v2_c
#              + (1 | group)
#              + (0 + v1_c | group) + (0 + v2_c | group)
#              [+ (1 | methodology)]

#' @param response response column (e.g. "log_value").
#' @param climate_vars character vector of climate variable stems.
#' @param group grouping factor (species), default "taxon_name".
#' @param methodology optional grouping factor added as a random intercept
#'   (the round-3 methodology term); NULL to omit.
#' @return a formula.
build_model_formula <- function(response = "log_value",
                                 climate_vars = c("log_annual_precipitation", "average_temperature"),
                                 group = "taxon_name",
                                 methodology = NULL) {
  mn <- paste0(climate_vars, "_mn")
  ce <- paste0(climate_vars, "_c")

  fixed <- c(mn, ce)
  re_intercept <- sprintf("(1 | %s)", group)
  re_slopes <- sprintf("(0 + %s | %s)", ce, group)

  rhs <- c(fixed, re_intercept, re_slopes)
  if (!is.null(methodology)) {
    rhs <- c(rhs, sprintf("(1 | %s)", methodology))
  }

  stats::as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}
