# The trait-climate model with woody and non-woody coefficients, in glmmTMB
# (fit_woodiness_model()): with `location`, the location model as glmmTMB
# can fit it (constant spread; the check on R/model/spread_rtmb.R); without,
# the "locations pooled" model, kept as a comparison
# (R/process/compare_pooled_location.R). .woodiness_data() and
# .group_coef_names() are shared with the RTMB location model.

woodiness_model_terms <- c("log_annual_precipitation_mn", "average_temperature_mn",
                           "log_annual_precipitation_c", "average_temperature_c")

#' Trait-climate model with woody and non-woody coefficients (glmmTMB).
#'
#' The model of `fit_trait_climate_model()` (species with more than `n_min`
#' records), with a separate intercept and climate slopes for woody and
#' non-woody species and the species random effects shared:
#' `log_value ~ 0 + group + group:(P_mn + T_mn + P_c + T_c) + (1 | taxon) +
#' (0 + P_c | taxon) + (0 + T_c | taxon)`. Species without a woody or
#' non-woody class are left out. A group with fewer than `min_group_species`
#' species is not fitted; the trait is then fitted on the other group alone.
#'
#' With `location`, the model of the location prototype (#41,
#' `R/model/local_range.R`): species and location variances per growth form
#' and a species x location term, so that the spread among species at one
#' location can be separated from the shift between locations of the same
#' climate:
#' `... + diag(0 + group | taxon) + (0 + P_c | taxon) + (0 + T_c | taxon)
#'  + diag(0 + group | location) + (1 | taxon:location)`.
#' `data` must then carry that column (`add_location_key()`).
#'
#' With `spread_by` as well, the species variance also depends on climate:
#' species are binned on their climate centre (`n_bins` quantile classes of
#' the species within each growth form, `add_spread_class()`) and each bin x
#' growth form gets its own species variance,
#' `diag(0 + spread_class | taxon)` in place of `diag(0 + group | taxon)`.
#'
#' @param data `unified`.
#' @param groups woody lookup (`build_woody_lookup()`).
#' @param location NULL (the pooled model) or the name of a location column.
#' @param spread_by NULL, "precipitation" or "temperature" (needs `location`).
#' @return list(trait, groups (fitted), n_species (per group, before the
#'   minimum), formula, fit, data, fixed, spread_bins), or NULL if no group
#'   can be fitted.
fit_woodiness_model <- function(data, trait, groups, n_min = 5, min_group_species = 20,
                                location = NULL, spread_by = NULL, n_bins = 4) {
  prep <- .woodiness_data(data, trait, groups, n_min, min_group_species)
  if (is.null(prep)) return(NULL)
  sub <- prep$data
  fitted <- prep$groups
  n_species <- prep$n_species
  fixed <- prep$fixed
  slopes <- "+ (0 + log_annual_precipitation_c | taxon_name) + (0 + average_temperature_c | taxon_name)"
  spread_bins <- NULL
  if (!is.null(spread_by)) {
    stopifnot(!is.null(location))
    classed <- add_spread_class(sub, spread_by, n_bins)
    sub <- classed$data
    spread_bins <- classed$bins
  }
  form <- if (is.null(location)) {
    stats::as.formula(paste("log_value ~", fixed, "+ (1 | taxon_name)", slopes))
  } else {
    # per-group variances need both groups; with one, plain intercepts
    by_group <- if (length(fitted) > 1) "diag(0 + woodiness_group |" else "(1 |"
    species_term <- if (is.null(spread_by)) paste(by_group, "taxon_name)") else
      "diag(0 + spread_class | taxon_name)"
    stats::as.formula(paste(
      "log_value ~", fixed, "+", species_term, slopes,
      sprintf("+ %s %s) + (1 | taxon_name:%s)", by_group, location, location)))
  }
  message("== ", trait, " (woodiness model: ", paste(fitted, collapse = ", "),
          if (!is.null(location)) paste0("; location = ", location),
          if (!is.null(spread_by)) paste0("; spread by ", spread_by), ") ==")
  fit <- glmmTMB::glmmTMB(form, data = sub)
  b <- glmmTMB::fixef(fit)$cond
  se <- sqrt(diag(stats::vcov(fit)$cond))[names(b)]
  list(trait = trait, groups = fitted, n_species = n_species, formula = form,
       location = location, spread_by = spread_by, spread_bins = spread_bins,
       fit = fit, data = sub,
       fixed = tibble::tibble(term = names(b), estimate = unname(b), std.error = unname(se)))
}

# The records a woodiness model is fitted to, and its fixed-effect formula:
# one trait, species with a woody or non-woody class and more than `n_min`
# records, growth forms with at least `min_group_species` species. Shared by
# fit_woodiness_model() and the RTMB spread model (R/model/spread_rtmb.R).
.woodiness_data <- function(data, trait, groups, n_min = 5, min_group_species = 20) {
  sub <- data |>
    dplyr::filter(trait_name == trait) |>
    dplyr::inner_join(groups |>
                        dplyr::filter(woodiness_group %in% c("woody", "non_woody")) |>
                        dplyr::select(taxon_name, woodiness_group),
                      by = "taxon_name")
  keep <- sub |> dplyr::count(taxon_name) |> dplyr::filter(n > n_min) |> dplyr::pull(taxon_name)
  sub <- dplyr::filter(sub, taxon_name %in% keep)
  n_species <- dplyr::count(dplyr::distinct(sub, taxon_name, woodiness_group), woodiness_group)
  fitted <- intersect(c("woody", "non_woody"),
                      n_species$woodiness_group[n_species$n >= min_group_species])
  if (length(fitted) == 0) return(NULL)
  sub <- dplyr::filter(sub, woodiness_group %in% fitted) |>
    dplyr::mutate(woodiness_group = factor(woodiness_group, fitted))
  fixed <- if (length(fitted) > 1) {
    paste(c("0 + woodiness_group", paste0("woodiness_group:", woodiness_model_terms)),
          collapse = " + ")
  } else {
    paste(woodiness_model_terms, collapse = " + ")
  }
  list(data = sub, groups = fitted, n_species = n_species, fixed = fixed)
}

# Coefficient names of one group: intercept, then the four climate terms.
.group_coef_names <- function(model, group) {
  if (length(model$groups) > 1) {
    c(paste0("woodiness_group", group),
      paste0("woodiness_group", group, ":", woodiness_model_terms))
  } else {
    c("(Intercept)", woodiness_model_terms)
  }
}
