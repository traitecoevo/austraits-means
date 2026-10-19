# Modelled community layers (#28): a fitted trait-climate model's net
# average trait value for a cell, averaging over between- and within-species
# effects, for woody and non-woody species separately. No species means are
# used.
#
# A record at climate x from a species centred at m has expected log trait
#   b0 + b_between * m + b_within * (x - m)
# (species random effects average to zero over species). Averaging over the
# species present in a cell gives
#   b0 + b_within * x + (b_between - b_within) * mean(m),
# so the cell value needs only the cell's climate x and the mean climate
# centre of the species found there (GBIF presence and the species medoids),
# for each climate variable. Species found in the cell need no trait data.
#
# The model used here gives woody and non-woody species their own intercept
# and climate slopes (fit_woodiness_model()), so each group's layer uses its
# own coefficients and its own species' mean centre; the layer for all species
# is the mix of the two, weighted by the number of species of each group in
# the cell.

woodiness_model_terms <- c("log_annual_precipitation_mn", "average_temperature_mn",
                           "log_annual_precipitation_c", "average_temperature_c")

#' Trait-climate model with woodiness, for the modelled community layers.
#'
#' The model of `fit_trait_climate_model()` (species with more than `n_min`
#' records), with a separate intercept and climate slopes for woody and
#' non-woody species and the species random effects shared:
#' `log_value ~ 0 + group + group:(P_mn + T_mn + P_c + T_c) + (1 | taxon) +
#' (0 + P_c | taxon) + (0 + T_c | taxon)`. Species without a woody or
#' non-woody class are left out. A group with fewer than `min_group_species`
#' species is not fitted; the trait is then fitted on the other group alone.
#'
#' @param data `unified`.
#' @param groups woody lookup (`build_woody_lookup()`).
#' @return list(trait, groups (fitted), n_species (per group, before the
#'   minimum), formula, fit, data, fixed), or NULL if no group can be fitted.
fit_woodiness_model <- function(data, trait, groups, n_min = 5, min_group_species = 20) {
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
  form <- stats::as.formula(paste(
    "log_value ~", fixed, "+ (1 | taxon_name)",
    "+ (0 + log_annual_precipitation_c | taxon_name) + (0 + average_temperature_c | taxon_name)"))
  message("== ", trait, " (woodiness model: ", paste(fitted, collapse = ", "), ") ==")
  fit <- glmmTMB::glmmTMB(form, data = sub)
  b <- glmmTMB::fixef(fit)$cond
  se <- sqrt(diag(stats::vcov(fit)$cond))[names(b)]
  list(trait = trait, groups = fitted, n_species = n_species, formula = form,
       fit = fit, data = sub,
       fixed = tibble::tibble(term = names(b), estimate = unname(b), std.error = unname(se)))
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

#' Mean climate of each grid cell (WorldClim averaged onto the layer grid,
#' cells touching Australian land only).
#'
#' @param rasters climate rasters (`load_climate_rasters()`).
#' @return tibble(cell_x, cell_y, annual_precipitation,
#'   log_annual_precipitation, average_temperature).
cell_climate <- function(rasters, settings = community_layer_settings()) {
  template <- community_grid(settings)
  clim <- terra::project(rasters, template, method = "average")
  names(clim) <- c("average_temperature", "annual_precipitation")
  australia <- terra::vect(sf::st_transform(ozmaps::ozmap_country, settings$crs))
  clim <- terra::mask(clim, australia, touches = TRUE)
  terra::as.data.frame(clim, xy = TRUE, na.rm = TRUE) |>
    tibble::as_tibble() |>
    dplyr::filter(annual_precipitation > 0) |>
    dplyr::transmute(
      cell_x = as.integer(floor(x / settings$resolution)),
      cell_y = as.integer(floor(y / settings$resolution)),
      annual_precipitation,
      log_annual_precipitation = log(annual_precipitation),
      average_temperature
    )
}

#' Mean climate centre of the woody and of the non-woody species present in
#' each cell.
#'
#' @param presence output of `gbif_presence_cells()` (GBIF names x cells).
#' @param medoids_path medoid parquet, keyed by GBIF name.
#' @param name_map,groups to place each GBIF species in woody / non-woody;
#'   species without a class are left out.
#' @return tibble(cell_x, cell_y, woodiness_group, n_species,
#'   centre_log_annual_precipitation, centre_average_temperature).
cell_species_centres <- function(presence, medoids_path, name_map, groups) {
  centres <- arrow::read_parquet(medoids_path) |>
    dplyr::filter(!is.na(annual_precipitation), annual_precipitation > 0) |>
    dplyr::transmute(gbif_taxon_name = taxon_name,
                     log_p = log(annual_precipitation), t = average_temperature)
  grp <- name_map |>
    dplyr::inner_join(dplyr::filter(groups, woodiness_group %in% c("woody", "non_woody")) |>
                        dplyr::select(taxon_name, woodiness_group), by = "taxon_name") |>
    dplyr::distinct(gbif_taxon_name, woodiness_group)
  presence |>
    dplyr::inner_join(centres, by = "gbif_taxon_name") |>
    dplyr::inner_join(grp, by = "gbif_taxon_name") |>
    dplyr::group_by(cell_x, cell_y, woodiness_group) |>
    dplyr::summarise(.groups = "drop",
                     n_species = dplyr::n(),
                     centre_log_annual_precipitation = mean(log_p),
                     centre_average_temperature = mean(t))
}

#' The woodiness model's net average per trait x cell x group.
#'
#' Woody and non-woody use their own coefficients and species centres; `all`
#' is their mix weighted by species counts (the mean over the species of both
#' groups present), and is empty for a trait whose model could not fit both
#' groups.
#'
#' @param models list of `fit_woodiness_model()` results.
#' @param climate output of `cell_climate()`.
#' @param centres output of `cell_species_centres()`.
#' @param level interval mass (delta method on the log scale, fixed-effect
#'   uncertainty only).
#' @return long tibble: trait, cell, group, `n_species`, cell climate, mean
#'   species centre, `community_mean` with `lower`/`upper`, `extrapolated`.
community_layers_modelled <- function(models, climate, centres, level = 0.94,
                                      settings = community_layer_settings()) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  cells <- centres |>
    dplyr::inner_join(climate, by = c("cell_x", "cell_y")) |>
    tidyr::pivot_wider(names_from = woodiness_group,
                       values_from = c(n_species, centre_log_annual_precipitation,
                                       centre_average_temperature),
                       values_fill = list(n_species = 0L))
  for (g in c("woody", "non_woody")) {
    for (col in paste0(c("n_species_", "centre_log_annual_precipitation_",
                         "centre_average_temperature_"), g)) {
      if (!col %in% names(cells)) cells[[col]] <- if (startsWith(col, "n_")) 0L else NA_real_
    }
  }

  purrr::map_dfr(purrr::compact(models), function(m) {
    b <- glmmTMB::fixef(m$fit)$cond
    V <- stats::vcov(m$fit)$cond[names(b), names(b)]
    zero <- matrix(0, nrow(cells), length(b), dimnames = list(NULL, names(b)))

    # one group's design matrix and extrapolation flag
    group_rows <- function(g) {
      cP <- cells[[paste0("centre_log_annual_precipitation_", g)]]
      cT <- cells[[paste0("centre_average_temperature_", g)]]
      X <- zero
      if (!g %in% m$groups) return(list(X = X, fitted = FALSE, extrapolated = NA))
      nm <- .group_coef_names(m, g)
      X[, nm[1]] <- 1
      X[, nm[2]] <- cP
      X[, nm[3]] <- cT
      X[, nm[4]] <- cells$log_annual_precipitation - cP
      X[, nm[5]] <- cells$average_temperature - cT
      X[is.na(X)] <- 0
      d <- m$data[m$data$woodiness_group == g, ]
      outside <- function(v, var) v < min(d[[var]]) | v > max(d[[var]])
      ext <- outside(cP, "log_annual_precipitation_mn") |
        outside(cT, "average_temperature_mn") |
        outside(X[, nm[4]], "log_annual_precipitation_c") |
        outside(X[, nm[5]], "average_temperature_c")
      list(X = X, fitted = TRUE, extrapolated = ext)
    }
    rows <- list(woody = group_rows("woody"), non_woody = group_rows("non_woody"))

    predict_group <- function(X, n, fitted, extrapolated, cP, cT, group) {
      eta <- as.numeric(X %*% b)
      se <- sqrt(rowSums((X %*% V) * X))
      ok <- fitted & n >= settings$min_species
      tibble::tibble(
        trait_name = m$trait, cell_x = cells$cell_x, cell_y = cells$cell_y,
        woodiness_group = group, n_species = as.integer(n),
        annual_precipitation = cells$annual_precipitation,
        log_annual_precipitation = cells$log_annual_precipitation,
        average_temperature = cells$average_temperature,
        centre_log_annual_precipitation = cP, centre_average_temperature = cT,
        community_mean = dplyr::if_else(ok, exp(eta), NA_real_),
        lower = dplyr::if_else(ok, exp(eta - z * se), NA_real_),
        upper = dplyr::if_else(ok, exp(eta + z * se), NA_real_),
        extrapolated = dplyr::if_else(ok, extrapolated, NA)
      ) |>
        dplyr::filter(n_species > 0)
    }

    n_w <- cells$n_species_woody
    n_n <- cells$n_species_non_woody
    n_all <- n_w + n_n
    wavg <- function(a, b) (n_w * dplyr::coalesce(a, 0) + n_n * dplyr::coalesce(b, 0)) / n_all
    X_all <- (n_w * rows$woody$X + n_n * rows$non_woody$X) / n_all
    # only where both groups are fitted, so 'all' never silently means woody
    all_fitted <- rows$woody$fitted & rows$non_woody$fitted
    ext_all <- (n_w > 0 & dplyr::coalesce(rows$woody$extrapolated, FALSE)) |
      (n_n > 0 & dplyr::coalesce(rows$non_woody$extrapolated, FALSE))

    dplyr::bind_rows(
      predict_group(rows$woody$X, n_w, rows$woody$fitted, rows$woody$extrapolated,
                    cells$centre_log_annual_precipitation_woody,
                    cells$centre_average_temperature_woody, "woody"),
      predict_group(rows$non_woody$X, n_n, rows$non_woody$fitted, rows$non_woody$extrapolated,
                    cells$centre_log_annual_precipitation_non_woody,
                    cells$centre_average_temperature_non_woody, "non_woody"),
      predict_group(X_all, n_all, all_fitted, ext_all,
                    wavg(cells$centre_log_annual_precipitation_woody,
                         cells$centre_log_annual_precipitation_non_woody),
                    wavg(cells$centre_average_temperature_woody,
                         cells$centre_average_temperature_non_woody), "all")
    )
  }) |>
    add_cell_centres(settings) |>
    dplyr::arrange(trait_name, woodiness_group, cell_y, cell_x)
}

#' Fixed effects of the woodiness models, one row per trait x term.
woodiness_model_summary <- function(models) {
  purrr::map_dfr(purrr::compact(models), function(m) {
    n <- tidyr::pivot_wider(m$n_species, names_from = woodiness_group, values_from = n,
                            names_prefix = "n_species_")
    dplyr::mutate(m$fixed, trait_name = m$trait, .before = 1,
                  groups_fitted = paste(m$groups, collapse = ", ")) |>
      dplyr::bind_cols(n)
  })
}

community_modelled_metadata <- function(data_basis, native_only = FALSE) {
  community_metadata(
    data_basis, paste0("modelled", if (native_only) ", native species"),
    paste("For each cell and group, the fitted trait-climate model's net average trait",
          "value over the species present: exp(b0 + b_within * x + (b_between - b_within) * m),",
          "where x is the cell's climate and m the mean climate centre (medoid) of the",
          "group's GBIF species recorded there, for log annual precipitation and mean annual",
          "temperature, and the coefficients are the group's own. Averages over between-",
          "and within-species effects; uses no species means. 'all' is the woody and",
          "non-woody values mixed in proportion to their species in the cell (log scale)."),
    groups = paste("woody and non-woody by the rule of #19, each with its own model",
                   "coefficients; 'all' is their species-weighted mix. A group the model",
                   "could not fit (fewer than 20 species with enough records) has no values,",
                   "nor then does 'all'"),
    extra = c(if (native_only) native_only_note(), list(
      model = paste("log_value ~ 0 + group + group:(P_mn + T_mn + P_c + T_c) + (1 | taxon)",
                    "+ (0 + P_c | taxon) + (0 + T_c | taxon), glmmTMB; group = woody / non-woody;",
                    "species with more than 5 records and a climate centre"),
      climate = "WorldClim 2.1 (bio1, bio12) averaged onto the grid",
      interval = "lower, upper: 94% interval, delta method on the log scale, fixed-effect uncertainty only",
      columns = list(n_species = "GBIF species of the group in the cell with a climate centre",
                     extrapolated = paste("TRUE where the group's mean centre or the cell-minus-centre",
                                          "difference lies outside the range its model was fitted on")))))
}

#' Map of the modelled layer beside the observed presence layer, woody and
#' non-woody, on one colour scale. Extrapolated cells are shown grey.
plot_community_modelled <- function(modelled, observed, trait, unit = NULL,
                                    species_label = "") {
  grp <- c(woody = "Woody", non_woody = "Non-woody")
  prep <- function(d) {
    dplyr::filter(d, trait_name == trait, woodiness_group %in% names(grp)) |>
      dplyr::mutate(woodiness_group = factor(grp[woodiness_group], grp))
  }
  m <- prep(modelled) |>
    dplyr::mutate(community_mean = dplyr::if_else(dplyr::coalesce(extrapolated, FALSE),
                                                  NA_real_, community_mean))
  o <- prep(observed)
  limits <- range(c(m$community_mean, o$community_mean), na.rm = TRUE)
  unit_label <- if (is.null(unit) || is.na(unit)) "" else paste0(" (", unit, ")")
  patchwork::wrap_plots(
    map_cells(m, "community_mean", layer_colours$mean, "Modelled: net average from the woodiness model",
              "log10", facet = ~woodiness_group, limits = limits),
    map_cells(o, "community_mean", layer_colours$mean,
              "Observed: species present (GBIF) x observed species means",
              "log10", facet = ~woodiness_group, limits = limits),
    ncol = 2, guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = paste0(trait, if (nzchar(species_label)) paste0(", ", species_label),
                     ", community mean", unit_label),
      subtitle = paste("10 km equal-area grid; grey: fewer than 10 species, or (modelled) outside",
                       "the fitted range, or a group the model could not fit"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                             plot.subtitle = ggplot2::element_text(size = 9, colour = "#6b6a66")))
}

#' Save one modelled-vs-observed map per trait; returns the file paths.
write_modelled_maps <- function(modelled, observed, units, dir, species_label = "") {
  unlist(lapply(sort(unique(modelled$trait_name)), function(tr) {
    p <- plot_community_modelled(modelled, observed, tr, units$unit[units$trait_name == tr],
                                 species_label)
    save_figure(p, file.path(dir, "figures", paste0(tr, "_modelled")), width = 9, height = 8)
  }), use.names = FALSE)
}
