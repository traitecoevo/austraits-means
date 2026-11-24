# Per-trait orchestration: fit the trait-climate model for each trait, extract
# species centres, and export them. Traits are independent, so the loop can run
# in parallel with furrr (multisession) when requested.

#' Fit + process + export one trait (the pooled model).
#'
#' @param out_dir where to write the trait's species centres; NULL writes
#'   nothing.
#' @param fits_dir if set and backend = "brms", cache each Stan fit to
#'   `<fits_dir>/<trait>_brms.rds` (reused unless model/data change).
#' @return a `trait_climate_fit` (centres attached as `$centres`; for brms,
#'   convergence diagnostics as `$diagnostics`).
run_one_trait <- function(data, trait, backend = "tmb", n_min = 5,
                          out_dir = "outputs/traits", fits_dir = NULL, ...) {
  sub <- dplyr::filter(data, trait_name == trait)

  file <- NULL
  if (!is.null(fits_dir) && backend == "brms") {
    dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)
    file <- file.path(fits_dir, paste0(clean_species_filename(trait), "_brms"))
  }

  fit <- fit_trait_climate_model(sub, trait = trait, backend = backend,
                                 n_min = n_min, file = file, ...)
  fit$centres <- species_centres(fit)
  # the pipeline passes out_dir = NULL: the working exports come from the
  # location model (export_location_centres()), and the pooled model is a
  # diagnostic that writes no files
  if (!is.null(out_dir))
    export_species_centres(fit$centres, trait, out_dir,
                           unit = sub$unit[1], backend = backend, data = fit$data)

  if (backend == "brms") {
    fit$diagnostics <- brms_diagnostics(fit$fit)
    d <- fit$diagnostics
    message(sprintf("  diagnostics: divergences %d | max Rhat %.3f | min bulk-ESS %.0f",
                    d$divergences, d$max_rhat, d$min_ess_bulk))
  }
  fit
}

#' Safe single-trait fit for use as a `targets` branch.
#'
#' Thin wrapper around `run_one_trait()` that returns `NULL` (with a message)
#' instead of erroring, so one un-fittable trait does not abort `tar_make()`.
#' Filter `NULL`s downstream (e.g. `purrr::compact()`).
fit_trait_target <- function(data, trait, backend = "tmb", n_min = 5,
                             out_dir = "outputs/traits", fits_dir = NULL, ...) {
  message("== ", trait, " (", backend, ") ==")
  tryCatch(
    run_one_trait(data, trait, backend = backend, n_min = n_min,
                  out_dir = out_dir, fits_dir = fits_dir, ...),
    error = function(e) {
      message("  failed: ", conditionMessage(e)); NULL
    }
  )
}

#' Fit all (or selected) traits and export per-trait species centres.
#'
#' @param data unified analysis dataset (output of `build_dataset()`).
#' @param traits character vector; NULL = all traits in `data`.
#' @param backend "tmb" or "brms".
#' @param n_min minimum records per species.
#' @param out_dir output directory for the exports.
#' @param parallel run traits concurrently with furrr (set up a plan first, e.g.
#'   `setup_parallel()`); otherwise sequential. For brms, prefer sequential with
#'   per-fit threading, since each fit already uses all chains/threads.
#' @param fits_dir if set (brms), cache each Stan fit under this directory.
#' @param ... passed to `fit_trait_climate_model()` (e.g. `threads`, `iter`).
#' @return named list of `trait_climate_fit` objects (one per trait).
run_trait_models <- function(data, traits = NULL, backend = "tmb", n_min = 5,
                             out_dir = "outputs/traits", parallel = FALSE,
                             fits_dir = NULL, ...) {
  if (is.null(traits)) traits <- sort(unique(data$trait_name))

  run <- function(tr) {
    message("== ", tr, " (", backend, ") ==")
    tryCatch(
      run_one_trait(data, tr, backend = backend, n_min = n_min,
                    out_dir = out_dir, fits_dir = fits_dir, ...),
      error = function(e) {
        message("  failed: ", conditionMessage(e)); NULL
      }
    )
  }

  fits <- if (parallel) {
    furrr::future_map(traits, run, .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("glmmTMB", "brms", "dplyr", "tibble", "tidyr", "broom.mixed")
    ))
  } else {
    lapply(traits, run)
  }
  stats::setNames(fits, traits)
}

#' One-row-per-trait summary of fitted fixed effects (between/within slopes).
summarise_trait_fits <- function(fits) {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  purrr::imap_dfr(fits, function(fit, tr) {
    fx <- fit$fixed |>
      dplyr::select(term, estimate) |>
      tidyr::pivot_wider(names_from = term, values_from = estimate)
    tibble::tibble(trait = tr, n_species = nrow(fit$ranef)) |>
      dplyr::bind_cols(fx)
  })
}

#' One-row-per-trait BRMS convergence summary (from `$diagnostics`).
summarise_trait_diagnostics <- function(fits) {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  fits <- fits[vapply(fits, \(f) !is.null(f$diagnostics), logical(1))]
  purrr::imap_dfr(fits, function(fit, tr) {
    tibble::tibble(trait = tr, n_species = nrow(fit$ranef)) |>
      dplyr::bind_cols(fit$diagnostics)
  })
}

# ---- the location model (one model for the modelled means and the community
# statistics; roadmap.md, 25 November 2026) ---------------------------------

#' Cell size (km) of the locations in the location model: the
#' `add_location_key(type = "cell")` grain of the products.
location_cell_km <- function() 50

#' Fit the location model for one trait, as a `targets` branch: the RTMB
#' spread model free (log species SD linear in the climate centre), then with
#' its climate slopes shrunk (fit_spread_shrunk()), with every species'
#' modelled mean at its climate centre.
#'
#' A growth form the model cannot fit is dropped (roadmap.md, 25 November
#' 2026), not propped up by borrowing from the other: if either fit fails or
#' does not converge (including a Hessian that is not positive definite), the
#' smallest growth form is dropped and the trait refitted without it. NULL,
#' with a message, if no growth form can be fitted.
#'
#' @param data `unified` with a `location` column (add_location_key()).
#' @param groups woody lookup.
#' @param min_group_species species a growth form needs to be tried.
#' @return the shrunk `rtmb_spread_fit`, with `backend = "rtmb"`, `centres`
#'   in the schema of species_centres() (location_species_centres()), `free`
#'   (the free fit's AIC, coefficients and convergence) and `dropped_groups`
#'   (growth forms dropped because the fit failed, with the reason).
fit_location_target <- function(data, trait, groups, n_min = 5, min_group_species = 5) {
  message("== ", trait, " (location model, RTMB) ==")
  dropped <- character(0)
  repeat {
    res <- tryCatch({
      free <- fit_spread_rtmb(data, trait, groups, n_min = n_min,
                              min_group_species = min_group_species)
      if (is.null(free)) return(NULL)
      if (!free$converged) stop("free fit did not converge (", free$message, ")")
      m <- fit_spread_shrunk(data, trait, groups, free = free, n_min = n_min,
                             min_group_species = min_group_species, centres = TRUE)
      if (!m$converged) stop("shrunk fit did not converge (", m$message, ")")
      m$free <- list(aic = free$aic, theta = free$spread$theta, cov = free$spread$cov,
                     converged = free$converged, seconds = free$seconds)
      m
    }, error = function(e) e)
    if (!inherits(res, "error")) break
    prep <- .woodiness_data(data, trait, groups, n_min, min_group_species)
    if (is.null(prep) || length(prep$groups) < 2) {
      message("  failed: ", conditionMessage(res))
      return(NULL)
    }
    n <- prep$n_species$n[match(prep$groups, prep$n_species$woodiness_group)]
    smallest <- prep$groups[which.min(n)]
    dropped <- c(dropped, stats::setNames(conditionMessage(res),
                                          paste0(smallest, " (", min(n), " species)")))
    message("  dropping ", names(dropped)[length(dropped)], ": ", conditionMessage(res))
    min_group_species <- min(n) + 1
  }
  res$backend <- "rtmb"
  res$centres <- location_species_centres(res)
  res$dropped_groups <- dropped
  res
}

#' Fixed effects of a location-model fit, one block per growth form, with the
#' pooled model's term names ("(Intercept)", "<var>_mn", "<var>_c").
#'
#' @return tibble(trait, woodiness_group, term, estimate, std.error,
#'   conf.low, conf.high), 95% intervals.
location_fixed_table <- function(m) {
  purrr::map(m$groups, function(g) {
    ft <- m$fixed_table[match(.group_coef_names(m, g), m$fixed_table$term), ]
    ft$term <- c("(Intercept)", woodiness_model_terms)
    dplyr::mutate(ft, trait = m$trait, woodiness_group = g, .before = 1)
  }) |> purrr::list_rbind()
}

#' The fixed-effect table of either model: a location-model fit's
#' (location_fixed_table()) or a pooled fit's `$fixed` (with `trait`).
fit_fixed_table <- function(f) {
  if (inherits(f, "rtmb_spread_fit")) location_fixed_table(f) else
    dplyr::mutate(f$fixed, trait = f$trait, .before = 1)
}

#' One row per trait x growth form: the location model's fixed effects and
#' the number of species (the `trait_summary` target).
summarise_location_fits <- function(fits) {
  purrr::map(purrr::compact(fits), function(m) {
    n_sp <- dplyr::distinct(m$centres, taxon_name, woodiness_group) |>
      dplyr::mutate(woodiness_group = as.character(woodiness_group)) |>
      dplyr::count(woodiness_group, name = "n_species")
    location_fixed_table(m) |>
      dplyr::select(trait, woodiness_group, term, estimate) |>
      tidyr::pivot_wider(names_from = term, values_from = estimate) |>
      dplyr::left_join(n_sp, by = "woodiness_group") |>
      dplyr::relocate(n_species, .after = woodiness_group)
  }) |> purrr::list_rbind()
}
