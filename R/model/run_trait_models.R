# Per-trait orchestration: fit the trait-climate model for each trait, extract
# species centres, and export them. Traits are independent, so the loop can run
# in parallel with furrr (multisession) when requested.

#' Fit + process + export one trait.
#'
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
