# Backend-agnostic fit for the within-group-centred trait-climate model.
#
# Two interchangeable backends share one formula (build_model_formula()):
#   - "tmb"  : glmmTMB, fast ML fit for iteration/testing.
#   - "brms" : full posterior via Stan/cmdstanr for the deliverable.
# Both return a standardized `trait_climate_fit` object so downstream code
# (species_centres(), reporting) is backend-blind.
#
# Bayesian guardrails (see the bayesian-workflow skill): explicit weakly
# informative priors, reproducible descriptive seed, 4 parallel chains,
# cmdstanr backend, non-centered parameterization (brms default for group
# effects). Convergence/PPC/LOO/priorsense checks are run around the fit in the
# analysis notebook, not hidden here; the raw fit is returned for that purpose.

#' Fit the trait-climate model.
#'
#' @param data record-level data with `<var>_mn`, `<var>_c` columns and `group`.
#' @param trait trait label (used in seed + stored on the object).
#' @param response response column.
#' @param climate_vars climate variable stems.
#' @param group grouping factor (species).
#' @param backend "tmb" or "brms".
#' @param n_min keep only groups with > n_min records.
#' @param methodology optional methodology grouping factor (round-3 hook).
#' @param ... extra args passed to the backend's fitting function.
#' @return object of class `trait_climate_fit`.
fit_trait_climate_model <- function(data,
                                    trait = "plant_height",
                                    response = "log_value",
                                    climate_vars = c("log_annual_precipitation", "average_temperature"),
                                    group = "taxon_name",
                                    backend = c("tmb", "brms"),
                                    n_min = 5,
                                    methodology = NULL,
                                    threads = NULL,
                                    file = NULL,
                                    ...) {
  backend <- match.arg(backend)

  # Keep groups with enough records to estimate a response
  keep <- data |>
    dplyr::count(.data[[group]]) |>
    dplyr::filter(n > n_min) |>
    dplyr::pull(.data[[group]])
  data_model <- data |> dplyr::filter(.data[[group]] %in% keep)

  form <- build_model_formula(response, climate_vars, group, methodology)

  # `threads`/`file` are BRMS-only (threading + on-disk fit cache); TMB ignores them.
  fit <- switch(backend,
    tmb  = .fit_tmb(form, data_model, ...),
    brms = .fit_brms(form, data_model, trait = trait,
                     threads = threads, file = file, ...)
  )

  structure(
    list(
      backend = backend,
      trait = trait,
      response = response,
      climate_vars = climate_vars,
      group = group,
      formula = form,
      fit = fit,
      data = data_model,
      fixed = .tidy_fixed(fit, backend),
      ranef = .tidy_ranef(fit, backend, group, climate_vars)
    ),
    class = "trait_climate_fit"
  )
}

# ---- backend fitters --------------------------------------------------------

.fit_tmb <- function(form, data, ...) {
  glmmTMB::glmmTMB(form, data = data, ...)
}

.fit_brms <- function(form, data, trait,
                      chains = 4, iter = 2000, warmup = 1000,
                      adapt_delta = 0.95, threads = NULL, file = NULL, ...) {
  # Weakly informative priors on the log scale; data dominate for large samples
  # (so the posterior is close to the TMB ML fit) while ruling out absurd values.
  priors <- c(
    brms::set_prior("normal(0, 5)", class = "b"),                # climate effects
    brms::set_prior("student_t(3, 0, 2.5)", class = "sd"),       # group-level SDs
    brms::set_prior("student_t(3, 0, 2.5)", class = "sigma")     # residual SD
  )

  # Within-chain threading (reduce_sum) for the large per-species fits.
  threading <- if (!is.null(threads) && threads > 1) brms::threading(threads) else NULL

  brms::brm(
    formula = form,
    data = data,
    family = stats::gaussian(),
    prior = priors,
    chains = chains, cores = chains,
    iter = iter, warmup = warmup,
    control = list(adapt_delta = adapt_delta),
    backend = "cmdstanr",
    threads = threading,
    seed = sum(utf8ToInt(paste0("trait-climate-", trait))),
    # On-disk cache: reuse the fit unless the model/data change.
    file = file,
    file_refit = if (is.null(file)) "never" else "on_change",
    ...
  )
}

#' Convergence diagnostics for a fitted brms model (efficient at scale).
#'
#' Summarises divergences and the worst R-hat / ESS over the *reported*
#' parameters (population effects, group SDs, sigma) via `summary()` — it does
#' not enumerate the thousands of individual random effects.
brms_diagnostics <- function(fit) {
  np <- brms::nuts_params(fit)
  ndiv <- sum(np$Value[np$Parameter == "divergent__"])

  s <- summary(fit)
  grab <- function(tab, col) if (!is.null(tab) && nrow(tab)) tab[, col] else numeric(0)
  pick <- function(col) c(
    grab(s$fixed, col), grab(s$spec_pars, col),
    unlist(lapply(s$random, grab, col), use.names = FALSE)
  )

  tibble::tibble(
    divergences = ndiv,
    max_rhat = max(pick("Rhat"), na.rm = TRUE),
    min_ess_bulk = min(pick("Bulk_ESS"), na.rm = TRUE),
    min_ess_tail = min(pick("Tail_ESS"), na.rm = TRUE)
  )
}

# ---- standardized extractors ------------------------------------------------

# Normalize a random/fixed term name to a standardized column name.
.std_term <- function(term, climate_vars) {
  term <- ifelse(term %in% c("(Intercept)", "Intercept"), "g_0", term)
  ifelse(term %in% paste0(climate_vars, "_c"), paste0("g_", term), term)
}

.tidy_fixed <- function(fit, backend) {
  if (backend == "tmb") {
    broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) |>
      dplyr::transmute(term, estimate, std.error, conf.low, conf.high)
  } else {
    fx <- brms::fixef(fit, probs = c(0.03, 0.97))
    terms <- rownames(fx)
    terms[terms == "Intercept"] <- "(Intercept)"  # align with TMB naming
    tibble::tibble(
      term = terms,
      estimate = fx[, "Estimate"],
      std.error = fx[, "Est.Error"],
      conf.low = fx[, 3],
      conf.high = fx[, 4]
    )
  }
}

.tidy_ranef <- function(fit, backend, group, climate_vars) {
  if (backend == "tmb") {
    re <- glmmTMB::ranef(fit)$cond[[group]]
    tib <- tibble::as_tibble(re, rownames = "taxon_name")
  } else {
    arr <- brms::ranef(fit)[[group]]      # [levels, stat, term]
    terms <- dimnames(arr)[[3]]
    tib <- tibble::tibble(taxon_name = dimnames(arr)[[1]])
    for (tm in terms) tib[[tm]] <- arr[, "Estimate", tm]
  }
  # standardize column names: intercept -> g_0, slope <v>_c -> g_<v>_c
  nm <- names(tib)
  nm[nm != "taxon_name"] <- .std_term(nm[nm != "taxon_name"], climate_vars)
  names(tib) <- nm
  tib
}
