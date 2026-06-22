# Predicted trait value at each species' climate centre (the medoid).
#
# With Option-B centring the medoid is exactly the point where every `<var>_c`
# is 0, so the species-mean trait at the centre is the model prediction at
# newdata with `_c = 0` and `_mn` = medoid. The within-species slopes drop out
# (multiplied by 0), so the centre depends only on the intercept (incl. the
# species random intercept) and the between-species mean terms.
#
# Backend-agnostic: TMB uses glmmTMB predict (delta-method SE); brms uses the
# posterior predictive expectation (94% credible interval).

#' @param fit_obj a `trait_climate_fit`.
#' @param level interval mass (default 0.94, per the Bayesian-workflow default).
#' @return tibble with one row per species: centre climate, predicted trait
#'   (original + log scale) with interval, sample size, and random effects.
species_centres <- function(fit_obj, level = 0.94) {
  stopifnot(inherits(fit_obj, "trait_climate_fit"))
  cv <- fit_obj$climate_vars
  group <- fit_obj$group
  mn_cols <- paste0(cv, "_mn")
  c_cols <- paste0(cv, "_c")

  # one row per species: medoid (_mn) held, deviations (_c) set to 0
  newdata <- fit_obj$data |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group, mn_cols))))
  for (cc in c_cols) newdata[[cc]] <- 0

  pred <- .predict_centre(fit_obj, newdata, level)

  # per-species context: sample size, taxonomy, medoid climate
  context <- fit_obj$data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group))) |>
    dplyr::summarise(
      .groups = "drop",
      n = dplyr::n(),
      dplyr::across(dplyr::any_of(c("family", "genus")), dplyr::first),
      dplyr::across(dplyr::any_of(c("annual_precipitation_centre",
                                    "average_temperature_centre",
                                    "log_annual_precipitation_centre")), dplyr::first)
    )

  newdata |>
    dplyr::select(dplyr::all_of(c(group, mn_cols))) |>
    dplyr::bind_cols(pred) |>
    dplyr::left_join(context, by = group) |>
    dplyr::left_join(fit_obj$ranef, by = group) |>
    dplyr::relocate(dplyr::any_of(c("family", "genus", "n")), .after = dplyr::all_of(group))
}

# Predict the centre value on the response (= log trait) scale, then exponentiate.
# Returns a tibble: log_estimate, log_conf.low/high and original-scale estimate,
# conf.low/high.
.predict_centre <- function(fit_obj, newdata, level) {
  a <- (1 - level) / 2
  probs <- c(a, 1 - a)

  if (fit_obj$backend == "tmb") {
    p <- stats::predict(fit_obj$fit, newdata = newdata,
                        se.fit = TRUE, re.form = NULL, type = "link")
    z <- stats::qnorm(1 - a)
    log_est <- p$fit
    log_lo <- p$fit - z * p$se.fit
    log_hi <- p$fit + z * p$se.fit
  } else {
    epred <- brms::posterior_epred(fit_obj$fit, newdata = newdata,
                                   re_formula = NULL)
    log_est <- colMeans(epred)
    q <- apply(epred, 2, stats::quantile, probs = probs)
    log_lo <- q[1, ]
    log_hi <- q[2, ]
  }

  tibble::tibble(
    log_estimate = as.numeric(log_est),
    log_conf.low = as.numeric(log_lo),
    log_conf.high = as.numeric(log_hi),
    estimate = exp(as.numeric(log_est)),
    conf.low = exp(as.numeric(log_lo)),
    conf.high = exp(as.numeric(log_hi))
  )
}
