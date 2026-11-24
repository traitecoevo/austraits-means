# Location model prototype (#41): the community variance at one location.
#
# The community statistics are, for a location and from its climate, the
# woody and non-woody community mean and community range (terminology.md).
# The community variance is the spread of the species means of one growth
# form at one location. A model with species effects only pools every record
# from the same climate, so the spread among its species also carries how
# locations with the same climate differ. The location model
# (fit_woodiness_model(location = ...)) separates the two.
#
# For the species i of growth form g recorded at location l, the model's
# expected species mean there is
#   f_i + u_i + s_i'(x_il - m_i) + L_lg + w_il + e_il,
# with f_i the fixed-effect prediction, u_i the species effect, s_i the
# species' climate slopes, L_lg the location effect (shared by the species
# of g there, so it does not spread them), w_il the species x location
# effect, and e_il the mean of n_il residuals. The expected spread of the
# observed species means at l is therefore
#   var_i(f_i) + sd_u,g^2 + sd_w^2 + mean_i(sum_k sd_s,k^2 d_ik^2) + mean_i(sigma^2 / n_il)
# where d_ik is the species' mean centred climate at l. The community
# variance is that spread less the last (sampling) term.

#' Add a `location` key to record-level data.
#'
#' Coarse cells are the intended grain: a location then pools several studies
#' and more species, so the location effect is less likely to absorb the
#' similarity of a few species picked together by one study (which would
#' understate the community range).
#'
#' @param grain "cell" (a square cell of the community grid's equal-area
#'   projection, the default), "site" (`dataset_id` x `location_id`) or
#'   "coord" (exact coordinates).
#' @param cell_km cell size for `grain = "cell"`.
add_location_key <- function(data, grain = c("cell", "site", "coord"), cell_km = 50) {
  grain <- match.arg(grain)
  location <- switch(grain,
    cell = {
      settings <- community_layer_settings()
      settings$resolution <- cell_km * 1000
      cell <- lonlat_to_cell(data$decimal_longitude, data$decimal_latitude, settings)
      paste(cell$cell_x, cell$cell_y, sep = "_")
    },
    site = paste(data$dataset_id, data$location_id, sep = "|"),
    coord = paste(round(data$decimal_latitude, 4), round(data$decimal_longitude, 4))
  )
  # not mutate(location = location): data masking would keep an existing column
  data$location <- location
  data
}

#' Variance components of a glmmTMB fit, one row per grouping x term, plus
#' the residual.
#'
#' @param fit a glmmTMB fit, or a model list with `$fit`.
#' @return tibble(grouping, term, sd, variance).
tidy_varcomp <- function(fit) {
  if (is.list(fit) && !inherits(fit, "glmmTMB")) fit <- fit$fit
  vc <- glmmTMB::VarCorr(fit)$cond
  re <- purrr::imap(vc, function(m, nm) {
    sd <- attr(m, "stddev")
    tibble::tibble(grouping = sub("\\.[0-9]+$", "", nm), term = names(sd), sd = unname(sd))
  }) |> purrr::list_rbind()
  dplyr::bind_rows(re, tibble::tibble(grouping = "residual", term = "sigma",
                                      sd = stats::sigma(fit))) |>
    dplyr::mutate(variance = sd^2)
}

# Variance of one grouping x term, 0 if the model has no such term.
.vc <- function(vc, grouping, term) {
  v <- vc$variance[vc$grouping == grouping & vc$term == term]
  if (length(v) == 0) 0 else v[[1]]
}

# Species-effect variance of a growth form: per group in the location model,
# shared in the published model.
.vc_species <- function(vc, group) {
  v <- .vc(vc, "taxon_name", paste0("woodiness_group", group))
  if (v == 0) .vc(vc, "taxon_name", "(Intercept)") else v
}

#' Observed and model-implied community variance at each location.
#'
#' For each location x growth form with at least `min_species` species, the
#' variance of the species' observed mean log values there, against the
#' spread the model implies for those same species (header). Works for the
#' published model and the location model alike; `vc` can be replaced by
#' known values (simulation).
#'
#' @param model output of `fit_woodiness_model()`; its `data` must carry
#'   `location` (`add_location_key()`).
#' @param vc variance components (`tidy_varcomp()`).
#' @param noise_vc variance components whose residual gives the sampling
#'   term subtracted from the observed variance; the location model's, so
#'   both models are scored against the same target.
#' @return tibble(location, woodiness_group, n_species, observed_variance,
#'   sampling_variance, observed_community_variance, fixed_variance,
#'   species_variance, species_location_variance, slope_variance,
#'   model_community_variance).
benchmark_community_variance <- function(model, vc = tidy_varcomp(model),
                                         noise_vc = vc, min_species = 5) {
  d <- model$data
  stopifnot("location" %in% names(d))
  d$fixed_pred <- stats::predict(model$fit, newdata = d, re.form = NA,
                                 allow.new.levels = TRUE)
  sd_p <- .vc(vc, "taxon_name", "log_annual_precipitation_c")
  sd_t <- .vc(vc, "taxon_name", "average_temperature_c")
  sd_w <- .vc(vc, "taxon_name:location", "(Intercept)")
  sigma2 <- .vc(noise_vc, "residual", "sigma")

  d |>
    dplyr::group_by(location, woodiness_group, taxon_name) |>
    dplyr::summarise(.groups = "drop", n = dplyr::n(), y = mean(log_value),
                     f = mean(fixed_pred),
                     d_p = mean(log_annual_precipitation_c),
                     d_t = mean(average_temperature_c)) |>
    dplyr::group_by(location, woodiness_group) |>
    dplyr::filter(dplyr::n() >= min_species) |>
    dplyr::summarise(.groups = "drop",
                     n_species = dplyr::n(),
                     observed_variance = stats::var(y),
                     sampling_variance = mean(sigma2 / n),
                     fixed_variance = stats::var(f),
                     slope_variance = mean(sd_p * d_p^2 + sd_t * d_t^2)) |>
    dplyr::mutate(
      woodiness_group = as.character(woodiness_group),
      observed_community_variance = observed_variance - sampling_variance,
      species_variance = vapply(woodiness_group, function(g) .vc_species(vc, g), numeric(1)),
      species_location_variance = sd_w,
      model_community_variance = fixed_variance + species_variance +
        species_location_variance + slope_variance)
}

#' Summary of `benchmark_community_variance()` per growth form: means over
#' locations, reported as SDs (natural-log units).
summarise_community_benchmark <- function(bench) {
  bench |>
    dplyr::group_by(woodiness_group) |>
    dplyr::summarise(.groups = "drop",
                     n_locations = dplyr::n(),
                     observed_sd = sqrt(mean(observed_variance)),
                     observed_community_sd = sqrt(max(mean(observed_community_variance), 0)),
                     model_community_sd = sqrt(mean(model_community_variance)),
                     ratio = model_community_sd / observed_community_sd)
}

#' Each species' modelled mean on the log scale: the model's value at the
#' species' climate centre (centred climate 0), with its species effect, and
#' the location and species x location effects at zero.
#' @return tibble(taxon_name, woodiness_group, log_centre).
species_centre_values <- function(model) {
  d <- model$data |>
    dplyr::distinct(taxon_name, woodiness_group, log_annual_precipitation_mn,
                    average_temperature_mn) |>
    dplyr::mutate(log_annual_precipitation_c = 0, average_temperature_c = 0)
  fixed <- stats::predict(model$fit, newdata = d |> dplyr::mutate(location = NA),
                          re.form = NA, allow.new.levels = TRUE)
  re <- glmmTMB::ranef(model$fit)$cond$taxon_name
  re <- tibble::as_tibble(re, rownames = "taxon_name")
  d |>
    dplyr::left_join(re, by = "taxon_name") |>
    dplyr::mutate(
      u = if ("(Intercept)" %in% names(re)) .data[["(Intercept)"]] else
        ifelse(woodiness_group == "woody", .data[["woodiness_groupwoody"]],
               .data[["woodiness_groupnon_woody"]]),
      log_centre = fixed + dplyr::coalesce(u, 0)) |>
    dplyr::transmute(taxon_name, woodiness_group = as.character(woodiness_group), log_centre)
}

#' The location prototype for one trait: the published woodiness model and
#' the location model, compared (roadmap.md, Later work; #41).
#'
#' @return list(trait, cell_km, n_records, n_species, n_locations, seconds,
#'   converged, aic, varcomp_old, varcomp_new, benchmark (per growth form, both
#'   models), species_shift (per growth form: SD and max of the change in the
#'   log species mean, and its correlation)), or NULL if no group fits.
run_location_prototype <- function(data, trait, groups, cell_km = 50, min_species = 5) {
  data <- add_location_key(dplyr::filter(data, trait_name == trait), "cell", cell_km = cell_km)
  old <- fit_woodiness_model(data, trait, groups)
  if (is.null(old)) return(NULL)
  t0 <- Sys.time()
  new <- fit_woodiness_model(data, trait, groups, location = "location")
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  vc_old <- tidy_varcomp(old)
  vc_new <- tidy_varcomp(new)
  bench <- dplyr::bind_rows(
    location = summarise_community_benchmark(
      benchmark_community_variance(new, noise_vc = vc_new, min_species = min_species)),
    published = summarise_community_benchmark(
      benchmark_community_variance(old, noise_vc = vc_new, min_species = min_species)),
    .id = "model")
  shift <- dplyr::inner_join(species_centre_values(old), species_centre_values(new),
                             by = c("taxon_name", "woodiness_group"), suffix = c("_old", "_new")) |>
    dplyr::group_by(woodiness_group) |>
    dplyr::summarise(.groups = "drop", n_species = dplyr::n(),
                     sd_change = stats::sd(log_centre_new - log_centre_old),
                     max_abs_change = max(abs(log_centre_new - log_centre_old)),
                     correlation = stats::cor(log_centre_old, log_centre_new))
  list(trait = trait, cell_km = cell_km, n_records = nrow(new$data),
       n_species = dplyr::n_distinct(new$data$taxon_name),
       n_locations = dplyr::n_distinct(new$data$location),
       seconds = secs, converged = isTRUE(new$fit$sdr$pdHess),
       aic = c(published = stats::AIC(old$fit), location = stats::AIC(new$fit)),
       varcomp_old = vc_old, varcomp_new = vc_new, benchmark = bench, species_shift = shift)
}

#' Community range from a community mean and variance on the log scale: the
#' central `level` interval of the species values, back-transformed.
community_range <- function(mean_log, variance_log, level = 0.90) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  sd <- sqrt(variance_log)
  tibble::tibble(community_mean = exp(mean_log), community_sd = sd,
                 range_lower = exp(mean_log - z * sd), range_upper = exp(mean_log + z * sd))
}

#' Simulated records in the unified schema, for testing the location model.
#'
#' Locations have climate and a location effect L; species of two growth
#' forms have climate centres, effects u and records at a few locations near
#' their centre. With `sorting > 0` a species is more likely to be recorded
#' where L is close to its own u (species sort among locations by the same
#' unmeasured factor), which breaks the independence the model assumes.
#'
#' @return list(data, groups (woody lookup), truth (variance components in
#'   the `tidy_varcomp()` format)).
simulate_location_data <- function(n_species = 300, n_locations = 150,
                                   locations_per_species = 4, records = 3,
                                   sd_u = 0.4, sd_L = 0.4, sd_w = 0.15,
                                   sigma = 0.2, sorting = 0, seed = 1) {
  set.seed(seed)
  loc <- tibble::tibble(location = paste0("L", seq_len(n_locations)),
                        x_p = stats::rnorm(n_locations, 6.5, 0.6),
                        x_t = stats::rnorm(n_locations, 18, 4),
                        L = stats::rnorm(n_locations, 0, sd_L))
  sp <- tibble::tibble(taxon_name = paste0("sp", seq_len(n_species)),
                       woodiness_group = rep(c("woody", "non_woody"), length.out = n_species),
                       m_p = stats::rnorm(n_species, 6.5, 0.5),
                       m_t = stats::rnorm(n_species, 18, 3.5),
                       u = stats::rnorm(n_species, 0, sd_u))
  b <- list(woody = c(1.0, 0.30, -0.04, 0.05, -0.01),
            non_woody = c(0.5, 0.20, -0.02, 0.03, 0.00))
  pick <- function(i) {
    dist <- sqrt(((loc$x_p - sp$m_p[i]) / 0.6)^2 + ((loc$x_t - sp$m_t[i]) / 4)^2)
    w <- exp(-dist^2) * exp(-sorting * (loc$L - sp$u[i])^2 / (2 * sd_L^2))
    sample(seq_len(n_locations), locations_per_species, prob = w + 1e-9)
  }
  pairs <- purrr::map(seq_len(n_species), function(i) tibble::tibble(i = i, j = pick(i))) |>
    purrr::list_rbind() |>
    dplyr::mutate(w = stats::rnorm(dplyr::n(), 0, sd_w))
  data <- pairs[rep(seq_len(nrow(pairs)), each = records), ] |>
    dplyr::mutate(
      taxon_name = sp$taxon_name[i], woodiness_group = sp$woodiness_group[i],
      location = loc$location[j],
      log_annual_precipitation_mn = sp$m_p[i], average_temperature_mn = sp$m_t[i],
      log_annual_precipitation_c = loc$x_p[j] - sp$m_p[i],
      average_temperature_c = loc$x_t[j] - sp$m_t[i])
  coef <- do.call(rbind, b[data$woodiness_group])
  X <- cbind(1, data$log_annual_precipitation_mn, data$average_temperature_mn,
             data$log_annual_precipitation_c, data$average_temperature_c)
  data <- data |>
    dplyr::mutate(trait_name = "sim",
                  log_value = rowSums(X * coef) + sp$u[i] + loc$L[j] + w +
                    stats::rnorm(dplyr::n(), 0, sigma),
                  value = exp(log_value)) |>
    dplyr::select(-i, -j, -w, -woodiness_group)
  groups <- dplyr::select(sp, taxon_name, woodiness_group)
  truth <- tibble::tribble(
    ~grouping, ~term, ~sd,
    "taxon_name", "woodiness_groupwoody", sd_u,
    "taxon_name", "woodiness_groupnon_woody", sd_u,
    "location", "woodiness_groupwoody", sd_L,
    "location", "woodiness_groupnon_woody", sd_L,
    "taxon_name:location", "(Intercept)", sd_w,
    "residual", "sigma", sigma) |>
    dplyr::mutate(variance = sd^2)
  list(data = data, groups = groups, truth = truth)
}
