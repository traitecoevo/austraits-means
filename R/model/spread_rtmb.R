# Location model with the between-species spread a continuous function of
# climate (#41), in RTMB: glmmTMB cannot make a random-effect variance depend
# on a covariate.
#
# Same model as fit_woodiness_model(location = ...) except the species
# effect: for species s of growth form g with climate centre (P_s, T_s),
# standardised to (p, t),
#   u_s ~ N(0, sd_u(s)^2),
#   log sd_u(s) = a_g + b1_g p + b2_g t + b3_g p^2 + b4_g t^2.
# The rest as before: growth-form intercepts and climate slopes (fixed),
# uncorrelated species random slopes on the centred climate, a location
# effect per growth form, a species x location effect and the residual.
# Random effects are integrated out by the Laplace approximation, as in
# glmmTMB, so for this Gaussian model the likelihood is exact and AICs are
# comparable with the glmmTMB fits on the same records.
#
# With `spread_terms = character(0)`, log sd_u(s) = a_g: the glmmTMB location
# model, which checks this code (the likelihoods must agree).

.spread_term_cols <- function(p, t, terms) {
  all <- cbind(p = p, t = t, p2 = p^2, t2 = t^2)
  all[, terms, drop = FALSE]
}

#' Fit the location model with the species spread a smooth function of the
#' species' climate centre.
#'
#' @param data records with a `location` column (`add_location_key()`).
#' @param spread_terms terms of log sd_u in the standardised centres: any of
#'   "p", "t", "p2", "t2" (precipitation, temperature, their squares).
#' @return an `rtmb_spread_fit`: list(trait, groups, data, fixed_pred, vc,
#'   species_variance, location_effects, spread (coefficients and their
#'   covariance), centre_scaling, logLik, df, aic, converged, seconds).
fit_spread_rtmb <- function(data, trait, groups, n_min = 5, min_group_species = 20,
                            spread_terms = c("p", "t", "p2", "t2")) {
  prep <- .woodiness_data(data, trait, groups, n_min, min_group_species)
  if (is.null(prep)) return(NULL)
  d <- prep$data
  stopifnot("location" %in% names(d))
  G <- length(prep$groups)
  X <- stats::model.matrix(stats::as.formula(paste("~", prep$fixed)), d)

  species <- d |>
    dplyr::distinct(taxon_name, woodiness_group, log_annual_precipitation_mn,
                    average_temperature_mn)
  scaling <- list(p_mean = mean(species$log_annual_precipitation_mn),
                  p_sd = stats::sd(species$log_annual_precipitation_mn),
                  t_mean = mean(species$average_temperature_mn),
                  t_sd = stats::sd(species$average_temperature_mn))
  p <- (species$log_annual_precipitation_mn - scaling$p_mean) / scaling$p_sd
  t <- (species$average_temperature_mn - scaling$t_mean) / scaling$t_sd
  g_sp <- as.integer(species$woodiness_group)
  # spread design: per growth form, an intercept and the climate terms
  B <- .spread_term_cols(p, t, spread_terms)
  Zs <- do.call(cbind, lapply(seq_len(G), function(gi) {
    ind <- as.numeric(g_sp == gi)
    cbind(ind, B * ind)
  }))
  colnames(Zs) <- unlist(lapply(prep$groups, function(g)
    paste0(g, ":", c("(Intercept)", spread_terms))))

  sp <- match(d$taxon_name, species$taxon_name)
  lg_key <- paste(d$location, d$woodiness_group, sep = "|")
  lg_levels <- unique(lg_key)
  lg <- match(lg_key, lg_levels)
  g_lg <- as.integer(factor(sub(".*\\|", "", lg_levels), prep$groups))
  k_key <- paste(d$taxon_name, d$location, sep = "|")
  k <- match(k_key, unique(k_key))

  dat <- list(y = d$log_value, X = X, Zs = Zs, sp = sp, lg = lg, g_lg = g_lg, k = k,
              Pc = d$log_annual_precipitation_c, Tc = d$average_temperature_c)
  beta0 <- stats::lm.fit(X, dat$y)$coefficients
  par <- list(beta = unname(beta0),
              # intercepts of log sd_u start at log(0.5), the climate terms at 0
              theta_u = ifelse(grepl("\\(Intercept\\)$", colnames(Zs)), log(0.5), 0),
              log_sd_P = log(0.2), log_sd_T = log(0.02), log_sd_L = rep(log(0.2), G),
              log_sd_w = log(0.2), log_sigma = log(0.3),
              u = numeric(nrow(species)), sP = numeric(nrow(species)), sT = numeric(nrow(species)),
              L = numeric(length(lg_levels)), w = numeric(max(k)))

  nll <- function(par) {
    RTMB::getAll(par)
    sd_u <- exp(dat$Zs %*% theta_u)
    eta <- dat$X %*% beta + u[dat$sp] + sP[dat$sp] * dat$Pc + sT[dat$sp] * dat$Tc +
      L[dat$lg] + w[dat$k]
    -sum(RTMB::dnorm(u, 0, sd_u, log = TRUE)) -
      sum(RTMB::dnorm(sP, 0, exp(log_sd_P), log = TRUE)) -
      sum(RTMB::dnorm(sT, 0, exp(log_sd_T), log = TRUE)) -
      sum(RTMB::dnorm(L, 0, exp(log_sd_L)[dat$g_lg], log = TRUE)) -
      sum(RTMB::dnorm(w, 0, exp(log_sd_w), log = TRUE)) -
      sum(RTMB::dnorm(dat$y, eta, exp(log_sigma), log = TRUE))
  }
  message("== ", trait, " (RTMB spread model: ", paste(prep$groups, collapse = ", "),
          "; spread terms: ", paste(spread_terms, collapse = ", "), ") ==")
  t0 <- Sys.time()
  obj <- RTMB::MakeADFun(nll, par, random = c("u", "sP", "sT", "L", "w"), silent = TRUE)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(eval.max = 2000, iter.max = 1000))
  # covariance of the fixed parameters from the Hessian of the Laplace
  # objective (sdreport would also compute every random effect's SD)
  hess <- stats::optimHess(opt$par, obj$fn, obj$gr)
  cov_fixed <- tryCatch(solve(hess), error = function(e) NULL)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  pl <- obj$env$parList(opt$par)

  fixed_names <- names(opt$par)
  theta_idx <- which(fixed_names == "theta_u")
  beta_hat <- pl$beta
  names(beta_hat) <- colnames(X)
  sd_u <- exp(as.vector(Zs %*% pl$theta_u))
  gname <- function(g) if (G > 1) paste0("woodiness_group", g) else "(Intercept)"
  vc <- dplyr::bind_rows(
    tibble::tibble(grouping = "taxon_name",
                   term = c("log_annual_precipitation_c", "average_temperature_c"),
                   sd = exp(c(pl$log_sd_P, pl$log_sd_T))),
    tibble::tibble(grouping = "location", term = vapply(prep$groups, gname, ""),
                   sd = exp(pl$log_sd_L)),
    tibble::tibble(grouping = "taxon_name:location", term = "(Intercept)", sd = exp(pl$log_sd_w)),
    tibble::tibble(grouping = "residual", term = "sigma", sd = exp(pl$log_sigma))) |>
    dplyr::mutate(variance = sd^2)
  loc_eff <- tibble::tibble(location = sub("\\|[^|]*$", "", lg_levels),
                            woodiness_group = prep$groups[g_lg],
                            location_effect = pl$L)
  theta <- pl$theta_u
  names(theta) <- colnames(Zs)
  cov_theta <- if (is.null(cov_fixed)) NULL else cov_fixed[theta_idx, theta_idx, drop = FALSE]
  if (!is.null(cov_theta)) dimnames(cov_theta) <- list(colnames(Zs), colnames(Zs))
  structure(list(
    trait = trait, groups = prep$groups, n_species = prep$n_species,
    data = d, fixed_pred = as.vector(X %*% beta_hat), fixed = beta_hat, vc = vc,
    species_variance = tibble::tibble(taxon_name = species$taxon_name, species_variance = sd_u^2),
    location_effects = loc_eff,
    spread = list(terms = spread_terms, theta = theta, cov = cov_theta),
    centre_scaling = scaling,
    logLik = -opt$objective, df = length(opt$par), aic = 2 * opt$objective + 2 * length(opt$par),
    converged = opt$convergence == 0 && !is.null(cov_fixed) &&
      all(is.finite(diag(cov_fixed))) && all(diag(cov_fixed) > 0),
    message = opt$message, seconds = secs),
    class = "rtmb_spread_fit")
}

#' The fitted species spread along one climate axis, the other held at the
#' growth form's median species centre, with a 90% interval (delta method on
#' log sd_u).
#'
#' @param model an `rtmb_spread_fit`.
#' @param axis "precipitation" or "temperature".
#' @return tibble(woodiness_group, axis, climate (mm or °C), species_sd,
#'   lower, upper).
spread_curve <- function(model, axis = c("precipitation", "temperature"), n = 60,
                         level = 0.90) {
  axis <- match.arg(axis)
  s <- model$centre_scaling
  z <- stats::qnorm(1 - (1 - level) / 2)
  species <- dplyr::distinct(model$data, taxon_name, woodiness_group,
                             log_annual_precipitation_mn, average_temperature_mn)
  purrr::map(model$groups, function(g) {
    sg <- dplyr::filter(species, woodiness_group == g)
    if (axis == "precipitation") {
      lp <- seq(stats::quantile(sg$log_annual_precipitation_mn, 0.02),
                stats::quantile(sg$log_annual_precipitation_mn, 0.98), length.out = n)
      tt <- rep(stats::median(sg$average_temperature_mn), n)
      climate <- exp(lp)
    } else {
      tt <- seq(stats::quantile(sg$average_temperature_mn, 0.02),
                stats::quantile(sg$average_temperature_mn, 0.98), length.out = n)
      lp <- rep(stats::median(sg$log_annual_precipitation_mn), n)
      climate <- tt
    }
    B <- .spread_term_cols((lp - s$p_mean) / s$p_sd, (tt - s$t_mean) / s$t_sd, model$spread$terms)
    cols <- paste0(g, ":", c("(Intercept)", model$spread$terms))
    Zg <- cbind(1, B)
    eta <- as.vector(Zg %*% model$spread$theta[cols])
    se <- if (is.null(model$spread$cov)) NA_real_ else
      sqrt(rowSums((Zg %*% model$spread$cov[cols, cols]) * Zg))
    tibble::tibble(woodiness_group = g, axis = axis, climate = climate,
                   species_sd = exp(eta), lower = exp(eta - z * se), upper = exp(eta + z * se))
  }) |> purrr::list_rbind()
}
