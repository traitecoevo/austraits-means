# Location model with the between-species spread a continuous function of
# climate (#41), in RTMB: glmmTMB cannot make a random-effect variance depend
# on a covariate.
#
# Same model as fit_woodiness_model(location = ...) except the species
# effect: for species s of growth form g with climate centre (P_s, T_s),
# standardised to (p, t),
#   u_s ~ N(0, sd_u(s)^2),
#   log sd_u(s) = a_g + b1_g p + b2_g t.
# Linear by default (decided 25 November 2026): the squares p^2 and t^2 can
# be added through `spread_terms`, but a quadratic in log SD grows or
# collapses fast outside the species' centres, the product predicts from
# climate anywhere, and the quadratic was only supported for leaf area and
# height, on top of a linear term that already was.
# The rest as before: growth-form intercepts and climate slopes (fixed),
# uncorrelated species random slopes on the centred climate, a location
# effect per growth form, a species x location effect and the residual.
# Random effects are integrated out by the Laplace approximation, as in
# glmmTMB, so for this Gaussian model the likelihood is exact and AICs are
# comparable with the glmmTMB fits on the same records.
#
# With `spread_terms = character(0)`, log sd_u(s) = a_g: the glmmTMB location
# model, which checks this code (the likelihoods must agree).
#
# Shrinkage (fit_spread_shrunk()): without it the free model overfits small
# groups (the free quadratic gave non-woody seed mass, 27 species, a hump to
# SD 2.5) and can send the SD to ~0 at the edge of the data. The climate
# coefficients b get a prior b ~ N(0, tau_d^2), one tau per polynomial degree
# d (linear p, t; second-order p^2, t^2, p t when used) shared by the growth
# forms, the taus estimated from the data: a ridge penalty with its strength
# chosen by marginal likelihood, as mgcv chooses smoothing parameters. A
# growth form whose species say little about how spread changes with climate
# is pulled toward constant spread, and tau_d -> 0 drops degree d.
#
# b is NOT integrated inside the same Laplace approximation as the species
# effects u. That was tried and fails on real data (NaN objective for LMA):
# a species seen in one cell can put its effect into the species x location
# term, leaving u_s ~ 0, and the joint density then grows without bound as
# sd_u(s) -> 0, so with a loose prior the joint mode in (u, b) runs into
# that funnel. Instead, empirical Bayes in two steps, integrating u first:
#   1. from the free fit's b-hat and covariance V (the Laplace marginal
#      likelihood is close to quadratic in b), tau maximises
#      N(b-hat; 0, V + T), T = diag(tau^2) (estimate_spread_shrinkage());
#   2. refit with b fixed-effect and the penalty b^2 / (2 tau^2): the
#      posterior mode given tau, its covariance from the penalised Hessian.

.spread_term_cols <- function(p, t, terms) {
  all <- cbind(p = p, t = t, p2 = p^2, t2 = t^2, pt = p * t)
  all[, terms, drop = FALSE]
}

# Polynomial degree of each spread term: the squares and the interaction are
# second order and share one tau.
.spread_degree <- function(term) ifelse(term %in% c("p2", "t2", "pt"), 2L, 1L)

#' Fit the location model with the species spread a smooth function of the
#' species' climate centre.
#'
#' @param data records with a `location` column (`add_location_key()`).
#' @param spread_terms terms of log sd_u in the standardised centres: any of
#'   "p", "t", "p2", "t2", "pt" (precipitation, temperature, their squares,
#'   their product).
#'   Linear by default; see the file header.
#' @param penalty_tau optional prior SD of each climate coefficient, one per
#'   element of `spread_terms` (the same for every growth form): adds the
#'   ridge penalty b^2 / (2 tau^2). Use fit_spread_shrunk() to choose it.
#' @param centres also return each species' modelled mean at its climate
#'   centre, with its `level` interval (rtmb_species_centres()).
#' @return an `rtmb_spread_fit`: list(trait, groups, data, fixed_pred, vc,
#'   species_variance, location_effects, spread (coefficients and their
#'   covariance), centre_scaling, logLik, df, aic, converged, seconds, and
#'   with `centres = TRUE` the species centres). With
#'   a penalty, the coefficients are the posterior mode, the covariance the
#'   posterior one, logLik is unpenalised at that mode and aic is NA (see
#'   fit_spread_shrunk()).
fit_spread_rtmb <- function(data, trait, groups, n_min = 5, min_group_species = 20,
                            spread_terms = c("p", "t"), penalty_tau = NULL,
                            centres = FALSE, level = 0.94) {
  penalised <- !is.null(penalty_tau)
  if (penalised) stopifnot(length(penalty_tau) == length(spread_terms), all(penalty_tau > 0))
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
  # spread design: per growth form an intercept (Za, coefficients a) and the
  # climate terms (Zb, coefficients b)
  B <- .spread_term_cols(p, t, spread_terms)
  K <- length(spread_terms)
  Za <- vapply(seq_len(G), function(gi) as.numeric(g_sp == gi), numeric(nrow(species)))
  Za <- matrix(Za, ncol = G)
  Zb <- do.call(cbind, lapply(seq_len(G), function(gi) B * as.numeric(g_sp == gi)))
  a_names <- paste0(prep$groups, ":(Intercept)")
  # (paste0 with a zero-length argument still returns the prefix)
  b_names <- if (K > 0) unlist(lapply(prep$groups, function(g) paste0(g, ":", spread_terms)))
  # the order of the coefficients in the result, as before: per growth form,
  # the intercept then its climate terms
  theta_names <- unlist(lapply(prep$groups, function(g)
    paste0(g, ":", c("(Intercept)", spread_terms))))
  # prior SD of each b, in the order of Zb's columns
  b_tau <- if (penalised) rep(unname(penalty_tau), G) else numeric(0)

  sp <- match(d$taxon_name, species$taxon_name)
  lg_key <- paste(d$location, d$woodiness_group, sep = "|")
  lg_levels <- unique(lg_key)
  lg <- match(lg_key, lg_levels)
  g_lg <- as.integer(factor(sub(".*\\|", "", lg_levels), prep$groups))
  k_key <- paste(d$taxon_name, d$location, sep = "|")
  k <- match(k_key, unique(k_key))

  dat <- list(y = d$log_value, X = X, Za = Za, Zb = Zb, b_tau = b_tau,
              sp = sp, lg = lg, g_lg = g_lg, k = k,
              Pc = d$log_annual_precipitation_c, Tc = d$average_temperature_c)
  beta0 <- stats::lm.fit(X, dat$y)$coefficients
  # the inner (Laplace) optimisation starts from the species' mean residuals
  # rather than 0, where the species density has no curvature in sd_u
  u0 <- as.vector(tapply(dat$y - as.vector(X %*% beta0), factor(sp, seq_len(nrow(species))), mean))
  par <- list(beta = unname(beta0),
              # intercepts of log sd_u start at log(0.5), the climate terms at 0
              a = rep(log(0.5), G),
              log_sd_P = log(0.2), log_sd_T = log(0.02), log_sd_L = rep(log(0.2), G),
              log_sd_w = log(0.2), log_sigma = log(0.3),
              u = u0, sP = numeric(nrow(species)), sT = numeric(nrow(species)),
              L = numeric(length(lg_levels)), w = numeric(max(k)))
  if (K > 0) par$b <- numeric(G * K)

  nll <- function(par) {
    RTMB::getAll(par)
    log_sd_u <- dat$Za %*% a
    if (K > 0) log_sd_u <- log_sd_u + dat$Zb %*% b
    eta <- dat$X %*% beta + u[dat$sp] + sP[dat$sp] * dat$Pc + sT[dat$sp] * dat$Tc +
      L[dat$lg] + w[dat$k]
    penalty <- if (penalised) sum(b^2 / (2 * dat$b_tau^2)) else 0
    penalty - sum(RTMB::dnorm(u, 0, exp(log_sd_u), log = TRUE)) -
      sum(RTMB::dnorm(sP, 0, exp(log_sd_P), log = TRUE)) -
      sum(RTMB::dnorm(sT, 0, exp(log_sd_T), log = TRUE)) -
      sum(RTMB::dnorm(L, 0, exp(log_sd_L)[dat$g_lg], log = TRUE)) -
      sum(RTMB::dnorm(w, 0, exp(log_sd_w), log = TRUE)) -
      sum(RTMB::dnorm(dat$y, eta, exp(log_sigma), log = TRUE))
  }
  message("== ", trait, " (RTMB spread model: ", paste(prep$groups, collapse = ", "),
          "; spread terms: ", paste(spread_terms, collapse = ", "),
          if (penalised) paste0("; penalty tau ", paste(signif(penalty_tau, 2), collapse = ", ")),
          ") ==")
  t0 <- Sys.time()
  obj <- RTMB::MakeADFun(nll, par, random = c("u", "sP", "sT", "L", "w"), silent = TRUE)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(eval.max = 2000, iter.max = 1000))
  pl <- obj$env$parList(opt$par)
  theta_hat <- c(stats::setNames(pl$a, a_names),
                 if (K > 0) stats::setNames(pl$b, b_names))[theta_names]
  # covariance of the fixed parameters from the Hessian of the Laplace
  # objective (sdreport would also compute every random effect's SD)
  hess <- stats::optimHess(opt$par, obj$fn, obj$gr)
  cov_fixed <- tryCatch(solve(hess), error = function(e) NULL)
  pd <- !is.null(cov_fixed) && all(is.finite(diag(cov_fixed))) && all(diag(cov_fixed) > 0)
  cov_theta <- NULL
  if (!is.null(cov_fixed)) {
    fixed_names <- names(opt$par)
    idx <- c(which(fixed_names == "a"), which(fixed_names == "b"))
    cov_theta <- cov_fixed[idx, idx, drop = FALSE]
    dimnames(cov_theta) <- list(c(a_names, b_names), c(a_names, b_names))
    cov_theta <- cov_theta[theta_names, theta_names, drop = FALSE]
  }
  sp_centres <- if (centres) {
    # the fixed-effect design of each species at its climate centre: its
    # growth form, its medoid, within-species deviations 0
    nd <- dplyr::mutate(species, log_annual_precipitation_c = 0, average_temperature_c = 0)
    X_sp <- stats::model.matrix(stats::as.formula(paste("~", prep$fixed)), nd)
    stopifnot(identical(colnames(X_sp), colnames(X)))
    .rtmb_species_centres(obj, opt, hess, X_sp, level) |>
      dplyr::mutate(taxon_name = species$taxon_name, woodiness_group = species$woodiness_group,
                    log_annual_precipitation_mn = species$log_annual_precipitation_mn,
                    average_temperature_mn = species$average_temperature_mn,
                    g_0 = pl$u, g_log_annual_precipitation_c = pl$sP,
                    g_average_temperature_c = pl$sT, .before = 1)
  }
  secs <- as.numeric(Sys.time() - t0, units = "secs")

  beta_hat <- pl$beta
  names(beta_hat) <- colnames(X)
  # fixed effects with 95% intervals, as glmmTMB's tidy table
  beta_se <- if (is.null(cov_fixed)) rep(NA_real_, length(beta_hat)) else
    sqrt(diag(cov_fixed)[names(opt$par) == "beta"])
  fixed_table <- tibble::tibble(term = names(beta_hat), estimate = unname(beta_hat),
                                std.error = beta_se,
                                conf.low = estimate - stats::qnorm(0.975) * std.error,
                                conf.high = estimate + stats::qnorm(0.975) * std.error)
  sd_u <- exp(as.vector(Za %*% pl$a + if (K > 0) Zb %*% pl$b else 0))
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
  structure(list(
    trait = trait, groups = prep$groups, n_species = prep$n_species,
    data = d, fixed_pred = as.vector(X %*% beta_hat), fixed = beta_hat,
    fixed_table = fixed_table, vc = vc,
    species_variance = tibble::tibble(taxon_name = species$taxon_name, species_variance = sd_u^2),
    location_effects = loc_eff,
    spread = list(terms = spread_terms, theta = theta_hat, cov = cov_theta,
                  penalty_tau = penalty_tau),
    centre_scaling = scaling, centres = sp_centres,
    logLik =-(opt$objective - if (penalised) sum(pl$b^2 / (2 * b_tau^2)) else 0),
    df = length(opt$par),
    aic = if (penalised) NA_real_ else 2 * opt$objective + 2 * length(opt$par),
    converged = opt$convergence == 0 && pd,
    message = opt$message, seconds = secs),
    class = "rtmb_spread_fit")
}

#' Each species' modelled mean at its climate centre, with a delta-method
#' interval, as glmmTMB's predict(se.fit = TRUE) gives for the pooled model.
#'
#' The mean is x_s beta + u_s: the fixed part at the species' medoid (`X_sp`,
#' one row per species, within-species deviations 0) plus its species effect;
#' location and species x location effects are 0, an average location. Its
#' variance, as in TMB's sdreport(), is the conditional variance of u_s given
#' the fixed parameters theta plus the propagated uncertainty of theta:
#'   var = (H^-1)_ss + g_s V g_s',  g_s = d(x_s beta + u-hat_s) / d theta,
#' H the Hessian of the joint likelihood in the random effects, V the
#' covariance of theta, and d u-hat / d theta = -H^-1 H_u,theta.
#' sdreport() would do this for an ADREPORTed vector with a dense solve of
#' (random effects x species), about 140k x 10k for height; here only the
#' diagonal of H^-1 (sparse selected inverse, inside sdreport()) and one
#' solve per fixed parameter are needed.
#'
#' @return tibble(log_estimate, log_se, log_conf.low, log_conf.high, estimate,
#'   conf.low, conf.high), one row per species in the fit's species order.
.rtmb_species_centres <- function(obj, opt, hess, X_sp, level) {
  sdr <- RTMB::sdreport(obj, par.fixed = opt$par, hessian.fixed = hess,
                        getJointPrecision = TRUE)
  Q <- sdr$jointPrecision                  # full parameter order
  r <- obj$env$random
  nonr <- setdiff(seq_len(nrow(Q)), r)
  stopifnot(identical(rownames(Q)[nonr], names(opt$par)))
  V <- sdr$cov.fixed
  # d u-hat / d theta for every random effect (n_random x n_fixed)
  J <- -as.matrix(Matrix::solve(Q[r, r], as.matrix(Q[r, nonr])))
  u_i <- which(rownames(Q)[r] == "u")
  stopifnot(length(u_i) == nrow(X_sp))
  Ju <- J[u_i, , drop = FALSE]
  # conditional variance of u given theta: sdreport's total less the theta term
  cond_var <- sdr$diag.cov.random[u_i] - rowSums((Ju %*% V) * Ju)
  g <- Ju
  b_i <- which(names(opt$par) == "beta")
  g[, b_i] <- g[, b_i] + X_sp
  se <- sqrt(cond_var + rowSums((g %*% V) * g))
  pl <- obj$env$parList(opt$par)
  est <- as.vector(X_sp %*% pl$beta) + pl$u
  z <- stats::qnorm(1 - (1 - level) / 2)
  tibble::tibble(log_estimate = est, log_se = se,
                 log_conf.low = est - z * se, log_conf.high = est + z * se,
                 estimate = exp(est), conf.low = exp(est - z * se),
                 conf.high = exp(est + z * se))
}

#' Empirical-Bayes shrinkage of the spread climate terms of a free fit.
#'
#' Treats the free fit's climate coefficients b-hat as b-hat ~ N(b, V), with
#' V their block of the fixed-parameter covariance, and b ~ N(0, T),
#' T = diag(tau_d^2) with one tau per polynomial degree d (shared by the
#' growth forms). tau maximises the marginal density N(b-hat; 0, V + T).
#' Each tau is estimated from only 2 x (terms of that degree) coefficients,
#' so it is noisy: a chance 3-SE coefficient keeps its degree partly in.
#' The same approximation gives log_bf, the log Bayes factor of the shrunk
#' model against constant spread (T -> 0): log N(b-hat; 0, V + T) -
#' log N(b-hat; 0, V).
#'
#' @param theta,cov a free fit's `spread$theta` and `spread$cov`.
#' @param tau_bounds search interval for each tau; a tau at the lower bound
#'   means that degree is dropped.
#' @return list(tau (per degree), tau_term (per spread term), at_lower (per
#'   degree), b_post, cov_post (the Gaussian-approximate posterior, for
#'   checking the refit), log_bf).
estimate_spread_shrinkage <- function(theta, cov, tau_bounds = c(1e-3, 10)) {
  b_names <- names(theta)[!grepl("\\(Intercept\\)$", names(theta))]
  b <- theta[b_names]
  V <- cov[b_names, b_names, drop = FALSE]
  term <- sub(".*:", "", b_names)
  degree <- .spread_degree(term)
  degrees <- sort(unique(degree))
  k <- length(b)
  log_dmvnorm0 <- function(x, S) {
    R <- chol(S)
    -0.5 * (k * log(2 * pi) + 2 * sum(log(diag(R))) + sum(backsolve(R, x, transpose = TRUE)^2))
  }
  Tmat <- function(lt) diag(exp(2 * lt[match(degree, degrees)]), k)
  fit <- stats::optim(rep(log(0.2), length(degrees)),
                      function(lt) -log_dmvnorm0(b, V + Tmat(lt)),
                      method = "L-BFGS-B", lower = log(tau_bounds[1]), upper = log(tau_bounds[2]))
  Tm <- Tmat(fit$par)
  cov_post <- solve(solve(V) + solve(Tm))
  b_post <- as.vector(cov_post %*% solve(V, b))
  names(b_post) <- b_names
  tau <- stats::setNames(exp(fit$par), paste0("degree ", degrees))
  term_levels <- unique(term)
  tau_term <- stats::setNames(tau[match(.spread_degree(term_levels), degrees)],
                              term_levels)
  list(tau = tau, tau_term = tau_term, b_post = b_post, cov_post = cov_post,
       at_lower = tau <= tau_bounds[1] * 1.01,
       log_bf = -fit$value - log_dmvnorm0(b, V))
}

#' Fit the spread model with its climate terms shrunk (see the file header):
#' choose tau from a free fit, drop the degrees whose tau is at the lower
#' bound, refit with the penalty.
#'
#' @param free a free `rtmb_spread_fit` from the same data, or the summary of
#'   one that run_spread_rtmb.R stores (with `theta` and `cov`).
#' @return the penalised `rtmb_spread_fit`, with `spread$shrinkage` (the
#'   estimate_spread_shrinkage() result, plus `edf`) and the AIC of a
#'   penalised fit: -2 logLik (unpenalised, at the mode) + 2 x (parameters
#'   other than b + effective df of b), edf = k - tr(C_bb P), C the
#'   posterior covariance and P = diag(1 / tau^2) the penalty. The taus are
#'   not counted, as in mgcv's AIC; on the scale of the free and constant
#'   fits' AICs.
fit_spread_shrunk <- function(data, trait, groups, free, n_min = 5, min_group_species = 20,
                              centres = FALSE, level = 0.94) {
  # a fit, or the summary run_spread_rtmb.R stores (theta and cov at the top)
  fs <- if (!is.null(free$spread)) free$spread else free
  stopifnot(!is.null(fs$theta), !is.null(fs$cov))
  sh <- estimate_spread_shrinkage(fs$theta, fs$cov)
  keep <- names(sh$tau_term)[!sh$at_lower[paste0("degree ", .spread_degree(names(sh$tau_term)))]]
  m <- fit_spread_rtmb(data, trait, groups, n_min, min_group_species,
                       spread_terms = keep,
                       penalty_tau = if (length(keep)) unname(sh$tau_term[keep]),
                       centres = centres, level = level)
  b_names <- names(m$spread$theta)[!grepl("\\(Intercept\\)$", names(m$spread$theta))]
  k <- length(b_names)
  sh$edf <- if (k == 0) 0 else
    k - sum(diag(m$spread$cov[b_names, b_names, drop = FALSE]) /
              sh$tau_term[sub(".*:", "", b_names)]^2)
  m$spread$shrinkage <- sh
  m$aic <- -2 * m$logLik + 2 * (m$df - k + sh$edf)
  m
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
