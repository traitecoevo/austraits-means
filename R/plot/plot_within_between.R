# Within- vs between-species climate responses of one trait's model, in the
# style of the earlier height report (Data_height_AVH.qmd, retired; in git
# history): the records (grey),
# each species' fitted line across the climate range of its records (black),
# the species' modelled means at their climate centre with the
# between-species line (red), and the average within-species line (blue).
# One panel per climate variable; the other variable is held at each
# species' centre (for the red and blue lines, at the mean centre).

.wb_axes <- list(
  precip = list(var = "log_annual_precipitation", other = "average_temperature",
                label = "Annual precipitation (mm)", log_x = TRUE, per = "per log mm"),
  temperature = list(var = "average_temperature", other = "log_annual_precipitation",
                     label = "Mean annual temperature (°C)", log_x = FALSE, per = "per °C")
)

#' Pieces of a fit for one growth form: its fixed effects with the pooled
#' model's term names, its records, and the species' random effects.
#'
#' A pooled `trait_climate_fit` has one set for all species (`group` NULL);
#' a location-model `rtmb_spread_fit` has one per growth form.
.wb_parts <- function(fit, group = NULL) {
  if (inherits(fit, "rtmb_spread_fit")) {
    nm <- .group_coef_names(fit, group)
    fixed <- fit$fixed_table[match(nm, fit$fixed_table$term), ]
    fixed$term <- c("(Intercept)", woodiness_model_terms)
    list(fixed = fixed,
         data = dplyr::filter(fit$data, woodiness_group == group),
         ranef = dplyr::select(dplyr::filter(fit$centres, woodiness_group == group),
                               taxon_name, dplyr::starts_with("g_")))
  } else {
    list(fixed = fit$fixed, data = fit$data, ranef = fit$ranef)
  }
}

#' The lines and points of the figure, on the model's (log) scale.
#'
#' @param fit a pooled `trait_climate_fit` (with `$data`, `$fixed`,
#'   `$ranef`) or a location-model `rtmb_spread_fit` with centres.
#' @param axis "precip" or "temperature".
#' @param group for a location-model fit, the growth form.
#' @return list(raw, species, means, between, within), each tibble(x, y)
#'   (species lines and means also `taxon_name`).
within_between_lines <- function(fit, axis = c("precip", "temperature"), group = NULL) {
  a <- .wb_axes[[match.arg(axis)]]
  parts <- .wb_parts(fit, group)
  b <- stats::setNames(parts$fixed$estimate, parts$fixed$term)
  v <- a$var
  mn <- paste0(v, "_mn")
  omn <- paste0(a$other, "_mn")
  slope_c <- b[[paste0(v, "_c")]]

  sp <- parts$data |>
    dplyr::group_by(taxon_name) |>
    dplyr::summarise(lo = min(.data[[v]]), hi = max(.data[[v]]),
                     m = dplyr::first(.data[[mn]]), om = dplyr::first(.data[[omn]])) |>
    dplyr::inner_join(parts$ranef, by = "taxon_name")
  # each species' expected log value at its centre (its modelled mean)
  centre <- b[["(Intercept)"]] + b[[mn]] * sp$m + b[[omn]] * sp$om + sp$g_0
  slope_sp <- slope_c + sp[[paste0("g_", v, "_c")]]

  mbar <- mean(sp$m)
  at_mean <- b[["(Intercept)"]] + b[[mn]] * mbar + b[[omn]] * mean(sp$om)
  xr <- range(sp$m)
  list(
    raw = tibble::tibble(x = parts$data[[v]], y = parts$data$log_value),
    species = dplyr::bind_rows(
      tibble::tibble(taxon_name = sp$taxon_name, x = sp$lo, y = centre + slope_sp * (sp$lo - sp$m)),
      tibble::tibble(taxon_name = sp$taxon_name, x = sp$hi, y = centre + slope_sp * (sp$hi - sp$m))),
    means = tibble::tibble(taxon_name = sp$taxon_name, x = sp$m, y = centre),
    between = tibble::tibble(x = xr, y = at_mean + b[[mn]] * (xr - mbar)),
    within = tibble::tibble(x = xr, y = at_mean + slope_c * (xr - mbar))
  )
}

wb_colours <- c("Raw data" = "grey", "Species means (turnover)" = "red",
                "Average within species" = "blue", "Within species" = "black")

#' One panel: the trait against one climate variable (for a location-model
#' fit, one growth form).
plot_within_between_panel <- function(fit, axis = c("precip", "temperature"), unit = NULL,
                                      group = NULL) {
  axis <- match.arg(axis)
  a <- .wb_axes[[axis]]
  L <- within_between_lines(fit, axis, group)
  tx <- if (a$log_x) exp else identity
  L <- lapply(L, \(d) dplyr::mutate(d, x = tx(x), y = exp(y)))
  parts <- .wb_parts(fit, group)
  fx <- dplyr::filter(parts$fixed, term %in% paste0(a$var, c("_mn", "_c")))
  form <- if (is.null(group)) "" else
    paste0(c(woody = "Woody", non_woody = "Non-woody")[[group]], " (",
           dplyr::n_distinct(parts$data$taxon_name), " species)\n")
  slope <- function(t) {
    r <- fx[fx$term == t, ]
    sprintf("%.3f (%.3f to %.3f)", r$estimate, r$conf.low, r$conf.high)
  }
  unit_label <- if (is.null(unit) || is.na(unit)) "" else paste0(" (", unit, ")")
  raw_alpha <- min(0.5, max(0.2, 3000 / nrow(L$raw)))

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x, y)) +
    ggplot2::geom_point(data = L$raw, ggplot2::aes(colour = "Raw data"),
                        alpha = raw_alpha, size = 0.7, stroke = 0) +
    ggplot2::geom_line(data = L$species, ggplot2::aes(group = taxon_name, colour = "Within species"),
                       alpha = 0.5, linewidth = 0.25) +
    ggplot2::geom_point(data = L$means, ggplot2::aes(colour = "Species means (turnover)"),
                        size = 1, alpha = 0.3, stroke = 0) +
    ggplot2::geom_line(data = L$between, ggplot2::aes(colour = "Species means (turnover)"),
                       linewidth = 0.9) +
    ggplot2::geom_line(data = L$within, ggplot2::aes(colour = "Average within species"),
                       linewidth = 0.9) +
    ggplot2::scale_colour_manual(
      values = wb_colours, limits = names(wb_colours),
      guide = ggplot2::guide_legend(override.aes = list(
        shape = c(16, 16, NA, NA), linetype = c(0, 1, 1, 1), alpha = 1, size = 1.5))) +
    ggplot2::scale_y_log10(labels = scales::label_number(drop0trailing = TRUE, big.mark = ",")) +
    ggplot2::labs(x = a$label, y = paste0(fit$trait, unit_label), colour = NULL,
                  subtitle = paste0(form, "Between species ", slope(paste0(a$var, "_mn")),
                                    ";\nwithin species ", slope(paste0(a$var, "_c")),
                                    " ", a$per, " (log trait)")) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8.5, colour = "#3d3d3a"),
                   legend.key.width = grid::unit(1.2, "cm"))
  if (a$log_x) p <- p + ggplot2::scale_x_log10(labels = scales::label_number(big.mark = ","))
  p
}

#' Both climate variables for one trait; for a location-model fit, one row
#' per growth form, each with its own coefficients.
plot_within_between <- function(fit, unit = NULL) {
  n_sp <- dplyr::n_distinct(fit$data$taxon_name)
  location <- inherits(fit, "rtmb_spread_fit")
  groups <- if (location) fit$groups else list(NULL)
  panels <- unlist(lapply(groups, function(g) list(
    plot_within_between_panel(fit, "precip", unit, g),
    plot_within_between_panel(fit, "temperature", unit, g))), recursive = FALSE)
  patchwork::wrap_plots(panels, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(fit$trait, ": within- and between-species responses to climate"),
      caption = paste0(scales::label_comma()(nrow(fit$data)), " records of ",
                       scales::label_comma()(n_sp), " species; ",
                       if (location) paste("location model (RTMB): species lines and means at an",
                                           "average location; ") else "pooled model (TMB); ",
                       "95% intervals; ",
                       "each panel holds the other climate variable at the species' centre"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                             plot.caption = ggplot2::element_text(colour = "#6b6a66")))
}
