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

#' The lines and points of the figure, on the model's (log) scale.
#'
#' @param fit a TMB `trait_climate_fit` (with `$data`, `$fixed`, `$ranef`).
#' @param axis "precip" or "temperature".
#' @return list(raw, species, means, between, within), each tibble(x, y)
#'   (species lines and means also `taxon_name`).
within_between_lines <- function(fit, axis = c("precip", "temperature")) {
  a <- .wb_axes[[match.arg(axis)]]
  b <- stats::setNames(fit$fixed$estimate, fit$fixed$term)
  v <- a$var
  mn <- paste0(v, "_mn")
  omn <- paste0(a$other, "_mn")
  slope_c <- b[[paste0(v, "_c")]]

  sp <- fit$data |>
    dplyr::group_by(taxon_name) |>
    dplyr::summarise(lo = min(.data[[v]]), hi = max(.data[[v]]),
                     m = dplyr::first(.data[[mn]]), om = dplyr::first(.data[[omn]])) |>
    dplyr::inner_join(fit$ranef, by = "taxon_name")
  # each species' expected log value at its centre (its modelled mean)
  centre <- b[["(Intercept)"]] + b[[mn]] * sp$m + b[[omn]] * sp$om + sp$g_0
  slope_sp <- slope_c + sp[[paste0("g_", v, "_c")]]

  mbar <- mean(sp$m)
  at_mean <- b[["(Intercept)"]] + b[[mn]] * mbar + b[[omn]] * mean(sp$om)
  xr <- range(sp$m)
  list(
    raw = tibble::tibble(x = fit$data[[v]], y = fit$data$log_value),
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

#' One panel: the trait against one climate variable.
plot_within_between_panel <- function(fit, axis = c("precip", "temperature"), unit = NULL) {
  axis <- match.arg(axis)
  a <- .wb_axes[[axis]]
  L <- within_between_lines(fit, axis)
  tx <- if (a$log_x) exp else identity
  L <- lapply(L, \(d) dplyr::mutate(d, x = tx(x), y = exp(y)))
  fx <- dplyr::filter(fit$fixed, term %in% paste0(a$var, c("_mn", "_c")))
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
                  subtitle = paste0("Between species ", slope(paste0(a$var, "_mn")),
                                    "; within species ", slope(paste0(a$var, "_c")),
                                    " ", a$per, " (log trait)")) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8.5, colour = "#3d3d3a"),
                   legend.key.width = grid::unit(1.2, "cm"))
  if (a$log_x) p <- p + ggplot2::scale_x_log10(labels = scales::label_number(big.mark = ","))
  p
}

#' Both climate variables for one trait.
plot_within_between <- function(fit, unit = NULL) {
  n_sp <- dplyr::n_distinct(fit$data$taxon_name)
  patchwork::wrap_plots(
    plot_within_between_panel(fit, "precip", unit),
    plot_within_between_panel(fit, "temperature", unit),
    ncol = 2, guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = paste0(fit$trait, ": within- and between-species responses to climate"),
      caption = paste0(scales::label_comma()(nrow(fit$data)), " records of ",
                       scales::label_comma()(n_sp), " species; TMB fit; 95% intervals; ",
                       "each panel holds the other climate variable at the species' centre"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 12),
                             plot.caption = ggplot2::element_text(colour = "#6b6a66")))
}
