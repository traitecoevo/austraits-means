# Across-species view of a fitted trait: predicted trait value at each species'
# climate centre (medoid) against the centre climate. Shows the between-species
# turnover the model estimates.

#' @param centres output of `species_centres()` for one trait.
#' @param trait trait label (axis/title).
#' @param x which centre climate on the x-axis: "precip" or "temperature".
#' @param log_y log-scale the trait axis (TRUE for positive trait values).
#' @param fit optional location-model fit: draws, per growth form, the
#'   model's distribution of species means against the climate centre (the
#'   between-species line, and ribbons holding the central 50% and 90% of
#'   species, from its climate-dependent species SD; the other climate
#'   variable at the form's median species centre) in place of a fitted trend.
plot_trait_centres <- function(centres, trait, x = c("precip", "temperature"),
                               log_y = TRUE, fit = NULL) {
  x <- match.arg(x)
  xv <- if (x == "precip") "annual_precipitation_centre" else "average_temperature_centre"
  xlab <- if (x == "precip") "Annual precipitation at centre (mm)" else "Mean annual temperature at centre (°C)"
  form_levels <- c(woody = "woody", non_woody = "non-woody")
  form_colours <- c(woody = "#1b7837", `non-woody` = "#b2182b")

  if (!is.null(fit)) {
    band <- trait_centre_band(fit, if (x == "precip") "precipitation" else "temperature") |>
      dplyr::mutate(woodiness_group = factor(form_levels[woodiness_group], form_levels))
    pts <- centres |>
      dplyr::mutate(woodiness_group = factor(form_levels[as.character(woodiness_group)], form_levels))
    p <- ggplot2::ggplot(mapping = ggplot2::aes(colour = woodiness_group, fill = woodiness_group)) +
      ggplot2::geom_ribbon(data = band, ggplot2::aes(climate, ymin = lower90, ymax = upper90),
                           alpha = 0.12, colour = NA) +
      ggplot2::geom_ribbon(data = band, ggplot2::aes(climate, ymin = lower50, ymax = upper50),
                           alpha = 0.18, colour = NA) +
      ggplot2::geom_point(data = pts, ggplot2::aes(.data[[xv]], estimate),
                          alpha = 0.35, stroke = 0) +
      ggplot2::geom_line(data = band, ggplot2::aes(climate, mean), linewidth = 1) +
      ggplot2::scale_colour_manual(values = form_colours, name = NULL, drop = TRUE) +
      ggplot2::scale_fill_manual(values = form_colours, name = NULL, drop = TRUE) +
      ggplot2::labs(x = xlab, y = paste0(trait, " at centre"), title = trait,
                    caption = paste(
                      "Points: each species' modelled mean (location model). Line and bands: the",
                      "model's mean and central 50% and 90% of species\nmeans at that climate centre,",
                      "the other climate variable at the median species centre. Species with few",
                      "records are pulled toward the line.")) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, colour = "#6b6a66", hjust = 0))
    if (x == "precip") p <- p + ggplot2::scale_x_log10()
    if (log_y) p <- p + ggplot2::scale_y_log10()
    return(p)
  }

  by_form <- "woodiness_group" %in% names(centres)
  p <- if (by_form) {
    # the location model: woody and non-woody species have their own
    # coefficients, so each gets its own points and trend
    centres |>
      dplyr::mutate(woodiness_group = factor(woodiness_group, c("woody", "non_woody"),
                                             c("woody", "non-woody"))) |>
      ggplot2::ggplot(ggplot2::aes(.data[[xv]], estimate, colour = woodiness_group)) +
      ggplot2::geom_point(alpha = 0.35, stroke = 0) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1) +
      ggplot2::scale_colour_manual(values = c(woody = "#1b7837", `non-woody` = "#b2182b"),
                                   name = NULL, drop = TRUE)
  } else {
    centres |>
      ggplot2::ggplot(ggplot2::aes(.data[[xv]], estimate)) +
      ggplot2::geom_point(alpha = 0.35, stroke = 0, colour = "#2C6E8E") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = "#fb0019", se = FALSE)
  }
  p <- p +
    ggplot2::labs(x = xlab, y = paste0(trait, " at centre"), title = trait) +
    ggplot2::theme_classic()

  if (x == "precip") p <- p + ggplot2::scale_x_log10()
  if (log_y) p <- p + ggplot2::scale_y_log10()
  p
}

#' The location model's distribution of species means along one climate
#' axis of the species' centre, per growth form: the between-species line
#' (fixed part at the centre, within-species terms 0) and the central 50% and
#' 90% of species means from the species SD there (spread_curve()), in trait
#' units.
#'
#' @return tibble(woodiness_group, climate, mean, lower50, upper50, lower90,
#'   upper90).
trait_centre_band <- function(fit, axis = c("precipitation", "temperature")) {
  axis <- match.arg(axis)
  sc <- spread_curve(fit, axis)
  med <- dplyr::distinct(fit$data, taxon_name, woodiness_group,
                         log_annual_precipitation_mn, average_temperature_mn) |>
    dplyr::group_by(woodiness_group) |>
    dplyr::summarise(p = stats::median(log_annual_precipitation_mn),
                     t = stats::median(average_temperature_mn)) |>
    dplyr::mutate(woodiness_group = as.character(woodiness_group))
  purrr::map(fit$groups, function(g) {
    b <- fit$fixed[.group_coef_names(fit, g)]
    s <- dplyr::filter(sc, woodiness_group == g)
    m <- med[med$woodiness_group == g, ]
    lp <- if (axis == "precipitation") log(s$climate) else m$p
    tt <- if (axis == "temperature") s$climate else m$t
    eta <- b[[1]] + b[[2]] * lp + b[[3]] * tt
    q <- function(pr) exp(eta + stats::qnorm(pr) * s$species_sd)
    tibble::tibble(woodiness_group = g, climate = s$climate, mean = exp(eta),
                   lower50 = q(0.25), upper50 = q(0.75), lower90 = q(0.05), upper90 = q(0.95))
  }) |> purrr::list_rbind()
}
