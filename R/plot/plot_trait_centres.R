# Across-species view of a fitted trait: predicted trait value at each species'
# climate centre (medoid) against the centre climate. Shows the between-species
# turnover the model estimates.

#' @param centres output of `species_centres()` for one trait.
#' @param trait trait label (axis/title).
#' @param x which centre climate on the x-axis: "precip" or "temperature".
#' @param log_y log-scale the trait axis (TRUE for positive trait values).
plot_trait_centres <- function(centres, trait, x = c("precip", "temperature"),
                               log_y = TRUE) {
  x <- match.arg(x)
  xv <- if (x == "precip") "annual_precipitation_centre" else "average_temperature_centre"
  xlab <- if (x == "precip") "Annual precipitation at centre (mm)" else "Mean annual temperature at centre (°C)"

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
