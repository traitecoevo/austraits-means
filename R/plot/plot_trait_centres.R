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

  p <- centres |>
    ggplot2::ggplot(ggplot2::aes(.data[[xv]], estimate)) +
    ggplot2::geom_point(alpha = 0.35, stroke = 0, colour = "#2C6E8E") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = "#fb0019", se = FALSE) +
    ggplot2::labs(x = xlab, y = paste0(trait, " at centre"), title = trait) +
    ggplot2::theme_classic()

  if (x == "precip") p <- p + ggplot2::scale_x_log10()
  if (log_y) p <- p + ggplot2::scale_y_log10()
  p
}
