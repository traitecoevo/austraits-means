plot_fitted <- function(data_model, data_fitted) {
  data_model |>
    ggplot(aes(annual_precipitation, estimate, group = taxon_name)) +
    # raw data
    geom_point(aes(y = value, color = "Raw data"), alpha = 0.5, data = data_model, stroke = 0) +
    # species lines
    geom_line(aes(color = "Species fitted lines"), data = data_fitted$ind_range, alpha = 0.5, linewidth = 0.25) +
    # gloabl mn
    geom_line(aes(color = "Global mean (between species)"), data = data_fitted$global_mn) +
    # geom_ribbon(data = data_fit5_global_mn, aes(ymin = conf.low, ymax = conf.high, group = taxon_name), col = "red", alpha = 0.1, show.legend = FALSE) +
    # global c
    geom_line(aes(color = "Global mean (within species)"), data = data_fitted$global_c) +
    # geom_ribbon(data = data_fit5_global_c, aes(ymin = conf.low, ymax = conf.high), col = "blue", alpha = 0.1, show.legend = FALSE) +
    scale_color_manual(values = c(
      "Raw data" = "grey",
      "Species fitted lines" = "black",
      "Global mean (between species)" = "red",
      "Global mean (within species)" = "blue"
    )) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    labs(y = "height", x = "annual rainfall", color = NULL)
}

plot_fitted_highlight_taxa <- function(data_model, data_fitted, taxa, color_values = NULL, fill_values = NULL, y_labels = waiver()) {
  if (is.null(color_values)) {
    color_values <- c(
      "All raw data" = "grey",
      "All species fitted" = "black",
      "Focal species data" = "orange",
      "Focal species mean" = "red",
      "Focal species CI" = "purple"
    )
  }

  if (is.null(fill_values)) {
    fill_values <- c("Focal species CI" = "purple")
  }

  data_model_taxa <-
    data_model |>
    filter(taxon_name %in% taxa) |>
    mutate(taxon2 = taxon_name)

  data_fit_ind_taxa <-
    data_fitted$ind |>
    filter(taxon_name %in% taxa) |>
    arrange(taxon_name, annual_precipitation) |>
    mutate(taxon2 = taxon_name)

  data_fit_taxa_mn <-
    data_fitted$mn |>
    filter(taxon_name %in% taxa) |>
    mutate(taxon2 = taxon_name)

  p_fit_taxa <-
    data_model |>
    ggplot(aes(annual_precipitation, estimate, group = taxon_name)) +
    # raw data
    geom_point(aes(y = value, color = "All raw data"), alpha = 0.5, stroke = 0) +
    # fitted lines
    geom_line(aes(color = "All species fitted"), data = data_fitted$ind, alpha = 0.5, linewidth = 0.25) +
    # highlight certain species
    geom_point(data = data_model_taxa, aes(y = value, group = taxon2, color = "Focal species data"), size = 1) +
    geom_point(data = data_fit_taxa_mn, aes(x = annual_precipitation, group = taxon2, color = "Focal species mean"), alpha = 0.5) +
    geom_line(data = data_fit_ind_taxa, aes(group = taxon2, color = "Focal species mean"), alpha = 1) +
    geom_ribbon(data = data_fit_ind_taxa, aes(ymin = conf.low, ymax = conf.high, group = taxon2, fill = "Focal species CI"), alpha = 0.2) +
    scale_color_manual(values = color_values) +
    scale_fill_manual(values = fill_values) +
    scale_x_log10() +
    scale_y_log10(labels = y_labels) +
    theme_classic() +
    labs(y = "Height (m)", x = "Annual rainfall (mm)", color = NULL, fill = NULL)

  p_fit_taxa
}
