process_fit <- function(fit) {
  message("\t\033[32mComputing global mean predictions by log_annual_precipitation_mn...\033[0m")
  data_fit_global_mn <-
    estimate_means(fit, by = "log_annual_precipitation_mn") |>
    as_tibble() |>
    rename(estimate = Mean, conf.low = CI_low, conf.high = CI_high) |>
    mutate(
      taxon_name = "average_mn",
      annual_precipitation = exp(log_annual_precipitation_mn),
      across(c("estimate", "conf.low", "conf.high"), ~ exp(.x)),
    )

  message("\t\033[32mComputing global centered predictions by log_annual_precipitation_c...\033[0m")
  data_fit_global_c <-
    estimate_means(fit, by = "log_annual_precipitation_c") |>
    as_tibble() |>
    rename(estimate = Mean, conf.low = CI_low, conf.high = CI_high) |>
    mutate(
      taxon_name = "average_c",
      annual_precipitation = exp(6.7 + log_annual_precipitation_c), # 6.7 is the mean log annual precipitation across all species
      across(c("estimate", "conf.low", "conf.high"), ~ exp(.x))
    )

  message("\t\033[32mComputing individual predictions for all observations...\033[0m")
  data_fit_ind <-
    marginaleffects::predictions(fit) |>
    as_tibble() |>
    select(identifier, taxon_name, decimal_latitude, decimal_longitude, annual_precipitation, log_annual_precipitation, log_annual_precipitation_mn, log_annual_precipitation_c, value, log_value, estimate, conf.low, conf.high) |>
    mutate(
      log_residual = estimate - log_value,
      residual = exp(estimate) - exp(log_value),
      across(c("estimate", "conf.low", "conf.high"), ~ exp(.x))
    )

  data_fit_ind_range <-
    # use the original data but add predictions from the fitted model
    # =only keep first and last record for each species to show the range of rainfall
    data_fit_ind |>
    group_by(taxon_name) |>
    arrange(log_annual_precipitation) |>
    slice(c(1, n())) |>
    ungroup()

  message("\t\033[32mComputing predictions at species mean precipitation...\033[0m")
  data_fit_mn <-
    marginaleffects::predictions(fit,
      newdata = data_fit_ind |> mutate(log_annual_precipitation_c = 0) |> select(log_annual_precipitation_mn, log_annual_precipitation_c, taxon_name) |>
        distinct()
    ) |>
    as_tibble() |>
    ## Keep the log estimates but also convert to original scale for plotting
    mutate(
      annual_precipitation = exp(log_annual_precipitation_mn),
      # Convert back to original scale
      across(c("estimate", "conf.low", "conf.high"), ~ exp(.x))
    )

  message("\t\033[32mDone processing fit!\033[0m")

  list(
    global_mn = data_fit_global_mn,
    global_c = data_fit_global_c,
    ind = data_fit_ind,
    ind_range = data_fit_ind_range,
    mn = data_fit_mn
  )
}

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
