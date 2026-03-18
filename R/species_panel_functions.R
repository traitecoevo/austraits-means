clean_species_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

species_panel_palette <- function() {
  list(
    avh = "#fb0019",
    predicted = "#ffd014",
    gbif = "#00536d",
    all_taxa = "black",
    error = "#fdc809",
    predicted_CI = "#f8ca8980",
    raw = "#B7B7B7"
  )
}

predict_species_surface <- function(fit_rt, taxon, taxon_avh, taxon_gbif, n_grid = 40) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package `glmmTMB` is required for prediction.")
  }

  p_rng <- range(taxon_gbif$annual_precipitation, na.rm = TRUE)
  t_rng <- range(taxon_gbif$average_temperature, na.rm = TRUE)

  if (!all(is.finite(c(p_rng, t_rng))) || diff(log(p_rng)) <= 0 || diff(t_rng) <= 0) {
    return(tibble::tibble(
      annual_precipitation = numeric(0),
      average_temperature = numeric(0),
      pred_height = numeric(0)
    ))
  }

  p_seq <- exp(seq(log(p_rng[1]), log(p_rng[2]), length.out = n_grid))
  t_seq <- seq(t_rng[1], t_rng[2], length.out = n_grid)

  log_p_mn <- unique(taxon_avh$log_annual_precipitation_mn) |> dplyr::first()
  t_mn <- unique(taxon_avh$average_temperature_mn) |> dplyr::first()

  newdata <- tidyr::expand_grid(
    annual_precipitation = p_seq,
    average_temperature = t_seq
  ) |>
    dplyr::mutate(
      log_annual_precipitation = log(annual_precipitation),
      log_annual_precipitation_mn = log_p_mn,
      average_temperature_mn = t_mn,
      log_annual_precipitation_c = log_annual_precipitation - log_annual_precipitation_mn,
      average_temperature_c = average_temperature - average_temperature_mn,
      taxon_name = taxon
    )

  newdata$pred_height <- stats::predict(fit_rt, newdata = newdata, type = "response", re.form = NULL)
  newdata
}

clip_surface_to_avh_hull <- function(surf, taxon_avh) {
  if (nrow(surf) == 0) {
    return(surf)
  }

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package `sf` is required for hull-based contour clipping.")
  }

  avh_pts <- taxon_avh |>
    dplyr::select(annual_precipitation, average_temperature) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(annual_precipitation), !is.na(average_temperature))

  if (nrow(avh_pts) < 3) {
    p_rng <- range(taxon_avh$annual_precipitation, na.rm = TRUE)
    t_rng <- range(taxon_avh$average_temperature, na.rm = TRUE)

    return(
      surf |>
        dplyr::filter(
          annual_precipitation >= p_rng[1],
          annual_precipitation <= p_rng[2],
          average_temperature >= t_rng[1],
          average_temperature <= t_rng[2]
        )
    )
  }

  hull_idx <- grDevices::chull(avh_pts$annual_precipitation, avh_pts$average_temperature)
  hull_coords <- avh_pts[hull_idx, ]
  hull_coords <- dplyr::bind_rows(hull_coords, hull_coords[1, , drop = FALSE])

  hull_poly <- sf::st_polygon(list(as.matrix(hull_coords))) |>
    sf::st_sfc()

  surf_sf <- sf::st_as_sf(surf, coords = c("annual_precipitation", "average_temperature"), remove = FALSE)
  inside <- sf::st_within(surf_sf, hull_poly, sparse = FALSE)[, 1]

  surf[inside, , drop = FALSE]
}

theme_panel_axes <- function() {
  ggplot2::theme(
    axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.4),
    panel.border = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "plain")
  )
}

sf_to_path_df <- function(x) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package `sf` is required for converting boundaries to paths.")
  }

  coords <- sf::st_coordinates(sf::st_cast(sf::st_geometry(x), "MULTILINESTRING"))
  coords <- as.data.frame(coords)

  group_cols <- intersect(c("L1", "L2", "L3"), names(coords))
  if (length(group_cols) == 0) {
    coords$group_id <- 1
  } else {
    coords$group_id <- apply(coords[, group_cols, drop = FALSE], 1, paste, collapse = "-")
  }

  tibble::as_tibble(coords) |>
    dplyr::rename(x = X, y = Y)
}

plot_panel_map <- function(taxon_gbif, taxon_avh, country_sf, states_sf) {
  pal <- species_panel_palette()

  country_path <- sf_to_path_df(country_sf)
  states_path <- sf_to_path_df(states_sf)

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = country_path,
      ggplot2::aes(x = x, y = y, group = group_id),
      color = "grey35",
      linewidth = 0.25
    ) +
    ggplot2::geom_path(
      data = states_path,
      ggplot2::aes(x = x, y = y, group = group_id),
      color = "grey50",
      linewidth = 0.2
    ) +
    ggplot2::geom_point(
      data = taxon_gbif,
      ggplot2::aes(decimal_longitude, decimal_latitude, color = source),
      alpha = 0.6,
      size = 0.6
    ) +
    ggplot2::geom_point(
      data = taxon_avh,
      ggplot2::aes(decimal_longitude, decimal_latitude, color = source),
      alpha = 1,
      size = 1
    ) +
    ggplot2::scale_color_manual(
      values = c("All Australian taxa" = pal$all_taxa, "GBIF" = pal$gbif, "AVH" = pal$avh),
      limits = c("All Australian taxa", "GBIF", "AVH"),
      name = "Points 1",
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::coord_cartesian(ylim = c(-45, -10), xlim = c(110, 155), expand = FALSE) +
    ggplot2::labs(title = "a) Species Distribution", x = "Longitude", y = "Latitude") +
    ggplot2::theme_minimal() +
    theme_panel_axes() +
    ggplot2::theme(legend.position = "none")
}

plot_panel_climate <- function(taxon_gbif, taxon_avh, surf, gbif_background) {
  pal <- species_panel_palette()
  surf <- clip_surface_to_avh_hull(surf, taxon_avh)

  gbif_background <- gbif_background |>
    dplyr::mutate(source = "All Australian taxa")

  ggplot2::ggplot() +
    # Background points (rasterized for speed)
    ggplot2::geom_point(
      data = gbif_background,
      ggplot2::aes(annual_precipitation, average_temperature, color = source),
      alpha = 0.3,
      size = 0.8,
      shape = 16
    ) +
    # Hex bin for GBIF density (faster than geom_point for large data)
    {if (nrow(taxon_gbif) > 100) ggplot2::geom_hex(data = taxon_gbif, ggplot2::aes(annual_precipitation, average_temperature), bins = 20) else ggplot2::geom_blank()} +
    # Focal species points
    ggplot2::geom_point(
      data = taxon_gbif,
      ggplot2::aes(annual_precipitation, average_temperature, color = source),
      alpha = 1,
      size = 1.2
    ) +
    ggplot2::geom_point(
      data = taxon_avh,
      ggplot2::aes(annual_precipitation, average_temperature, color = source),
      alpha = 1,
      size = 1.2
    ) +
    # Contours (only if surface has data)
    {if (nrow(surf) > 0) ggplot2::geom_contour(data = surf, ggplot2::aes(annual_precipitation, average_temperature, z = pred_height), color = pal$predicted, linewidth = 0.4, bins = 10)} +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(
      values = c(
        "All Australian taxa" = pal$raw,
        "GBIF" = pal$gbif,
        "AVH" = pal$avh),
      name = "Points"
    ) +
    ggplot2::scale_fill_distiller(palette = "Blues", direction = 1, trans = "sqrt") +
    ggplot2::labs(
      title = "b) Environmental space",
      x = "Annual precipitation (mm)",
      y = "Average temperature (deg C)",
      fill = "GBIF count"
    ) +
    ggplot2::theme_minimal() +
    theme_panel_axes()
}

plot_panel_state <- function(taxon_avh) {
  pal <- species_panel_palette()
  avh_state <- taxon_avh |>
    dplyr::mutate(
      state = ifelse(is.na(state) | state == "", "Other", state),
      state = dplyr::case_when(
        state %in% c("New South Wales", "NSW") ~ "NSW",
        state %in% c("Victoria", "VIC") ~ "VIC",
        state %in% c("Queensland", "QLD") ~ "QLD",
        state %in% c("South Australia", "SA") ~ "SA",
        state %in% c("Western Australia", "WA") ~ "WA",
        state %in% c("Tasmania", "TAS") ~ "TAS",
        state %in% c("Northern Territory", "NT") ~ "NT",
        state %in% c("Australian Capital Territory", "ACT") ~ "ACT",
        TRUE ~ state
      )
    )

  state_levels <- sort(c("ACT", "NSW", "NT", "QLD", "SA", "TAS", "VIC", "WA", "Other"))
  avh_state <- avh_state |>
    dplyr::mutate(state = factor(state, levels = state_levels))

  # Model-by-state means on log scale, then convert to height scale with 95% CI.
  observed_levels <- avh_state |>
    dplyr::filter(!is.na(state)) |>
    dplyr::pull(state) |>
    unique() |>
    as.character()

  if (length(observed_levels) >= 2) {
    mod <- stats::lm(log(value) ~ state, data = avh_state)
  } else {
    mod <- stats::lm(log(value) ~ 1, data = avh_state)
  }

  state_est <- purrr::map_dfr(
    state_levels,
    \(s) {
      if (!(s %in% observed_levels)) {
        return(tibble::tibble(
          state = factor(s, levels = state_levels),
          log_est = NA_real_,
          log_lo = NA_real_,
          log_hi = NA_real_,
          est = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_
        ))
      }

      pred <- stats::predict(
        mod,
        newdata = if (length(observed_levels) >= 2) {
          data.frame(state = factor(s, levels = observed_levels))
        } else {
          data.frame(dummy = 1)
        },
        se.fit = TRUE
      )

      tibble::tibble(
        state = factor(s, levels = state_levels),
        log_est = as.numeric(pred$fit),
        log_lo = as.numeric(pred$fit - stats::qnorm(0.975) * pred$se.fit),
        log_hi = as.numeric(pred$fit + stats::qnorm(0.975) * pred$se.fit),
        est = exp(log_est),
        conf.low = exp(log_lo),
        conf.high = exp(log_hi)
      )
    }
  )

  ggplot2::ggplot(avh_state, ggplot2::aes(x = state, y = value)) +
    ggplot2::geom_jitter(width = 0.15, alpha = 1, size = 1.8, color = pal$avh) +
    ggplot2::geom_errorbar(
      data = state_est,
      inherit.aes = FALSE,
      mapping = ggplot2::aes(x = state, ymin = conf.low, ymax = conf.high),
      width = 0.15,
      color = pal$error,
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = state_est,
      inherit.aes = FALSE,
      mapping = ggplot2::aes(x = state, y = est),
      color = pal$predicted,
      size = 3.4
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "c) Difference by State", x = "State", y = "Height (m)") +
    ggplot2::theme_minimal() +
    theme_panel_axes() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey88", linewidth = 0.35)
    )
}

plot_panel_rain <- function(taxon, data_fitted_bg, taxon_avh, rain_bg_plot) {
  pal <- species_panel_palette()

  # Use pre-computed predictions from the fitted model (avoids expensive glmmTMB predict() call)
  data_fit_ind_taxa <- data_fitted_bg$ind |>
    dplyr::filter(taxon_name == taxon) |>
    dplyr::arrange(annual_precipitation) |>
    dplyr::mutate(taxon2 = taxon_name)

  slope_fit <- stats::lm(log_value ~ log(annual_precipitation), data = data_fit_ind_taxa)
  slope <- stats::coef(slope_fit)[[2]]
  slope_ci <- stats::confint(slope_fit)[2, ]

  # R² between observed and model-fitted log height
  r2 <- stats::cor(data_fit_ind_taxa$log_value, log(data_fit_ind_taxa$estimate), use = "complete.obs")^2

  data_fit_taxa_mn <- data_fitted_bg$mn |>
    dplyr::filter(taxon_name == taxon) |>
    dplyr::mutate(taxon2 = taxon_name)

  rain_bg_plot +
    ggplot2::geom_point(
      data = data_fit_ind_taxa,
      ggplot2::aes(annual_precipitation, value, group = taxon2, color = "Focal species data"),
      size = 1
    ) +
    ggplot2::geom_point(
      data = data_fit_taxa_mn,
      ggplot2::aes(annual_precipitation, estimate, group = taxon2, color = "Focal species mean"),
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      data = data_fit_ind_taxa,
      ggplot2::aes(annual_precipitation, estimate, group = taxon2, color = "Focal species mean"),
      alpha = 1
    ) +
    ggplot2::geom_ribbon(
      data = data_fit_ind_taxa,
      ggplot2::aes(annual_precipitation, ymin = conf.low, ymax = conf.high, group = taxon2),
      fill = pal$predicted,
      alpha = 0.25
    ) +
    ggplot2::labs(
      title = "d) Rain fall response",
      subtitle = sprintf("Slope = %.2f [%.1f-%.1f], R2 = %.2f", slope, slope_ci[1], slope_ci[2], r2)
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

build_rain_background_plot <- function(data_model_bg, data_fitted_bg) {
  pal <- species_panel_palette()
  color_values <- c(
    "All raw data" = pal$raw,
    "All species fitted" = pal$all_taxa,
    "Focal species data" = pal$avh,
    "Focal species mean" = pal$predicted,
    "Focal species CI" = pal$predicted_CI
  )
  fill_values <- c("Focal species CI" = pal$predicted_CI)

  ggplot2::ggplot(data_model_bg, ggplot2::aes(annual_precipitation, group = taxon_name)) +
    # Use ind_range (2 rows per species: min/max precip) instead of ind (all points) — much less data
    ggplot2::geom_point(
      ggplot2::aes(y = value, color = "All raw data"),
      alpha = 0.5, stroke = 0, size = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = estimate, color = "All species fitted"),
      data = data_fitted_bg$ind_range,
      alpha = 0.5, linewidth = 0.25
    ) +
    ggplot2::scale_color_manual(values = color_values, drop = FALSE) +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10(labels = scales::label_number(accuracy = 0.01, trim = TRUE)) +
    ggplot2::theme_classic() +
    ggplot2::labs(y = "Height (m)", x = "Annual rainfall (mm)", color = NULL, fill = NULL)
}

build_species_plot <- function(
  taxon,
  fit_rt,
  data_model_bg,
  data_fitted_bg,
  data_traits_climate,
  gbif_path,
  out_dir,
  fit_r,
  country_sf,
  states_sf,
  gbif_background,
  rain_bg_plot
) {
  taxon_avh <- data_traits_climate |>
    dplyr::filter(taxon_name == taxon) |>
    dplyr::select(
      taxon_name, decimal_longitude, decimal_latitude,
      value, annual_precipitation, average_temperature,
      log_value,
      log_annual_precipitation_mn, log_annual_precipitation_c,
      average_temperature_mn, state
    ) |>
    dplyr::filter(
      !is.na(decimal_longitude), !is.na(decimal_latitude),
      !is.na(annual_precipitation), !is.na(average_temperature)
    ) |>
    dplyr::mutate(source = "AVH")

  if (nrow(taxon_avh) < 3) {
    message("Skipping ", taxon, " (insufficient AVH records)")
    return(invisible(NULL))
  }

  # Each worker opens GBIF data lazily from parquet
  gbif_full <- arrow::open_dataset(gbif_path) |>
    dplyr::collect()

  taxon_gbif <- gbif_full |>
    dplyr::filter(taxon_name == taxon) |>
    dplyr::select(
      taxon_name, decimal_longitude, decimal_latitude,
      annual_precipitation, average_temperature
    ) |>
    dplyr::distinct() |>
    dplyr::filter(
      !is.na(decimal_longitude), !is.na(decimal_latitude),
      !is.na(annual_precipitation), !is.na(average_temperature)
    ) |>
      dplyr::mutate(source = "GBIF")

  if (nrow(taxon_gbif) < 10) {
    message("Skipping ", taxon, " (insufficient GBIF records)")
    return(invisible(NULL))
  }

  surf <- predict_species_surface(fit_rt, taxon, taxon_avh, taxon_gbif)

  p1 <- plot_panel_map(taxon_gbif, taxon_avh, country_sf, states_sf)

  p2 <- plot_panel_climate(taxon_gbif, taxon_avh, surf, gbif_background)
  
  p3 <- plot_panel_state(taxon_avh)
  
  p4 <- plot_panel_rain(taxon, data_fitted_bg, taxon_avh, rain_bg_plot)

  title_text <- sprintf("%s (AVH n=%s, GBIF n=%s)", taxon, nrow(taxon_avh), nrow(taxon_gbif))

  p_all <-
    patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2, byrow = TRUE) +
    patchwork::plot_layout(widths = c(1, 1), heights = c(1, 1), guides = "collect") +
    patchwork::plot_annotation(
      title = title_text,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 18, face = "bold"))
    )

  p_all <- p_all &
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      legend.position = "right"
    )

  file_out <- file.path(out_dir, paste0(clean_species_filename(taxon), ".png"))
  ggplot2::ggsave(file_out, p_all, width = 14, height = 10, dpi = 120)

  message("Saved ", file_out)
  invisible(NULL)
}

generate_species_panels <- function(
  data_traits_climate,
  fits,
  taxa,
  gbif_path = "data-processed/data_gbif_cleaned.parquet",
  out_dir_species = "outputs/height/species",
  workers_n = max(1, min(8, parallelly::availableCores() - 1)),
  future_strategy = c("auto", "multisession", "sequential", "multicore"),
  max_global_size = 4 * 1024^3,
  seed = 123,
  n_background = 200000
) {
  future_strategy <- match.arg(future_strategy)

  if (missing(taxa) || length(taxa) == 0) {
    stop("Please provide a non-empty character vector to `taxa`.")
  }

  taxa <- unique(as.character(taxa))

  dir.create(out_dir_species, recursive = TRUE, showWarnings = FALSE)

  country_sf <- ozmaps::ozmap_country
  states_sf <- ozmaps::ozmap_states

  gbif_background <- arrow::open_dataset(gbif_path) |>
    dplyr::select(annual_precipitation, average_temperature) |>
    dplyr::filter(!is.na(annual_precipitation), !is.na(average_temperature)) |>
    dplyr::collect()
  
  gbif_background <- gbif_background |>
    dplyr::slice_sample(n = min(n_background, nrow(gbif_background)))

  set.seed(seed)
  data_model_bg <- fits$data_model

  data_fitted_bg <- list(
    ind = fits$data_fitted_r$ind,
    ind_range = fits$data_fitted_r$ind_range,
    mn = fits$data_fitted_r$mn
  )

  taxa_test <- taxa

  # Build background rain plot once — reused across all species
  rain_bg_plot <- build_rain_background_plot(data_model_bg, data_fitted_bg)

  options(future.globals.maxSize = max_global_size)

  if (future_strategy == "auto") {
    # Forked multicore is unstable on macOS with graphics/ObjC-loaded libs.
    if (identical(Sys.info()[["sysname"]], "Darwin")) {
      future_strategy <- "multisession"
    } else if (.Platform$OS.type == "unix") {
      future_strategy <- "multicore"
    } else {
      future_strategy <- "multisession"
    }
  }

  if (future_strategy == "sequential") {
    future::plan(future::sequential)
  } else if (future_strategy == "multisession") {
    if (identical(Sys.info()[["sysname"]], "Darwin")) {
      options(future.fork.enable = FALSE)
    }
    future::plan(future::multisession, workers = workers_n)
  } else if (future_strategy == "multicore") {
    future::plan(future::multicore, workers = workers_n)
  }

  on.exit(future::plan(future::sequential), add = TRUE)

  furrr::future_walk(
    taxa_test,
    \(taxon) {
      tryCatch(
        build_species_plot(
          taxon = taxon,
          fit_rt = fits$fit_rt,
          fit_r = fits$fit_r,
          data_model_bg = data_model_bg,
          data_fitted_bg = data_fitted_bg,
          data_traits_climate = data_traits_climate,
          gbif_path = gbif_path,
          out_dir = out_dir_species,
          country_sf = country_sf,
          states_sf = states_sf,
          gbif_background = gbif_background,
          rain_bg_plot = rain_bg_plot
        ),
        error = \(e) message("Failed for ", taxon, ": ", conditionMessage(e))
      )
    },
    .options = furrr::furrr_options(seed = TRUE, packages = c("glmmTMB", "ggplot2", "dplyr", "tidyr", "tibble", "purrr", "arrow", "patchwork", "sf"))
  )

  n_png <- list.files(out_dir_species, pattern = "\\.png$", full.names = FALSE) |> length()
  message("Generated ", n_png, " png files in ", out_dir_species)
  invisible(n_png)
}
