library(tidyverse)
library(arrow)
library(glmmTMB)
library(patchwork)
library(tictoc)

source("R/init.R")

fits <- readRDS("data-processed/fits_height_climate.rds")
data_traits_climate <- arrow::read_parquet("data-processed/data_height_climate.parquet")

# Get 4 test species
taxa <- fits$random_effects_r |>
  arrange(desc(g_c)) |>
  slice_head(n = 4) |>
  pull(taxon_name)

print(paste("Testing with species:", paste(taxa, collapse = ", ")))

# Setup (same as in generate_species_panels)
gbif_path <- "data-processed/data_gbif_cleaned.parquet"
out_dir_species <- "outputs/height/species"
dir.create(out_dir_species, recursive = TRUE, showWarnings = FALSE)

country_sf <- ozmaps::ozmap_country
states_sf <- ozmaps::ozmap_states

gbif_background <- arrow::open_dataset(gbif_path) |>
  dplyr::select(annual_precipitation, average_temperature) |>
  dplyr::filter(!is.na(annual_precipitation), !is.na(average_temperature)) |>
  dplyr::collect()

gbif_background <- gbif_background |>
  dplyr::slice_sample(n = min(200000, nrow(gbif_background)))

data_model_bg <- fits$data_model
data_fitted_bg <- list(
  ind = fits$data_fitted_r$ind,
  ind_range = fits$data_fitted_r$ind_range,
  mn = fits$data_fitted_r$mn
)

# Build background rain plot once before the loop
tic("Build rain background (one-off)")
rain_bg_plot <- build_rain_background_plot(data_model_bg, data_fitted_bg)
toc()

# Run sequentially with timing
tic("Total time for 4 species")
for (i in seq_along(taxa)) {
  taxon <- taxa[i]
  tic(paste("Species", i, "-", taxon))
  
  tic("  - GBIF data loading")
  taxon_gbif <- arrow::open_dataset(gbif_path) |>
    dplyr::filter(taxon_name == taxon) |>
    dplyr::select(
      taxon_name, decimal_longitude, decimal_latitude,
      annual_precipitation, average_temperature
    ) |>
    dplyr::collect() |>
    dplyr::distinct() |>
    dplyr::filter(
      !is.na(decimal_longitude), !is.na(decimal_latitude),
      !is.na(annual_precipitation), !is.na(average_temperature)
    ) |>
    dplyr::mutate(source = "GBIF")
  toc()
  
  taxon_avh <- data_traits_climate |>
    dplyr::filter(taxon_name == taxon) |>
    dplyr::select(
      taxon_name, decimal_longitude, decimal_latitude,
      value,log_value, annual_precipitation, average_temperature,
      log_annual_precipitation, log_annual_precipitation_mn, log_annual_precipitation_c,
      average_temperature_mn, average_temperature_c, state
    ) |>
    dplyr::filter(
      !is.na(decimal_longitude), !is.na(decimal_latitude),
      !is.na(annual_precipitation), !is.na(average_temperature)
    ) |>
    dplyr::mutate(source = "AVH")
  
  if (nrow(taxon_avh) < 3 || nrow(taxon_gbif) < 10) {
    message("Skipping (insufficient data)")
    next
  }
  
  
  tic("  - Plot 1 (map)")
  p1 <- plot_panel_map(taxon_gbif, taxon_avh, country_sf, states_sf)
  toc()
  
  tic("  - Plot 2 (climate)")
  p2 <- plot_panel_climate(taxon_gbif, taxon_avh, gbif_background)
  toc()
  
  tic("  - Plot 4 (state)")
  p4 <- plot_panel_state(taxon_avh)
  toc()

  tic("  - Plot 5 (rain)")
  p5 <- plot_panel_rain(taxon, data_fitted_bg, taxon_avh, rain_bg_plot)
  toc()

  tic("  - Plot 6 (temp placeholder)")
  p6 <- plot_panel_temperature_placeholder()
  toc()

  tic("  - Plot 7 (height heatmap)")
  p7 <- plot_panel_height_heatmap(taxon_gbif, taxon_avh, fits$fit_rt)
  toc()
  
  tic("  - Combine plots and save")
  title_text <- sprintf("%s (AVH n=%s, GBIF n=%s)", taxon, nrow(taxon_avh), nrow(taxon_gbif))
  p_all <-
    patchwork::wrap_plots(p1, p2, p4, p7, p5, p6, ncol = 2, byrow = TRUE) +
    patchwork::plot_layout(widths = c(1, 1), heights = c(1, 1, 1), guides = "collect") +
    patchwork::plot_annotation(
      title = title_text,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 18, face = "bold"))
    )
  
  p_all <- p_all &
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      legend.position = "right"
    )
  
  file_out <- file.path(out_dir_species, paste0(clean_species_filename(taxon), ".png"))
  ggplot2::ggsave(file_out, p_all, width = 14, height = 15, dpi = 120)
  toc()
  
  toc()
}
toc()
