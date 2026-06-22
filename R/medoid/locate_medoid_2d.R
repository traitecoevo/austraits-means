
# X is an n x 2 matrix or data.frame of points (columns: x, y)
locate_medoid_2d <- function(data_taxon, method = "euclidean", cull_percentage = 0.01, max_cull = 30) {
  
  data_taxon <- data_taxon |> 
    filter(!is.na(annual_precipitation), !is.na(average_temperature))
  
  # distance matrix
  D <- data_taxon |>
    select(log_annual_precipitation_centered, average_temperature_centered) |>
    dist(method = "euclidean") |> 
    as.matrix()

  # Drop the points furthest away from most other points (in environmental space), 
  # as some of these are likely errors
  if(nrow(data_taxon) > 10) {
    n_cull <- (nrow(data_taxon) * cull_percentage) |> ceiling() |> min(max_cull)
    s <- D |> rowSums()
    idx <- order(s, decreasing = TRUE)[seq_len(n_cull)]

    # Update dataset
    data_culled <- data_taxon[idx,]
    
    data_taxon <- data_taxon[-idx,]
    D <- D[-idx, -idx]
   } else {
      data_culled <- data_taxon |> filter(taxon_name == "nothing")
   }
  
  # Of reminaing points, medoid is the one closts to most points
  s <- D |> rowSums()
  idx <- s |> which.min() |> first()

  # return that point
  list(
    medoid = data_taxon[idx, , drop = FALSE],
    data_taxon = data_taxon,
    data_culled = data_culled
  )
}

plot_taxon <- function(taxon, data_taxon, medoid_pt, data_culled, file = NULL) {


  p_main <- 
    ggplot(data_taxon, aes(annual_precipitation, average_temperature)) +
    geom_point(color = "#2C6E8E") +
    geom_point(data = data_culled, col = "orange", alpha = 0.5) +
    geom_point(data = medoid_pt, color = "red", size = 5) +
    geom_hline(yintercept = medoid_pt$average_temperature, color = "red", linetype = "dashed") +
    geom_vline(xintercept = medoid_pt$annual_precipitation, color = "red", linetype = "dashed") +
    scale_x_log10() +
    labs(
      x = "Annual precipitation (mm)",
      y = "Average temperature (deg C)"
    ) +
    theme_minimal()

  x_limits <- 10^layer_scales(p_main)$x$range$range
  y_limits <- layer_scales(p_main)$y$range$range

  p_top <- ggplot(data_taxon, aes(annual_precipitation)) +
    geom_histogram(fill = "#2C6E8E", bins = 30) +
    labs(y = "Frequency") +
    geom_vline(xintercept = medoid_pt$annual_precipitation, color = "red", linetype = "dashed") +
    scale_x_log10(limits = x_limits) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )

  p_right <- ggplot(data_taxon, aes(average_temperature)) +
    geom_histogram(fill = "#2C6E8E", bins = 30) +
    labs(y = "Frequency") +
    geom_vline(xintercept = medoid_pt$average_temperature, color = "red", linetype = "dashed") +
    scale_x_continuous(limits = y_limits) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    )

  subtitle <- sprintf("Precip = %s mm, Temp = %0.2f degC\nLat = %s deg, Long = %s deg", 
    medoid_pt$annual_precipitation,
    medoid_pt$average_temperature,
    medoid_pt$decimal_latitude,
    medoid_pt$decimal_longitude
  )

  title <- sprintf("%s (n=%s)", taxon, nrow(data_taxon))

  data_ozmaps <- ozmaps::ozmap_country
  
  p_map <- ggplot() +
    geom_sf(data = data_ozmaps, fill = NA, color = "black", linewidth = 0.2) +
    geom_point(data = data_taxon, aes(x = decimal_longitude, y = decimal_latitude), alpha = 0.5, col = "#2C6E8E") +
    geom_point(data = data_culled, aes(x = decimal_longitude, y = decimal_latitude), alpha = 0.5, col = "orange") +
    geom_point(data = medoid_pt, aes(x = decimal_longitude, y = decimal_latitude), color = "red", size = 3) +
    coord_sf(ylim = c(-45, -10), xlim = c(110, 155), expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = "", y = "") +
    theme_minimal()

  p_all <- p_map / (p_top + plot_spacer() + p_main + p_right +
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4)))

  if (!is.null(file)) {
    ggsave(file, p_all, width = 8, height = 10)
  }
}

process_taxon_medoid <- function(taxon, output_dir = "~/Downloads/climate_centres/", overwrite = FALSE, max_n = 20000, cull_percentage = 0.01, data_file = "data-processed/data_gbif_cleaned.parquet") {

  output_data <- file.path(output_dir, "data", sprintf("%s.rds", taxon))
  output_plot <- file.path(output_dir, "figures", sprintf("%s.png", taxon))

  if(!overwrite & all(file.exists(c(output_data, output_plot)))) return(invisible());
  
  data_taxon <- arrow::open_dataset(data_file) |>
    filter(taxon_name == taxon) |>
    select(taxon_name, decimal_longitude, decimal_latitude, annual_precipitation, average_temperature, log_annual_precipitation_centered, average_temperature_centered, ) |>
    collect() |> 
    distinct() |> 
    filter(!is.na(annual_precipitation), !is.na(average_temperature))
  
  # Take a subsample of very large datasets. Need as Calculating distance matrix crashes is memory intensive and crashes my MBP for large numbers of points (e.g 50k)

  if(nrow(data_taxon) > max_n) {
    set.seed(123)
    data_taxon <- data_taxon |>
      sample_n(max_n)
   }

  # Calculate the medoid, also identify outliers
  med_out <- locate_medoid_2d(data_taxon, cull_percentage = cull_percentage)
  
  saveRDS(med_out$medoid, file = output_data)
  plot_taxon(taxon, med_out$data_taxon, med_out$medoid,  med_out$data_culled, file = output_plot)

  rm(data_taxon, med_out)
  return(invisible())
}
