library(tidyverse)
library(arrow)
library(glmmTMB)
library(patchwork)
library(furrr)
library(future)
library(parallelly)
library(ozmaps)
library(sf)
library(scales)

source("R/init.R")

fits <- readRDS("data-processed/fits_height_climate.rds")
data_traits_climate <- arrow::read_parquet("data-processed/data_height_climate.parquet")

taxa <- fits$random_effects_r |>
  dplyr::pull(taxon_name) |> 
  unique() |> sample(50)

generate_species_panels(
  data_traits_climate = data_traits_climate,
  fits = fits,
  taxa = taxa[1:4],
  gbif_path = "data-processed/data_gbif_cleaned.parquet",
  out_dir_species = "outputs/height/species",
  future_strategy = "multisession"
)
