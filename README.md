# Vesk-austraits-means

Quarto and R workflow for analysing Australian plant height variation with AVH occurrence data and climate covariates.

## File tree

```text
.
├── Analysis.qmd
├── Data_height_AVH.qmd
├── Data_height_AVH_extra.qmd
├── Dataset_construction.qmd
├── climate_centres.qmd
├── trait_species_means.qmd
├── demo.qmd
├── Falster-trait_means.Rproj
├── R/
│   ├── data/
│   │   └── load_data_height_avh.R
│   ├── locate_medoid_2d.R
│   ├── process_fitted_model.R
│   └── species_means.R
├── data/
│   ├── australia.tif
│   ├── climate/
│   └── data-ausplants/
├── data-processed/
│   ├── climate_centres/
│   ├── data_height_climate.parquet
│   ├── data_traits_climate-height.parquet
│   ├── data_traits_climate-height-woody.csv
│   ├── data_gbif_cleaned.parquet
│   └── fits_height_climate.rds
├── outputs/
│   ├── fits/
│   ├── height/
│   └── mixed models/
├── austraits-means.code-workspace
└── .gitignore
```
