# Vesk-austraits-means

Estimating **national means of plant traits at the centre of each species'
climate distribution** for Australia, from AusTraits + AVH herbarium records and
WorldClim covariates. Part of the ARDC-funded AusTraits programme; also the basis
for a paper reporting these species means.

A species' trait value varies across its range, and trait data are sampled
unevenly across that range. We therefore estimate, per trait, the species mean at
the **centre of the species' climate distribution** (a 2D precipitation–temperature
medoid), using a hierarchical within-group-centred model. Two interchangeable
fitting backends are available: `tmb` (fast) and `brms` (full posterior).

See `CLAUDE.md` for conventions, `refactor_plan.md` for the design/roadmap, and
`demo.qmd` for the statistical motivation.

## Workflow

Two Quarto notebooks, which only orchestrate functions in `R/`:

1. `dataset_construction.qmd` — climate covariates → GBIF + medoid climate centres
   (+ reasonableness checks) → trait sources (AVH height, AusTraits) → cached
   **unified record-level dataset**.
2. `analysis.qmd` — loads the unified dataset → fits one model per trait →
   exports per-species centre estimates and figures.

Run order: `dataset_construction.qmd` then `analysis.qmd`. Expensive steps are
cached behind `if (update || !file.exists(...))`; set `update <- TRUE` to rebuild.

## Layout

```text
.
├── dataset_construction.qmd     # data pipeline -> unified dataset
├── analysis.qmd                 # per-trait models -> estimates + figures
├── trait_species_means.qmd      # round-2 simple means (low-coverage traits, WIP)
├── demo.qmd                     # within-group-centring method illustration
├── CLAUDE.md  refactor_plan.md
├── R/
│   ├── init.R                   # source() this to load the whole library
│   ├── data/                    # build_dataset, extract_climate, assign_state, load_data_height_avh
│   ├── medoid/                  # locate_medoid_2d, check_medoid
│   ├── model/                   # model_formula, fit_trait_climate_model, run_trait_models
│   ├── process/                 # process_fit, species_centres, export_centres
│   ├── plot/                    # plot_fitted, plot_trait_centres
│   ├── panels/                  # species diagnostic panels
│   ├── means/                   # simple_means (round-2 fallback) + species_means
│   ├── utils/                   # save_figure (png+pdf), parallel_setup (multisession)
│   └── scripts/                 # run-scripts (NOT sourced by init.R)
├── data/                        # inputs (git-ignored): AVH, GBIF, WorldClim, AusTraits
├── data-processed/              # cached intermediates (git-ignored)
├── outputs/                     # figures + per-trait response CSVs (git-ignored)
└── archive/                     # superseded notebooks (retired NVIS report, etc.)
```
