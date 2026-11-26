# austraits-means

Species-level and community-level trait means for the Australian flora, built from [AusTraits](https://austraits.org) with WorldClim climate and GBIF occurrences. Part of the AusTraits programme.

The repository publishes, in a Zenodo record of its own and for a short list of traits (terms fixed in [`terminology.md`](terminology.md)):

1. **Modelled means**: per species × trait, the geometric mean at the centre of the species' climate distribution, from a hierarchical model in which each species' climate is centred on its medoid (the most central point of its GBIF occurrences in precipitation–temperature space; `docs/within-species-centring.qmd` explains why). Fitted with TMB (`glmmTMB`, maximum likelihood); a Bayesian `brms` backend is optional.
2. **Modelled community statistics**: for a location, predicted from its climate, the woody and non-woody community mean and community range. In development: the model they need is the next step in [`roadmap.md`](roadmap.md).

The **observed means** (simple species means) are published by austraits.build with each AusTraits release, and read here. The *combined species means* table (the modelled mean where a species has one, otherwise its observed mean), and the community statistics (first version) are built but not published: they are inputs or work in progress.

[`roadmap.md`](roadmap.md) records the milestones, the decisions and their reasons, and the invariants; [`AGENTS.md`](AGENTS.md) the conventions for working in the code.

## The products

Every species means table carries `estimate_type` (`observed_mean` or `modelled_mean`) and `method` (`geometric` or `arithmetic` for an observed mean; `rtmb`, `tmb` or `brms` for a modelled mean), and comes as four files: `.parquet`, `.csv.gz`, `_dictionary.csv` (every column defined) and `_metadata.json` (provenance: source build, commit, data basis). Each product carries the version tag of the AusTraits release it derives from.

What the pipeline writes, published or not (status in [`roadmap.md`](roadmap.md)):

| Product | Files (under `outputs/`) | Documentation |
|---|---|---|
| Observed means, geometric method | Published with each AusTraits release by austraits.build (`austraits-<version>-observed-means.*`); read here, not written | [`docs/observed-means-geometric-method.md`](docs/observed-means-geometric-method.md) |
| Modelled means (the location model, RTMB; `method = rtmb`) | `species_means/<build>/modelled_means.*` | `terminology.md`, `docs/local-range-model.qmd`, `docs/within-species-centring.qmd`, `reports/analysis.qmd` (rendered as `reports/analysis.html`) |
| Forecasting tables (draft, unpublished): a named species' trait value at any climate, a species without trait data from its growth form and climate centre, and a growth form's community statistics at any climate, each with uncertainty | `forecast/<build>/`: `species_climate_responses.*` (per trait × species: the mean at the climate centre, the two within-species slopes and their covariance, the records' climate range), `site_variances.*` (per trait × growth form: what to add for a new site or record), `community_climate_surface.*` (community statistics on a climate grid), `growth_form_responses.*` (per trait × growth form: what forecasting a species without trait data needs) and `species_climate_centres.*` (the 22,632 species that can be forecast that way: climate centre and growth form) | `reports/forecasting.qmd` (how to use them, with checks against the data), `R/process/forecast.R` |
| Combined species means (modelled where available, otherwise observed; internal input) | `species_means/<build>/species_means_combined.*` | `terminology.md` |
| Community statistics, 50 km equal-area grid, from each cell's climate alone (draft, unpublished) | `community/<build>/`: `community_statistics.parquet` (one row per trait × cell × growth form: community mean and SD with 94% intervals, 90% range, sorting factor, extrapolation flags), `sorting_factors.csv`, `geotiff/`, `figures/<trait>_community_statistics`, `community_statistics_metadata.json`, `check/` (against each cell's own GBIF species), `woody_lookup.csv` | `terminology.md`, `docs/local-range-model.qmd`, `R/community/community_climate.R` |

`<build>` is the version of the AusTraits release the products came from: currently `8.0.0` ([10.5281/zenodo.24543275](https://doi.org/10.5281/zenodo.24543275)). Alongside, under `outputs/`: `figures/` (the figure set), `comparisons/` (observed vs modelled; the pooled and locations-pooled models against the location model, `pooled_vs_location_*`; and the refit against the June 2026 baseline, the earlier outputs kept in `docs/baselines/2026-06/`), `checks/` (the medoid checks), and `traits/`, the pipeline's working per-trait files, from the location model (each trait's modelled means as `<trait>_climate_response.csv.gz` with a `_meta.csv` dictionary, the fixed-effect summary per trait and growth form `trait_model_summary.csv`, `data_basis.json`, and per-trait figures: `figures/<trait>` with the species' modelled means against climate, and `figures/<trait>_within_between` with, per growth form, the records, each species' within-species line, and the between-species and average within-species lines for rainfall and temperature); these are not versioned products. The two rendered reports, `analysis.html` and `dataset_construction.html`, are written to `reports/`, next to their sources (git-ignored).

To read one:

```r
means <- arrow::read_parquet("outputs/species_means/8.0.0/modelled_means.parquet")
readr::read_csv("outputs/species_means/8.0.0/modelled_means_dictionary.csv")
```

## Building

You need R 4.6. Steps, from the repository root:

1. Install the pinned packages: `renv::restore()`.
2. Put the inputs in place (all large and git-ignored):

| Input | Where | How to get it |
|---|---|---|
| AusTraits 8.0.0: the database and its observed means | `data/austraits/austraits-8.0.0.rds` and `data/austraits/austraits-8.0.0-observed-means.parquet` (set by `austraits_rds` and `austraits_observed_means` in `_targets.R`) | from the release's Zenodo record, [10.5281/zenodo.24543275](https://doi.org/10.5281/zenodo.24543275); the pipeline checks each file's MD5 against the record. Commands, and how to use a development build instead, in [`data/austraits/README.md`](data/austraits/README.md) |
| WorldClim 2.1 at 2.5′, BIO1 and BIO12 | `data/climate/wc2.1_2.5m_bio_1.tif`, `data/climate/wc2.1_2.5m_bio_12.tif` | [worldclim.org](https://www.worldclim.org) |
| Australia extent raster | `data/australia.tif` | any raster covering 112.9–153.7° E, 43.8–9.1° S; it is used only for its extent, to crop the climate layers |
| Cleaned GBIF occurrences | `data-processed/data_gbif_cleaned_<key>.parquet` (the extract in use is set by `gbif_cleaned_file` in `_targets.R`: GBIF download [10.15468/dl.766erw](https://doi.org/10.15468/dl.766erw), 15 October 2026) | `Rscript R/scripts/build_gbif.R <key> <folder>` fetches a GBIF download by key (no credentials needed), cleans it, keeps only CC0 and CC BY records, aligns names to the pinned APC release (from the provider's name, GBIF's name as fallback) and attaches climate, writing `data-processed/data_gbif_cleaned_<key>.parquet`; downloads used are listed in `gbif_downloads` (`R/data/gbif.R`) |
| Species medoids | `data-processed/climate_centres_<key>.parquet` (set by `medoids_file` in `_targets.R`) | `Rscript R/scripts/build_medoids.R <key> <workers>` after the GBIF step (about 9 minutes for 22,900 taxa on 10 workers) |

Only those two files in `data-processed/` are read; the others there are legacy caches.

3. Build and test:

```r
targets::tar_make()                   # build everything that is out of date
testthat::test_dir("tests/testthat")  # unit tests
```

With the inputs in place, `tar_make()` takes about 40 minutes from an empty targets store, half of it the location model for plant height (which needs about 9 GB of memory). It rebuilds only what changed and renders the three reports at the end. `targets::tar_visnetwork()` draws the dependency graph; `targets::tar_outdated()` lists what a make would rebuild. The `brms` backend is opt-in (`run_brms <- TRUE` in `_targets.R`; diagnostics go to `outputs/traits/brms_diagnostics.csv`, fits are cached in `data-processed/brms_fits/`).

**Pipeline tiers.** The cleaned GBIF occurrences (~1.2 GB) and the 22k-species medoids are tracked as file inputs: `tar_make()` watches their hashes but does not recompute them; rebuild them deliberately with `R/scripts/build_gbif.R` and `R/scripts/build_medoids.R`. Everything else is computed by the pipeline, from the AusTraits release onward: flattening, the name crosswalk, the unified record table, the models, the products, the layers, the figures and the reports.

## Layout

```text
.
├── _targets.R                 # the pipeline (orchestration only; calls R/ functions)
├── renv.lock                  # pinned dependency versions
├── roadmap.md  terminology.md  AGENTS.md
├── reports/                   # the three reports, rendered by the pipeline next to their .qmd:
│                              #   analysis.qmd (models, shifts against June 2026, observed vs modelled),
│                              #   dataset_construction.qmd (data basis, name matching, medoid checks),
│                              #   forecasting.qmd (using the forecasting tables, with checks)
├── docs/                      # method and design documents, the June 2026 baseline,
│                              #   within-species-centring.qmd (why within-species centring),
│                              #   collection-method-notes.md (working notes for #3)
├── R/
│   ├── init.R                 # source() this to load the whole library
│   ├── data/                  # records, climate, crosswalk, medoid lookup, roll-up
│   ├── medoid/                # medoid location and checks
│   ├── model/                 # the location model (RTMB), pooled models, per-trait runs
│   ├── process/               # species centres, exports, comparisons
│   ├── means/                 # observed means (read from the release), reliability, writers
│   ├── community/             # woody lookup, community statistics
│   ├── plot/  panels/  utils/ # figures (save_figure: PNG + PDF), helpers
│   └── scripts/               # run-scripts (not sourced by init.R)
├── tests/testthat/            # tests
├── data/  data-processed/     # inputs, git-ignored except data/austraits/README.md (provenance)
├── _targets/                  # targets store (git-ignored)
└── outputs/                   # products, figures, checks (git-ignored)
```

## Licence

The code and the data products in this repository are licensed under [CC BY 4.0](LICENSE). The input data keep their own licences and citation requirements: AusTraits, AVH, GBIF/ALA and WorldClim.

## Contributing

Issues across the AusTraits family are tracked on one board, [AusTraits #9](https://github.com/orgs/traitecoevo/projects/9), and new issues are auto-added. Please read the [issue & labelling guide](https://github.com/traitecoevo/austraits-meta/blob/master/governance/issue-guide.md) in [`austraits-meta`](https://github.com/traitecoevo/austraits-meta) — the family's cross-package knowledge and governance hub — before filing.

## Acknowledgements

AusTraits is made possible by contributions from our partner organisations — the
[University of New South Wales](https://www.unsw.edu.au/),
[Western Sydney University](https://www.westernsydney.edu.au/),
[Botanic Gardens of Sydney](https://www.botanicgardens.org.au/),
[the University of Melbourne](https://www.unimelb.edu.au/),
the [Atlas of Living Australia](https://www.ala.org.au/), and the Australian Government
[Department of Climate Change, Energy, the Environment and Water](https://www.dcceew.gov.au) — and
from our [advisory board, data contributors, and past partners](https://austraits.org/team/team-partners.html).

AusTraits is a co-investment partnership with the
[Australian Research Data Commons](https://ardc.edu.au/) (ARDC) through the Planet Research Data
Commons ([DOI: 10.3565/nyk4-4r91](https://doi.org/10.3565/nyk4-4r91)). The ARDC is enabled by the
Australian Government's [National Collaborative Research Infrastructure Strategy](https://www.education.gov.au/ncris)
(NCRIS).
