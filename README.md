# austraits-means

Species-level and community-level trait means for the Australian flora, built from [AusTraits](https://austraits.org) with WorldClim climate and GBIF occurrences. Part of the AusTraits programme.

The repository produces three data products, with the terms fixed in [`terminology.md`](terminology.md):

1. **Observed means**: per species × trait, computed directly from the records. Two methods: the *portal method*, which reproduces exactly the species means the AusTraits data portal shows for release 7.0.0 (an arithmetic mean), and its successor, the *geometric method*.
2. **Modelled means**: per species × trait, the geometric mean at the centre of the species' climate distribution, from a hierarchical model in which each species' climate is centred on its medoid (the most central point of its GBIF occurrences in precipitation–temperature space; `demo.qmd` explains why). Fitted with TMB (`glmmTMB`, maximum likelihood); a Bayesian `brms` backend is optional.
3. **Community layers** (first draft, on provisional choices): the community mean and community variance of the species means present in each 0.25° cell, for woody and non-woody plants, plus a model-based surface.

The *combined species means* table is not a fourth product: it merges products 1 and 2 (the modelled mean where a species has one, otherwise its geometric observed mean) and is the input to the community layers.

[`roadmap.md`](roadmap.md) records the milestones, the decisions and their reasons, and the invariants; [`AGENTS.md`](AGENTS.md) the conventions for working in the code.

## The products

Every species means table carries `estimate_type` (`observed_mean` or `modelled_mean`) and `method` (`portal_arithmetic`, `geometric`, `tmb` or `brms`), and comes as four files: `.parquet`, `.csv.gz`, `_dictionary.csv` (every column defined) and `_metadata.json` (provenance: source build, commit, data basis). Each product carries the version tag of the AusTraits release or build it derives from.

| Product | Files (under `outputs/`) | Documentation |
|---|---|---|
| Observed means, portal method, for AusTraits 7.0.0 | `species_means/7.0.0/observed_means_portal_arithmetic.*` | [`docs/observed-means-portal-method.md`](docs/observed-means-portal-method.md) |
| Observed means, geometric method | `species_means/<build>/observed_means_geometric.*` | [`docs/observed-means-geometric-method.md`](docs/observed-means-geometric-method.md) |
| Modelled means (TMB) | `species_means/<build>/modelled_means_tmb.*` | `terminology.md`, `demo.qmd`, `analysis.qmd` (rendered as `analysis.html`) |
| Combined species means (modelled where available, otherwise observed) | `species_means/<build>/species_means_combined.*` | `terminology.md` |
| Community layers (10 km equal-area grid) | `community/<build>/<product>/`: `community_layers.parquet` (one row per trait × cell × group; `woodiness_group` is `woody`, `non_woody` or `all`), `geotiff/`, `figures/`, `community_layers_metadata.json`; `woody_lookup.csv` alongside. Products: `observed_local` (each species' records in the cell, then the species), `observed_presence` (GBIF species present × observed means), `modelled` (the trait model's net average for the cell, all species) | [`docs/community-layers-design.md`](docs/community-layers-design.md) |

`<build>` is the version string of the austraits.build development build the products came from: currently `7.0.0.9000`, built from austraits.build commit `d535de9`. Alongside, under `outputs/`: `figures/` (the figure set), `comparisons/` (observed vs modelled; portal vs geometric observed means, which differ in data version as well as method, 7.0.0 against the development build; and the refit against the June 2026 baseline, the earlier outputs kept in `docs/baselines/2026-06/`), `checks/` (the portal reproduction check and the medoid checks), and `traits/`, the pipeline's working per-trait files (each trait's modelled means as `<trait>_climate_response.csv.gz` with a `_meta.csv` dictionary, the fixed-effect summary `trait_model_summary.csv`, and `data_basis.json`); these are not versioned products. The two rendered reports, `analysis.html` and `dataset_construction.html`, are written to the repository root.

To read one:

```r
means <- arrow::read_parquet("outputs/species_means/7.0.0.9000/modelled_means_tmb.parquet")
readr::read_csv("outputs/species_means/7.0.0.9000/modelled_means_tmb_dictionary.csv")
```

## Building

You need R 4.6. Steps, from the repository root:

1. Install the pinned packages: `renv::restore()`.
2. Put the inputs in place (all large and git-ignored):

| Input | Where | How to get it |
|---|---|---|
| AusTraits 7.0.0 release | `data/austraits/austraits-7.0.0.rds` | `austraits::load_austraits(version = "7.0.0", path = "data/austraits")`; the pipeline checks its MD5 |
| austraits.build development build | `data/austraits/austraits-7.0.0.9000-d535de9.rds` (the file name is set by `austraits_dev_rds` in `_targets.R`) | build the austraits.build repository at that commit and copy its export; steps in [`data-raw/austraits/README.md`](data-raw/austraits/README.md) |
| WorldClim 2.1 at 2.5′, BIO1 and BIO12 | `data/climate/wc2.1_2.5m_bio_1.tif`, `data/climate/wc2.1_2.5m_bio_12.tif` | [worldclim.org](https://www.worldclim.org) |
| Australia extent raster | `data/australia.tif` | any raster covering 112.9–153.7° E, 43.8–9.1° S; it is used only for its extent, to crop the climate layers |
| Cleaned GBIF occurrences | `data-processed/data_gbif_cleaned_<key>.parquet` (the extract in use is set by `gbif_cleaned_file` in `_targets.R`: GBIF download [10.15468/dl.766erw](https://doi.org/10.15468/dl.766erw), 15 October 2026) | `Rscript R/scripts/build_gbif.R <key> <folder>` fetches a GBIF download by key (no credentials needed), cleans it, keeps only CC0 and CC BY records, aligns names to the pinned APC release (from the provider's name, GBIF's name as fallback) and attaches climate, writing `data-processed/data_gbif_cleaned_<key>.parquet`; downloads used are listed in `gbif_downloads` (`R/data/gbif.R`) |
| Species medoids | `data-processed/climate_centres_<key>.parquet` (set by `medoids_file` in `_targets.R`) | `Rscript R/scripts/build_medoids.R <key> <workers>` after the GBIF step (about 9 minutes for 22,900 taxa on 10 workers) |
| The portal's deployed species means (optional) | `data/portal/austraits-7.0.0-full/austraits-species-averages.parquet` | a copy of the file the AusTraits data portal built for 7.0.0 (`inst/extdata/austraits/austraits-7.0.0-full/` in the portal's own data build; see [`data-raw/portal/README.md`](data-raw/portal/README.md)); used only by the reproduction check; without it the check's targets report an error and the rest of the pipeline still builds |

Only those two files in `data-processed/` are read; the others there are legacy caches.

3. Build and test:

```r
targets::tar_make()                   # build everything that is out of date
testthat::test_dir("tests/testthat")  # tests of the portal method, and the reproduction result
```

With the inputs in place, `tar_make()` takes about 15 minutes from an empty targets store. It rebuilds only what changed and renders the two reports at the end. The portal reproduction check is part of the pipeline (`portal_reproduction` target, written to `outputs/checks/portal_reproduction_*.csv`); the test reads that target's result and is skipped if the target has never been built. `targets::tar_visnetwork()` draws the dependency graph; `targets::tar_outdated()` lists what a make would rebuild. The `brms` backend is opt-in (`run_brms <- TRUE` in `_targets.R`; diagnostics go to `outputs/traits/brms_diagnostics.csv`, fits are cached in `data-processed/brms_fits/`). A product's Zenodo deposit is uploaded as a draft, never published, by `R/scripts/create_zenodo_draft.R`.

**Pipeline tiers.** The cleaned GBIF occurrences (~1.2 GB) and the 22k-species medoids are tracked as file inputs: `tar_make()` watches their hashes but does not recompute them; rebuild them deliberately with `R/scripts/build_gbif.R` and `R/scripts/build_medoids.R`. Everything else is computed by the pipeline, from the AusTraits builds onward: flattening, the name crosswalk, the unified record table, the models, the products, the layers, the figures and the reports.

## Layout

```text
.
├── _targets.R                 # the pipeline (orchestration only; calls R/ functions)
├── renv.lock                  # pinned dependency versions
├── analysis.qmd               # report: models, the refit against June 2026, observed vs modelled
├── dataset_construction.qmd   # report: data basis, name matching, medoid checks
├── demo.qmd                   # why within-species centring
├── climate_centres.qmd        # exploratory medoid notebook; replaced by R/scripts/build_medoids.R
├── Data_height_AVH.qmd        # retired NVIS height report (reads the old AVH caches)
├── trait_species_means.qmd    # legacy notebook; predates the products and uses retired functions
├── methods.md                 # working notes on collection-method factors (#3)
├── roadmap.md  terminology.md  AGENTS.md
├── docs/                      # method and design documents; the June 2026 baseline
├── data-raw/                  # small pinned inputs with provenance (portal core-trait list)
├── R/
│   ├── init.R                 # source() this to load the whole library
│   ├── data/                  # records, climate, crosswalk, medoid lookup, roll-up
│   ├── medoid/                # medoid location and checks
│   ├── model/                 # model formula, fitting (tmb | brms), per-trait runs
│   ├── process/               # species centres, exports, comparisons
│   ├── means/                 # observed means (portal, geometric), reliability, writers
│   ├── community/             # woody lookup, community layers, model surface
│   ├── release/               # Zenodo deposit
│   ├── plot/  panels/  utils/ # figures (save_figure: PNG + PDF), helpers
│   └── scripts/               # run-scripts (not sourced by init.R)
├── tests/testthat/            # tests
├── data/  data-processed/     # inputs (git-ignored)
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
