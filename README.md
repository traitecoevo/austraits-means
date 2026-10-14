# austraits-means

Estimating **national means of plant traits at the centre of each species'
climate distribution** for Australia, from AusTraits + AVH herbarium records and
WorldClim covariates. Part of the ARDC-funded AusTraits programme; also the basis
for a paper reporting these species means.

A species' trait value varies across its range, and trait data are sampled
unevenly across that range. We therefore estimate, per trait, the species mean at
the **centre of the species' climate distribution** (a 2D precipitation–temperature
medoid), using a hierarchical within-group-centred model. Two interchangeable
fitting backends are available: `tmb` (fast) and `brms` (full posterior).

See `AGENTS.md` for conventions, `roadmap.md` for the roadmap and decisions log, and
`demo.qmd` for the statistical motivation.

## Workflow

Orchestration is a [`targets`](https://docs.ropensci.org/targets/) pipeline
(`_targets.R`); all logic lives in functions under `R/`. Dependencies are pinned
with [`renv`](https://rstudio.github.io/renv/) (`renv::restore()` to reproduce
the library).

```r
renv::restore()         # first time: install the pinned package versions
targets::tar_make()     # build everything that is out of date
targets::tar_visnetwork()  # view the dependency graph
targets::tar_outdated()    # preview what a make would rebuild
```

`tar_make()` rebuilds only what changed (content-hash invalidation — no manual
`update` toggles) and renders the two reports at the end:

- `dataset_construction.qmd` — provenance of the data-building inputs + medoid
  reasonableness checks.
- `analysis.qmd` — the write-up: per-trait fitted effects, across-species
  turnover figures, and the round-2 combined best-estimate product.

Both notebooks are **reporting layers**: they `tar_load()` pipeline results and
narrate them, they do not compute.

**Pipeline tiers.** The heavy, rarely rebuilt artefacts (cleaned GBIF ~1.2 GB,
the 22k-species medoids, the AVH↔GBIF crosswalk, the AusTraits flatten) are
tracked as **file inputs** — `tar_make()` watches their hashes but does not
recompute them. Regenerate them deliberately with the builders in
`R/data/build_sources.R` (`build_gbif_climate()`, `build_austraits_climate()`),
`load_data()` (height), and the producer notebooks at the root
(`climate_centres.qmd` for medoids, `Data_height_AVH.qmd` for the crosswalk). Everything from the **unified
record-level dataset** forward is computed and rebuilt incrementally.

The full posterior (`brms`) backend is opt-in: set `run_brms <- TRUE` in
`_targets.R` (height stays on `tmb` — too large for NUTS).

## Layout

```text
.
├── _targets.R                   # the pipeline (orchestration; calls R/ functions)
├── renv.lock                    # pinned dependency versions
├── dataset_construction.qmd     # report: data provenance + medoid checks
├── analysis.qmd                 # report: per-trait models -> estimates + figures
├── trait_species_means.qmd      # round-2 simple means (low-coverage traits, WIP)
├── demo.qmd                     # within-group-centring method illustration
├── climate_centres.qmd          # producer: species climate medoids (run deliberately)
├── Data_height_AVH.qmd          # producer: AVH<->GBIF crosswalk (+ retired NVIS height report)
├── AGENTS.md  terminology.md  roadmap.md
├── R/
│   ├── init.R                   # source() this to load the whole library
│   ├── data/                    # build_dataset, build_sources, extract_climate, assign_state, load_data_height_avh
│   ├── medoid/                  # locate_medoid_2d, check_medoid
│   ├── model/                   # model_formula, fit_trait_climate_model, run_trait_models
│   ├── process/                 # process_fit, species_centres, export_centres
│   ├── plot/                    # plot_fitted, plot_trait_centres
│   ├── panels/                  # species diagnostic panels
│   ├── means/                   # simple_means (round-2 fallback) + observed_means_portal
│   ├── utils/                   # save_figure (png+pdf), parallel_setup (multisession)
│   └── scripts/                 # run-scripts (NOT sourced by init.R)
├── _targets/                    # targets store (git-ignored; rebuilt by tar_make)
├── data/                        # inputs (git-ignored): AVH, GBIF, WorldClim, AusTraits
├── data-processed/              # cached intermediates + seeded inputs (git-ignored)
└── outputs/                     # figures + per-trait response CSVs (git-ignored)
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
