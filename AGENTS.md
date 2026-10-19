# AGENTS.md

Guidance for working in this repository, for agents and contributors alike (`.claude/CLAUDE.md` defers here). Family-wide conventions are at the end: [AusTraits family — cross-package context](#austraits-family--cross-package-context).

## What this project does

Builds species-level and community-level trait means for the Australian flora from AusTraits, with WorldClim climate and GBIF occurrences. Three products (`terminology.md`, `roadmap.md`):

1. **Observed means**, computed directly from the records: the *portal method* (`method = portal_arithmetic`), which reproduces exactly the species means the AusTraits data portal shows for release 7.0.0, and its successor the *geometric method* (`method = geometric`), on the austraits.build development build.
2. **Modelled means** (`method = tmb` or `brms`): the mean of each trait at the centre of a species' climate distribution (the 2D precipitation–temperature medoid of its occurrences), from a hierarchical within-group-centred ("Mundlak") model.
3. **Community layers**: the community mean and community variance of the species means present in each spatial unit, woody and non-woody, first on provisional choices (`docs/community-layers-design.md`).

Two framings:
1. **National species-trait estimates** for Australia, under the ARDC-funded
   AusTraits programme — a consistent, defensible species mean per trait.
2. **A paper** reporting these species means and giving credence to the
   estimates (uncertainty intervals, within/between decomposition, full
   posterior inference via the BRMS backend).

(An earlier height analysis for the NVIS team was a reportable subproject; that framing is retired. The report survives as `Data_height_AVH.qmd`, which reads the old AVH caches and no longer feeds the pipeline.)

The method and its statistical motivation (why naive mixed models bias the within-species slope when group means correlate with predictors) are in `demo.qmd`; the case for full posterior inference was in `archive/notes_why_bayesian.md` (rough working notes, now in git history only).
The roadmap, milestones and decisions log are in **`roadmap.md`** — read it before making structural changes.

**Terminology is fixed in `terminology.md`**: the two published estimates are the *observed mean* and the *modelled mean* (both geometric from the geometric method onward; collectively "species means"), and community-level summaries are the *community mean* and *community variance* (split woody / non-woody). Use those terms in code, exports, docs and the paper; never "average". Exports carry `estimate_type` and `method` and name values `value`, `log_value`, `value_lower`, `value_upper`.

## Workflow: `targets` pipeline + two reporting notebooks

The analysis is orchestrated by a **`targets` pipeline** (`_targets.R`); `targets::tar_make()` builds the DAG and only rebuilds what changed. All logic lives in `R/` functions — `_targets.R` and the notebooks only call them. Tests: `testthat::test_dir("tests/testthat")` from the root.

The two Quarto notebooks are **reporting layers** (rendered by `tar_quarto` targets): they `tar_load()` pipeline results and narrate them, they do not compute.

- `dataset_construction.qmd` — the data basis, name matching and medoid reasonableness checks.
- `analysis.qmd` — fitted fixed effects, the refit against the June 2026 baseline, observed vs modelled means, across-species turnover, the products.

Sections of `_targets.R`, in order: tracked inputs; the development-build records, crosswalk, medoid lookup and `unified`; the per-trait TMB fits and their exports (modelled means); the portal-method product for 7.0.0, its reproduction check and its Zenodo package; the geometric observed means; the woody lookup, combined species means and community layers (variants via `tar_map`) and the model surface; the comparisons; the figure set; the reports.

**Pipeline tiers** (see the header comment in `_targets.R`):
- *Tier B — tracked file inputs.* The cleaned GBIF occurrences (`data_gbif_cleaned.parquet`, ~1.2GB) and the 22k-species O(n²) medoids (`climate_centres.parquet`) are tracked by hash but **not** recomputed by `tar_make()`. Regenerate them deliberately with `R/scripts/build_gbif.R` (builders in `R/data/gbif.R`: download by key, cleaning ported from cornsters_australian_annuals, CC0/CC BY records only, APC alignment, climate) and `R/scripts/build_medoids.R` (`build_climate_medoids()`, parallel over taxa). The AusTraits rds files (`data/austraits/`) are file inputs too; how the development build is made is in `data-raw/austraits/README.md`.
- *Tier A — computed.* Everything else: flattening the AusTraits builds, the name crosswalk (APC release pinned as `apc_version`), the unified record table, the fits, the products, the community layers, the figures and the rendered reports.

BRMS is opt-in: set `run_brms <- TRUE` in `_targets.R` to add the
`brms_fits`/`brms_diag` targets (height excluded — TMB only).

Supporting files: `demo.qmd` (method illustration — the justification for the within-group-centred approach), `trait_species_means.qmd` (legacy notebook; predates the products and calls functions since retired).

### Producer notebooks

None remain. The species medoids, once built by `climate_centres.qmd`, come from `build_climate_medoids()` (`R/data/build_sources.R`) via `R/scripts/build_medoids.R` since #33; the notebook is kept as an exploratory record. `Data_height_AVH.qmd` is the retired NVIS height report and no longer feeds the pipeline. (The rest of the old `archive/` — scratch notes and `notes_why_bayesian.md` — was removed in `eecdbf2` and survives only in git history.)

## R/ library layout

All logic lives in functions under `R/`, loaded via **`source("R/init.R")`**
(sources every function file; deliberately excludes `R/scripts/`).

```
R/
  init.R       # source this once; loads the whole library
  utils/       save_figure.R (png+pdf), parallel_setup.R (multisession)
  data/        build_dataset.R (records, value-type rules, roll-up, medoid
               lookup, name matching), build_sources.R (crosswalk, AusTraits
               flatten, GBIF builder, data basis), extract_climate.R,
               assign_state.R
  medoid/      locate_medoid_2d.R, check_medoid.R
  model/       model_formula.R, fit_trait_climate_model.R (tmb | brms),
               run_trait_models.R
  process/     process_fit.R, species_centres.R, export_centres.R (modelled
               means product), compare_baseline.R, compare_means.R
  plot/        plot_fitted.R, plot_trait_centres.R, figure_set.R
  means/       observed_means_portal.R, observed_means_geometric.R,
               reliability.R, species_means_io.R, compare_portal.R
  community/   woody_lookup.R, community_layers.R, community_modelled.R
  release/     zenodo.R (draft deposits; never publishes)
  panels/      species_panel_functions.R
  scripts/     generate_species_panels.R, benchmark_species_panels.R,
               combine_species_panels.R, create_zenodo_draft.R
                                           # RUN-SCRIPTS — execute on source,
                                           # NOT loaded by init.R
```

`R/scripts/*` are runnable from the project root (`Rscript R/scripts/x.R`); they
do real work when sourced, so they must never be sourced by `init.R`.

## Conventions

- **All logic in `R/` functions; `_targets.R` + notebooks only call, narrate.**
- **Caching is handled by `targets`** (content-hash invalidation), not manual
  `if (update || !file.exists)` guards — those are retired. Add a step by adding
  a `tar_target()`; `tar_make()` rebuilds only the affected targets and their
  descendants (`tar_outdated()` previews this; `tar_visnetwork()` draws the DAG).
  The targets store is `_targets/` (git-ignored); seeded inputs still live in
  `data-processed/`, deliverables in `outputs/`.
- **`renv`** locks dependencies (`renv.lock`); `renv::snapshot()` after adding a
  package, `renv::restore()` to reproduce the library.
- **Parallelism: `furrr` + `setup_parallel()` → `plan(multisession)`.** Never
  `multicore` — forking is unstable/disabled on macOS/RStudio. `setup_parallel()`
  reads true core count via `parallelly::availableCores(methods = "system")` to
  work around a low `mc.cores` in `.Rprofile`.
- **Figures: save as both PNG and PDF via `save_figure(plot, path, w, h)`** (PNG
  for embedding/preview, PDF as vector deliverable).

## Dependencies (`renv`)

Package versions are pinned with `renv` so the pipeline is reproducible.

- **Activation is automatic.** `.Rprofile` sources `renv/activate.R`, so every R
  session in this project uses the project library at
  `renv/library/…` (not the system library). `.Rprofile`, `renv.lock`, and
  `renv/activate.R` are committed; `renv/library/` is git-ignored.
- **Reproduce the library:** `renv::restore()` (installs the locked versions).
- **After adding/using a new package:** install it into the project
  (`renv::install("pkg")`), then `renv::snapshot()` to record it in `renv.lock`.
  Check `renv::status()` if things feel out of sync.
- **GitHub-sourced packages:** `austraits` and `APCalign` are pinned to their
  `traitecoevo/…` GitHub remotes in the lockfile (not CRAN) — `restore()` fetches
  them from GitHub.
- **Implicit deps gotcha:** `renv` discovers dependencies by scanning code for
  `library()` / `pkg::fun` / `requireNamespace()`. Packages used *implicitly*
  (here `quarto`, needed by `tar_quarto`, and `crew` for the optional parallel
  controller) are not auto-detected, so `_targets.R` references them via
  `requireNamespace(...)` guards purely so `renv` tracks them. If you add another
  implicitly-used package, add a similar reference or it will be dropped from the
  lockfile on the next `snapshot()`.

## Modelling decisions (locked — see the decisions log in roadmap.md)

- **Centring on the medoid (Option B):** `x_c = obs − medoid`, `x_mn = medoid`
  per climate var; consistent across traits.
- **Per-trait models** now (TMB default; BRMS postponed to the paper
  milestone); joint multivariate / gap-filling deferred.
- **Two backends** behind one interface: `fit_trait_climate_model(backend =
  "tmb" | "brms")`. Consult the `bayesian-workflow` skill before writing
  Stan/brms code.
- **No medoid → dropped** from the climate-response models; such species still
  get an observed mean.
- **Records:** field records with `value_type ∈ {raw, mean, median}`, except
  where `dataset_value_type_rules()` says otherwise: the AVH height dataset
  `AVH_2026` carries one label bound per specimen and enters by its minimum
  and maximum records.
- **Taxa:** subspecies, varieties and forms are rolled up to species
  (`roll_up_to_species()`), for the modelled and geometric observed means;
  the portal-method product is not rolled up (one of its kept quirks).
- **`dataset_id`:** height records are all `AVH_2026` (the per-collector
  grouping of the old AVH CSV is not carried by the build; see #3).

## Running BRMS at scale

TMB is the default (seconds per trait). BRMS gives the full posterior but is slow,
so it is opt-in: set `run_brms <- TRUE` in `_targets.R` to add the
`brms_fits` / `brms_diag` targets, then `tar_make()`. The targets wrap
`run_trait_models(..., backend = "brms")` with the settings below (height
excluded; on-disk fit cache in `data-processed/brms_fits/`; `threads ≈ cores/4`;
`iter = 2000, warmup = 1000, refresh = 250, silent = 1`), and
`summarise_trait_diagnostics()` produces the per-trait divergence / Rhat / ESS
table. The same call works interactively for ad-hoc fits:

```r
brms_traits <- setdiff(unique(data_unified$trait_name), "plant_height")  # height stays on TMB
fits <- run_trait_models(data_unified, traits = brms_traits, backend = "brms",
                         fits_dir = "data-processed/brms_fits",  # on-disk fit cache
                         threads = 4, iter = 2000, warmup = 1000,
                         refresh = 50, silent = 1)  # keep warmup progress visible
summarise_trait_diagnostics(fits)   # divergences, max Rhat, min ESS per trait
```

- **Height is TMB-only, excluded from BRMS.** In the unified dataset height is
  ~332k records / ~10,300 species (n>5) → ~31,000 random-effect parameters. NUTS
  there runs many hours: the per-iteration cost stays flat (~11s) across warmup, so
  mass-matrix adaptation never rescues it, with real divergence/memory risk. Not
  worth the full posterior for one trait — use TMB's asymptotic intervals. (The old
  "~114s for 712 species" benchmark predates source unification; it is obsolete.)
- **Threading** (`threads = n`, reduce_sum) speeds the large fits; 4 chains run in
  parallel so total cores ≈ `4 * threads` — size `threads ≈ cores / 4`. Note
  reduce_sum scales sublinearly (~60% efficiency) because the gradient still
  carries all RE parameters serially; it helps but is not 1:1 with thread count.
- **Visibility:** keep `refresh = 50, silent = 1`. With `refresh = 0`, warmup emits
  nothing and nothing is checkpointed until completion, so a long fit gives no read
  on progress. (To check a live fit externally: `ps … | grep method=sample` shows
  per-chain args; post-warmup draw counts appear as non-comment rows in the cmdstan
  CSVs under the session tempdir, but warmup draws are not written.)
- **Caching** (`fits_dir`) writes `<trait>_brms.rds`; reused unless model/data
  change (reload ≈ 1s vs minutes to fit). Only written on completion — an
  interrupted fit checkpoints nothing.
- **Diagnostics** are captured per fit (`fit$diagnostics`) and flagged to console;
  always check them. Benchmark: ~44–53s for 142–309-species traits; **larger traits
  need more iterations** — `leaf_area` did not converge at `iter = 1500` (Rhat
  1.115). Refit flagged traits with `iter = 4000`.
- One-time `cmdstanr::rebuild_cmdstan()` clears stale precompiled-header (PCH)
  warnings the first time threading turns on.
- `data-processed/brms_fits/` is git-ignored (fits are large).

## Climate layers

WorldClim 2.5 arc-min, cropped to Australia (`load_climate_rasters()`):
- BIO_1  → `average_temperature` (°C)
- BIO_12 → `annual_precipitation` (mm/yr); modelled on the log scale.

## Data locations (git-ignored; large)

- `data/` — a directory of inputs (only `data/data-ausplants` is a symlink, into a shared drive):
  `austraits/` (the AusTraits 7.0.0 release rds and the development-build export), `portal/` (the portal's deployed species means, for the reproduction check), `climate/` (WorldClim), `australia.tif`; the older AVH height CSVs are no longer read. Raw GBIF downloads live outside the repository; only the cleaned extract is used.
- `data-processed/` — seeded pipeline inputs tracked by `targets` as file inputs: `data_gbif_cleaned.parquet`, `climate_centres.parquet` (medoids), plus `brms_fits/` (BRMS fit cache). `data_height_climate.parquet`, `data_traits_climate.parquet`, `gbif_avh_name_crosswalk.parquet` and `data_traits_climate_unified.parquet` are no longer read by the pipeline (left in place; don't delete without asking).
- `outputs/` — `species_means/<version>/` (the products), `community/<version>/` (layers), `traits/` (per-trait modelled means, fixed-effect summary, `data_basis.json`), `comparisons/`, `checks/`, `figures/`, `zenodo/` (deposit packages).

## Gotchas

- **Taxonomy across releases:** the AusTraits build and the GBIF extract were
  aligned to different APC releases. `build_name_crosswalk()` aligns both name
  sets to one pinned release (`apc_version` in `_targets.R`) and
  `build_medoid_lookup()` uses it to bridge trait names to GBIF medoids; open
  issue #1 tracks aligning the GBIF side itself (release candidate).
- **Global `*.gz` git ignore:** the maintainer's global gitignore ignores
  `*.gz`, so a tracked `.csv.gz` (such as `docs/baselines/2026-06/species_centres.csv.gz`)
  needs `git add -f`.
- **Medoid memory:** the 2D distance matrix in `locate_medoid_2d()` is O(n²);
  `process_taxon_medoid()` subsamples species with >20k occurrences.
- **`.Rprofile` may set `mc.cores` low** — use `setup_parallel()`.
- Rerunning the full pipeline is expensive but acceptable; outputs are
  regenerated rather than preserved bit-for-bit across refactors.

---

## AusTraits family — cross-package context

`austraits-means` is part of the **AusTraits family** (a subset of the [`traitecoevo`](https://github.com/traitecoevo) org) — here, species-level trait means (observed and climate-centred modelled) built downstream of the released AusTraits database, with taxonomy aligned by `APCalign`. Family-wide concerns are documented centrally in **[austraits-meta](https://github.com/traitecoevo/austraits-meta)** — don't restate them here, read them there:

- **Start with [`AGENTS.md`](https://github.com/traitecoevo/austraits-meta/blob/master/AGENTS.md)** — pipeline order, who owns what, dependency direction, source-of-truth rules, cross-boundary artifacts, gotchas.
- **[`dependencies.yml`](https://github.com/traitecoevo/austraits-meta/blob/master/dependencies.yml)** — machine-readable package graph + cross-boundary artifacts.
- **[`governance/`](https://github.com/traitecoevo/austraits-meta/tree/master/governance)** — label taxonomy, board #9 conventions, release playbooks, triage.

Source-of-truth rules that bite here: trait definitions and units come from **APD**, and taxonomy from **APCalign** — consume them, never redefine them in this repo.

**Filing issues:** the whole family is tracked on one board, [AusTraits #9](https://github.com/orgs/traitecoevo/projects/9) (new issues auto-add to it). Follow the [issue & labelling guide](https://github.com/traitecoevo/austraits-meta/blob/master/governance/issue-guide.md): pick one work-type label (`bug` / `task` / `epic`); Status and Priority are set on the board, not as labels.

**No PRs in this repo — push straight to `master`** (maintainer decision, 2026-10-14). This is a deliberate exception to the family's feature-branch + PR discipline; don't open PRs here unless asked. The family's commit-message shape still applies to the commits themselves: subject ≤50 characters as typed, body ≤10 lines saying why and what changed. See [`commit-messages.md`](https://github.com/traitecoevo/austraits-meta/blob/master/governance/commit-messages.md).

> austraits-meta is hand-maintained prose — a map, not ground truth. Verify specifics against the actual repos.
