# CLAUDE.md

Guidance for working in this repository.

## What this project does

Estimates the **mean of each plant trait at the centre of a species' climate
distribution** (mean annual precipitation + mean annual temperature), from
AusTraits + AVH data, using a hierarchical within-group-centred ("Mundlak")
model.

Two framings:
1. **National species-trait estimates** for Australia, under the ARDC-funded
   AusTraits programme — a consistent, defensible species mean per trait.
2. **A paper** reporting these species means and giving credence to the
   estimates (uncertainty intervals, within/between decomposition, full
   posterior inference via the BRMS backend).

(An earlier height analysis for the NVIS team was a reportable subproject; that
framing is retired — see `archive/Data_height_AVH.qmd`.)

The method and its statistical motivation (why naive mixed models bias the
within-species slope when group means correlate with predictors) are in
`demo.qmd`; the case for full posterior inference is in
`archive/notes_why_bayesian.md` (rough working notes).
The full roadmap is in **`refactor_plan.md`** — read it before making structural
changes.

## Notebook structure

**Exactly two Quarto notebooks** drive the analysis; both only call functions
from `R/`:

- `dataset_construction.qmd` — climate covariates → GBIF + medoid centres
  (+ reasonableness check) → trait sources (AVH height, AusTraits) → unified
  record-level table.
- `analysis.qmd` — per-trait model fitting → processing → exports → figures +
  write-up.

Supporting files: `demo.qmd` (method illustration — the justification for the
within-group-centred approach), `trait_species_means.qmd` (round-2 simple means
for low-coverage traits — WIP). Superseded/scratch material lives in `archive/`
(the retired NVIS height report, the old medoid exploration, and
`notes_why_bayesian.md` — rough rationale notes for the round-2 Bayesian work).

## R/ library layout

All logic lives in functions under `R/`, loaded via **`source("R/init.R")`**
(sources every function file; deliberately excludes `R/scripts/`).

```
R/
  init.R       # source this once; loads the whole library
  utils/       save_figure.R (png+pdf), parallel_setup.R (multisession)
  data/        build_dataset.R, extract_climate.R, assign_state.R,
               load_data_height_avh.R
  medoid/      locate_medoid_2d.R, check_medoid.R
  model/       model_formula.R, fit_trait_climate_model.R (tmb | brms),
               run_trait_models.R
  process/     process_fit.R, species_centres.R, export_centres.R
  plot/        plot_fitted.R, plot_trait_centres.R
  means/       simple_means.R (round-2 fallback), species_means.R
  panels/      species_panel_functions.R
  scripts/     generate_species_panels.R, benchmark_species_panels.R,
               combine_species_panels.R   # RUN-SCRIPTS — execute on source,
                                           # NOT loaded by init.R
```

`R/scripts/*` are runnable from the project root (`Rscript R/scripts/x.R`); they
do real work when sourced, so they must never be sourced by `init.R`.

## Conventions

- **All logic in `R/` functions; notebooks only call, parameterize, narrate.**
- **Cache expensive steps with `if` guards**, one `update <- FALSE` toggle per
  notebook:
  ```r
  if (update || !file.exists(out)) { result <- expensive(); write(result, out) }
  result <- read(out)
  ```
  Cached artefacts live in `data-processed/` (parquet/rds), outputs in `outputs/`.
- **Parallelism: `furrr` + `setup_parallel()` → `plan(multisession)`.** Never
  `multicore` — forking is unstable/disabled on macOS/RStudio. `setup_parallel()`
  reads true core count via `parallelly::availableCores(methods = "system")` to
  work around a low `mc.cores` in `.Rprofile`.
- **Figures: save as both PNG and PDF via `save_figure(plot, path, w, h)`** (PNG
  for embedding/preview, PDF as vector deliverable).

## Modelling decisions (locked — see refactor_plan.md §5)

- **Centring on the medoid (Option B):** `x_c = obs − medoid`, `x_mn = medoid`
  per climate var; consistent across traits.
- **Per-trait models** now (TMB default, BRMS for the deliverable); joint
  multivariate / gap-filling deferred to round 2.
- **Two backends** behind one interface: `fit_trait_climate_model(backend =
  "tmb" | "brms")`. Build BRMS validated against TMB on a subset, then scale
  (`cmdstanr`, non-centered parameterization). Consult the `bayesian-workflow`
  skill before writing Stan/brms code.
- **No medoid → dropped** from the climate-response models (routed to round-2
  means). Medoid quality is gated by a reasonableness check.
- **AusTraits records:** field, `value_type ∈ {raw, mean, median}`.
- **`dataset_id`:** one per AVH collector (+ AusTraits study id) — carried for
  the future methodology-factor extension.

## Running BRMS at scale

TMB is the default (seconds per trait). BRMS gives the full posterior but is slow,
so it is opt-in:

```r
brms_traits <- setdiff(unique(data_unified$trait_name), "plant_height")  # height stays on TMB
fits <- run_trait_models(data_unified, traits = brms_traits, backend = "brms",
                         fits_dir = "data-processed/brms_fits",  # on-disk fit cache
                         threads = 4, iter = 2000, warmup = 1000,
                         refresh = 50, silent = 1)  # keep warmup progress visible
summarise_trait_diagnostics(fits)   # divergences, max Rhat, min ESS per trait
```

- **Height is TMB-only, excluded from BRMS.** In the unified dataset height is
  ~306k records / ~10,200 species (n>5) → ~30,000 random-effect parameters. NUTS
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

- `data/` — inputs (AVH CSV, GBIF/ALA, WorldClim `data/climate/`, `australia.tif`,
  AusTraits `data/data-ausplants/`). `data/` and `downloads/` are symlinks into
  the sibling `Cornsters-australian_annuals` project.
- `data-processed/` — cached intermediates: `data_gbif_cleaned.parquet`,
  `climate_centres.parquet` (medoids), `gbif_avh_name_crosswalk.parquet`, the
  source caches `data_height_climate.parquet` + `data_traits_climate.parquet`,
  the unified `data_traits_climate_unified.parquet`, and `brms_fits/` (BRMS cache).
- `outputs/traits/` — per-trait `<trait>_climate_response.csv.gz` (+ `_meta.csv`),
  `species_trait_means_combined.csv.gz`, `trait_model_summary.csv`, figures.

## Gotchas

- **Taxonomy across releases:** AVH, AusTraits and GBIF were aligned to
  different APC releases. Joining medoids/occurrences to trait records needs a
  crosswalk (align all to current APC with `APCalign`); `build_medoid_lookup()`
  bridges AVH/AusTraits names to GBIF medoids via the cached crosswalk. See also
  the crosswalk build in `archive/Data_height_AVH.qmd`, and open issue #1 (APC
  versions out of sync).
- **Medoid memory:** the 2D distance matrix in `locate_medoid_2d()` is O(n²);
  `process_taxon_medoid()` subsamples species with >20k occurrences.
- **`.Rprofile` may set `mc.cores` low** — use `setup_parallel()`.
- Rerunning the full pipeline is expensive but acceptable; outputs are
  regenerated rather than preserved bit-for-bit across refactors.
