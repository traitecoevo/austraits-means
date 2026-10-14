# AGENTS.md

Guidance for working in this repository, for agents and contributors alike (`.claude/CLAUDE.md` defers here). Family-wide conventions are at the end: [AusTraits family — cross-package context](#austraits-family--cross-package-context).

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

(An earlier height analysis for the NVIS team was a reportable subproject; that framing is retired. The report was `archive/Data_height_AVH.qmd`, removed from the tree in `eecdbf2` — see [Retired notebooks](#retired-notebooks).)

The method and its statistical motivation (why naive mixed models bias the within-species slope when group means correlate with predictors) are in `demo.qmd`; the case for full posterior inference was in `archive/notes_why_bayesian.md` (rough working notes, now in git history only).
The full roadmap is in **`refactor_plan.md`** — read it before making structural changes.

**Terminology is fixed in `terminology.md`**: the two published estimates are the *observed mean* and the *modelled mean* (both geometric; collectively "species means"), and community-level summaries are the *community mean* and *community variance* (split woody / non-woody). Use those terms in code, exports, docs and the paper; never "average".

## Workflow: `targets` pipeline + two reporting notebooks

The analysis is orchestrated by a **`targets` pipeline** (`_targets.R`);
`targets::tar_make()` builds the DAG and only rebuilds what changed. All logic
still lives in `R/` functions — `_targets.R` and the notebooks only call them.

The two Quarto notebooks are now **reporting layers** (rendered by `tar_quarto`
targets): they `tar_load()` pipeline results and narrate them, they do not
compute.

- `dataset_construction.qmd` — provenance of the data-construction inputs +
  medoid reasonableness checks (loads `unified`, `medoid_gbif_check`).
- `analysis.qmd` — the write-up: fitted fixed effects, across-species turnover
  figures, the round-2 combined product (loads `trait_summary`, `tmb_fits`,
  `combined`).

**Pipeline tiers** (see the header comment in `_targets.R`):
- *Tier B — tracked file inputs.* The heavy, rarely rebuilt artefacts
  (`data_gbif_cleaned.parquet` ~1.2GB, the 22k-species O(n²) medoids
  `climate_centres.parquet`, the AVH↔GBIF crosswalk, the AusTraits flatten) are
  tracked by hash but **not** recomputed by `tar_make()`. Regenerate them
  deliberately with the builders in `R/data/build_sources.R`
  (`build_gbif_climate()`, `build_austraits_climate()`), `load_data()` (height),
  and the retired notebooks (medoids, crosswalk) — see [Retired notebooks](#retired-notebooks).
- *Tier A — computed.* Everything from `unified` (the linchpin, was
  `data_traits_climate_unified.parquet`) forward: per-trait TMB fits (one branch
  per trait via `tar_map`), exports, the combined product, figures, and the two
  rendered reports. This is the part `targets` rebuilds incrementally.

BRMS is opt-in: set `run_brms <- TRUE` in `_targets.R` to add the
`brms_fits`/`brms_diag` targets (height excluded — TMB only).

Supporting files: `demo.qmd` (method illustration — the justification for the within-group-centred approach), `trait_species_means.qmd` (round-2 simple means for low-coverage traits — WIP).

### Retired notebooks

`archive/` was removed from the tree in `eecdbf2`. It held the retired NVIS height report, `notes_why_bayesian.md`, and — still needed to regenerate two Tier B inputs — the **medoid producer** `climate_centres.qmd` and the **AVH↔GBIF crosswalk producer** `Data_height_AVH.qmd`. Recover any of them from history:

```sh
git show daeaa36:archive/climate_centres.qmd > climate_centres.qmd
git show daeaa36:archive/Data_height_AVH.qmd > Data_height_AVH.qmd
```

## R/ library layout

All logic lives in functions under `R/`, loaded via **`source("R/init.R")`**
(sources every function file; deliberately excludes `R/scripts/`).

```
R/
  init.R       # source this once; loads the whole library
  utils/       save_figure.R (png+pdf), parallel_setup.R (multisession)
  data/        build_dataset.R, build_sources.R (regen seeded inputs),
               extract_climate.R, assign_state.R, load_data_height_avh.R
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
- `data-processed/` — seeded pipeline inputs (tracked by `targets` as file
  inputs, not recomputed by `tar_make`): `data_gbif_cleaned.parquet`,
  `climate_centres.parquet` (medoids), `gbif_avh_name_crosswalk.parquet`, and the
  source caches `data_height_climate.parquet` + `data_traits_climate.parquet`;
  plus `brms_fits/` (BRMS fit cache). The unified record-level table is now a
  `targets` object (`unified`) in the `_targets/` store, not a standalone parquet
  (`data_traits_climate_unified.parquet` is the legacy/orphaned cache).
- `outputs/traits/` — per-trait `<trait>_climate_response.csv.gz` (+ `_meta.csv`),
  `species_trait_means_combined.csv.gz`, `trait_model_summary.csv`, figures.

## Gotchas

- **Taxonomy across releases:** AVH, AusTraits and GBIF were aligned to
  different APC releases. Joining medoids/occurrences to trait records needs a
  crosswalk (align all to current APC with `APCalign`); `build_medoid_lookup()`
  bridges AVH/AusTraits names to GBIF medoids via the cached crosswalk. See also the crosswalk build in the retired `Data_height_AVH.qmd` ([Retired notebooks](#retired-notebooks)), and open issue #1 (APC versions out of sync).
- **Medoid memory:** the 2D distance matrix in `locate_medoid_2d()` is O(n²);
  `process_taxon_medoid()` subsamples species with >20k occurrences.
- **`.Rprofile` may set `mc.cores` low** — use `setup_parallel()`.
- Rerunning the full pipeline is expensive but acceptable; outputs are
  regenerated rather than preserved bit-for-bit across refactors.

---

## AusTraits family — cross-package context

`austraits-means` is part of the **AusTraits family** (a subset of the [`traitecoevo`](https://github.com/traitecoevo) org) — here, species-level trait means (observed and climate-centred modelled) built downstream of the released AusTraits database, with taxonomy aligned by `APCalign`. Family-wide concerns are documented centrally in **[austraits-meta](https://github.com/traitecoevo/austraits-meta)** — don't restate them here, read them there:

- **Start with [`AGENTS.md`](https://github.com/traitecoevo/austraits-meta/blob/main/AGENTS.md)** — pipeline order, who owns what, dependency direction, source-of-truth rules, cross-boundary artifacts, gotchas.
- **[`dependencies.yml`](https://github.com/traitecoevo/austraits-meta/blob/main/dependencies.yml)** — machine-readable package graph + cross-boundary artifacts.
- **[`governance/`](https://github.com/traitecoevo/austraits-meta/tree/main/governance)** — label taxonomy, board #9 conventions, release playbooks, triage.

Source-of-truth rules that bite here: trait definitions and units come from **APD**, and taxonomy from **APCalign** — consume them, never redefine them in this repo.

**Filing issues:** the whole family is tracked on one board, [AusTraits #9](https://github.com/orgs/traitecoevo/projects/9) (new issues auto-add to it). Follow the [issue & labelling guide](https://github.com/traitecoevo/austraits-meta/blob/main/governance/issue-guide.md): pick one work-type label (`bug` / `task` / `epic`); Status and Priority are set on the board, not as labels.

**Commit messages:** every family repo squash-merges, so the **PR title and body become the permanent commit message**. Keep the subject ≤50 characters as typed and the body ≤10 lines; put the working detail — what you tried, benchmarks, test counts, rejected alternatives — in the **first PR comment** instead. Full convention: [`commit-messages.md`](https://github.com/traitecoevo/austraits-meta/blob/main/governance/commit-messages.md).

> austraits-meta is hand-maintained prose — a map, not ground truth. Verify specifics against the actual repos.
