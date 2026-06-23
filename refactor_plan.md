# Refactor & extension plan: trait means at species climate centres

Status: **draft for review**. Decisions locked: medoid-centering; per-trait models (joint deferred); vertical-slice-first build order; **drop species with no medoid** (→ round-2 means); **AusTraits field raw/mean/median only**; **BRMS benchmarked on a subset before scaling**; **medoid reasonableness checked** as a data-quality gate.

## 1. Goal

Estimate the mean of each trait at the **centre of each species' climate distribution** (mean annual precipitation + mean annual temperature), from AusTraits + AVH data, using a hierarchical within-group-centred model. Provide **two interchangeable fitting backends**:

- **TMB** (`glmmTMB`) — fast fits for iteration/testing.
- **BRMS** — full posterior for the deliverable.

A later round adds a methodology factor (by `dataset_id`) to partition variance due to differing collection methods, and a joint multivariate model with gap-filling.

## 2. What exists today (assessment)

| Piece | File | State |
|---|---|---|
| Height record table + climate + centring | `R/data/load_data_height_avh.R` | Works, but AVH-only and centres on **trait-sample mean** |
| Fit → tidy predictions | `R/process_fitted_model.R` (`process_fit`) | Works, but **glmmTMB-only** and **precipitation-hardcoded** (`global_c_center`, `annual_precipitation`) |
| 2D climate medoid | `R/locate_medoid_2d.R` | Works; run for ~22k spp → `data-processed/climate_centres.parquet` |
| Species panels | `R/species_panel_functions.R`, `generate_species_panels.R`, `combine_species_panels.R` | Works |
| Round-2 simple means | `R/species_means.R` | Draft |
| Height report | `Data_height_AVH.qmd` | Polished; the template for everything |
| Multi-trait data build | `dataset_construction.qmd` | Partly height-specific; separate from AVH pipeline |
| Old exploratory analysis | `analysis.qmd`, `Data_height_AVH_extra.qmd`, `Sldies.qmd` | Scratch — to archive |
| Bayesian/TMB demos | `demo.qmd` | Reference snippets only (brms/rstan on toy data) |

**Key debt:** (1) no backend abstraction; (2) medoid centres computed but **not wired into fitting**; (3) AVH (height) and AusTraits (other traits) are separate pipelines and must merge — AusTraits carries no height; (4) `process_fit` is single-variable; (5) substantial duplication across `load_data_height_avh.R` / `dataset_construction.qmd` / `analysis.qmd`.

## 2a. Conventions (apply throughout)

> **Update (targets migration):** orchestration now lives in a **`targets`
> pipeline** (`_targets.R`); `tar_make()` rebuilds only what changed. The two
> `.qmd` files are now *reporting layers* rendered by `tar_quarto` targets (they
> `tar_load()` results, they do not compute). The manual `if (update || …)`
> caching below is **retired** — `targets` content-hashing replaces it. Heavy
> upstream artefacts (GBIF clean, medoids, crosswalk, AusTraits flatten) are
> tracked as *file inputs* and regenerated via `R/data/build_sources.R`; see the
> `_targets.R` header and `CLAUDE.md` for the current model. `renv.lock` pins
> dependencies. The notes below are kept for design rationale.

- **Two Quarto files** as reporting layers (rendered by the pipeline):
  - `dataset_construction.qmd` — provenance of the data-building inputs + the medoid reasonableness checks.
  - `analysis.qmd` — modelling write-up, figures, and the round-2 combined product. The current `Data_height_AVH.qmd` / `climate_centres.qmd` content is archived (they remain the producers for the seeded medoid/crosswalk inputs).
- **All logic lives in reusable functions in `R/`.** The `.qmd` files and `_targets.R` only call functions, set parameters, and narrate.
- ~~**Cache expensive steps with `if` guards**: `if (update || !file.exists(out)) { … write } ; read`, one `update <- FALSE` toggle per notebook.~~ Replaced by `targets` (see the update note above).
- **Parallelism via `furrr` + `plan(multisession, workers = n)`** (not `multicore` — forking is unreliable/disabled on macOS/RStudio). Helper to set workers sensibly from `parallelly::availableCores()`.
- **Figures saved as both PNG and PDF** via a `save_figure(plot, name, width, height)` helper (PNG for embedding/preview, PDF as the vector deliverable).

## 3. Target architecture

A layered `R/` library with a single source of truth for each concern, consumed by the **two** Quarto notebooks above.

```
R/
  data/
    align_taxonomy.R      # APC alignment + crosswalk (AVH / AusTraits / GBIF -> current APC)
    extract_climate.R     # shared WorldClim bio_1 + bio_12 raster extraction
    assign_state.R        # assign_state_fast (lifted out of load_data_height_avh.R)
    build_dataset.R       # MERGE AVH height + AusTraits traits -> unified long table
  model/
    model_formula.R       # formula builder: 1-2 climate vars, uncorrelated REs, optional methodology
    fit_trait_climate_model.R  # backend = "tmb" | "brms"; returns standardized fit object
  process/
    process_fit.R         # backend-agnostic tidy predictions (MAT + MAP, any trait)
    species_centres.R     # predicted trait value at each species' medoid + CI
  plot/   plot_fitted.R, plot_fitted_highlight_taxa.R, plot_fitted_progression.R
  panels/ species_panel_functions.R, generate_species_panels.R, combine_species_panels.R
  medoid/ locate_medoid_2d.R, process_taxon_medoid.R, check_medoid.R
  means/  species_means.R          # round-2 fallback for low-coverage traits
  utils/  save_figure.R (png + pdf), parallel_setup.R (multisession workers)
```

The two notebooks call these and nothing else:

```
dataset_construction.qmd : extract_climate -> align_taxonomy -> medoid (+ check) -> build_dataset
                           (each step cached behind if (update || !file.exists()))
analysis.qmd             : load unified data -> fit_trait_climate_model (per trait, furrr)
                           -> process_fit -> species_centres -> exports + figures (png+pdf)
```

### 3.1 Unified record-level dataset (`build_dataset.R`)

One row per trait observation, common schema across sources:

`source` (AVH/AusTraits) · `dataset_id` (study; methodology hook) · `taxon_name` · `original_name` · `family` · `genus` · `trait_name` · `value` · `log_value` · `unit` · `value_type` · `decimal_longitude` · `decimal_latitude` · `state` · `annual_precipitation` · `average_temperature` · `log_annual_precipitation` · then the centre + centred columns (below).

- **Height** from the AVH CSV (`trait_name = "plant_height"`, `source = "AVH"`).
- **Other traits** from AusTraits 7.0.0 flattened, `basis_of_record == "field"`, geo-referenced, `value_type` in {raw, mean, median}. Flora min/max records are excluded here and handled by the round-2 means path.
- Climate via shared `extract_climate.R` (identical layers across sources — already consistent).
- Taxonomy aligned to current APC; crosswalk so the **medoid join works for both sources** (extends the existing AVH↔GBIF crosswalk to AusTraits).
- Cache → `data-processed/data_traits_climate_unified.parquet` with an `update` flag, matching the existing parquet-cache convention.

### 3.2 Centring (medoid, Option B)

Per species, join the medoid from `climate_centres.parquet`:
```
log_annual_precipitation_mn = medoid log precip       # per-species constant
log_annual_precipitation_c  = log_annual_precipitation - medoid
average_temperature_mn      = medoid temp
average_temperature_c       = average_temperature - medoid
```
**No medoid** (no/few GBIF records, or a name mismatch the crosswalk can't resolve) → the species is **dropped from the climate-response models** and routed to the round-2 means path instead. This keeps the centre definition consistent.

### 3.2a Medoid reasonableness check (data-quality gate)

Before a species' medoid is used for centring, verify it actually represents the climate centre and sits sensibly relative to the trait records we're centring on it:

- **Inside its own cloud** — medoid distance to the GBIF occurrence centroid, relative to the occurrence spread (medoid should be near the dense core, not pulled to an edge). The existing `plot_taxon()` already shows medoid vs cloud vs culled outliers for spot checks.
- **Agreement with simple summaries** — medoid MAP/MAT vs the median/mean of occurrences; large divergence flags a skewed/multimodal range or residual outliers.
- **Within the trait-data climate span** — because we set `x_c = obs − medoid`, flag species whose medoid falls far outside the range of their *trait* records' climate (there, all `x_c` share a sign and the centre prediction is an extrapolation).
- Emit a per-species `medoid_ok` flag + a summary plot/table; review before the full run. Species failing the check are reported and excluded (→ round-2 means).

### 3.3 Backend-agnostic fit (`fit_trait_climate_model`)

```r
fit_trait_climate_model(data,
  trait, climate_vars = c("log_annual_precipitation", "average_temperature"),
  backend = c("tmb", "brms"), n_min = 5, methodology = NULL, ...)
```

`model_formula.R` builds, for the chosen climate vars:
```
log_value ~ v1_mn + v2_mn + v1_c + v2_c
            + (1 | taxon_name) + (0 + v1_c | taxon_name) + (0 + v2_c | taxon_name)
```
(uncorrelated REs, no climate interaction — matches the convergence-driven choices in `Data_height_AVH.qmd`). Methodology hook: optional `+ (1 | dataset_id)` or fixed methodology covariates, off by default.

Returns a **standardized fit object** (S3) so downstream code is backend-blind:
`backend` · `formula` · `data` · `fixed` (term/estimate/conf.low/conf.high) · `ranef` (taxon_name, g_0, g_c per var) · `family` · `predict_at(newdata)` closure · raw `fit`. TMB path uses `glmmTMB` + `marginaleffects`; BRMS path uses posterior draws / `fitted`.

**BRMS specifics:** consult the `bayesian-workflow` skill before writing — priors, prior/posterior predictive checks, divergences, R-hat/ESS, **non-centered parameterization** (essential for thousands of species), LOO. Compute is non-trivial: a hierarchical model over ~10k species is a serious Stan job. Approach: **get BRMS working and validated against TMB on a species subset during the vertical slice, benchmark runtime, then scale** to full data with `cmdstanr` + non-centered parameterization + threading; full multi-trait BRMS runs may be batched.

### 3.4 Processing & outputs

- `process_fit.R` generalized to MAT+MAP and arbitrary trait names (remove hardcoded `global_c_center` / `annual_precipitation`).
- `species_centres.R` → predicted trait value at each species' medoid centre + CI → export `…/<trait>_climate_response.csv.gz` + a `_meta.csv` data dictionary, mirroring the height deliverable.
- All produced inside `analysis.qmd`, looping over traits; figures written via `save_figure()` (png+pdf). The NVIS write-up from `Data_height_AVH.qmd` becomes the height section of `analysis.qmd` (its prose preserved), generalized across traits.

## 4. Phased delivery (for review)

**Phase 0 — Scaffolding (no behaviour change)**
- Create `CLAUDE.md` (pipeline order, data locations, centring convention, backend abstraction, two-file rule, caching/furrr/figure conventions, gotchas).
- Add `R/utils/` helpers: `save_figure()` (png+pdf) and `parallel_setup()` (multisession workers).
- Reorganize `R/` into the layered structure; extract shared `extract_climate.R` / `assign_state.R`; de-duplicate. Pure move/refactor — existing height results must reproduce bit-for-bit.

**Phase 1 — Vertical slice: height on the new architecture**
- `build_dataset.R` for height (AVH) with medoid centring + taxonomy crosswalk; drop species with no medoid.
- **Medoid reasonableness check** (§3.2a): per-species `medoid_ok` flag + summary plot/table; exclude failures.
- `model_formula.R` + `fit_trait_climate_model(backend = "tmb")` reproducing current height fixed effects.
- Generalize `process_fit.R`; add `species_centres.R` at the medoid.
- Add `backend = "brms"` on a **species subset**; **confirm TMB ≈ BRMS**; benchmark runtime.
- **Validate medoid- vs trait-mean-centring** on height (params + species centres).
- Regenerate height outputs/report through the new path.

**Phase 2 — Generalize to all traits**
- Extend `build_dataset.R` to AusTraits multi-trait and merge height → full unified table (in `dataset_construction.qmd`).
- Per-trait loop in `analysis.qmd` (TMB default; BRMS opt-in) with `n_min` filter, parallelized with furrr.
- Per-trait exports + figures (png+pdf).

**Phase 3 — Round 2 + methodology hook + cleanup**
- `species_means.R` fallback for low-coverage traits (and species dropped for lacking a medoid).
- Plumb `dataset_id` (one per AVH collector; AusTraits study id) so the methodology term can be switched on in the model.
- Collapse everything into the two notebooks; archive `Data_height_AVH.qmd`, `Data_height_AVH_extra.qmd`, `climate_centres.qmd`, `Sldies.qmd`, and the old `analysis.qmd` content; finalize `CLAUDE.md`.

## 5. Resolved decisions

- **Centring:** medoid (Option B) — `x_c = obs − medoid`, `x_mn = medoid`.
- **Traits:** per-trait models now; joint multivariate/gap-filling in round 2.
- **Build order:** vertical slice (height) first, then generalize.
- **No medoid:** drop from models → round-2 means.
- **AusTraits records:** field, `value_type` ∈ {raw, mean, median}.
- **BRMS:** validated vs TMB (centres r=1.000); production path set up — cmdstanr + within-chain threading + on-disk fit cache (`fits_dir`) + per-fit diagnostics. Benchmarked; larger traits need `iter` ≥ 4000 to converge. Opt-in in `analysis.qmd`.
- **Medoid quality:** reasonableness check + `medoid_ok` flag (§3.2a).
- **AVH methodology:** one `dataset_id` per collector.
- **Structure:** `targets` pipeline (`_targets.R`) orchestrates; the two notebooks (`dataset_construction.qmd`, `analysis.qmd`) are `tar_quarto` reporting layers; all logic in `R/` functions; caching via `targets` content-hashing (the old `if`-guarded `update` toggles are retired); `renv.lock` pins deps; figures saved png + pdf.
- **Trait list:** all traits — plant height (AVH) + the AusTraits set (LMA, leaf area, leaf N/P per mass & area, seed mass, wood density). `leaf_delta13C` deferred (needs a discrimination transform before inclusion).
- **Framing:** NVIS dropped (it was a reportable subproject; the height report is archived). Two framings now: (1) national species-trait estimates under the ARDC AusTraits programme; (2) a paper reporting the species means and giving credence to the estimates.

## 6. Open questions

_None outstanding for rounds 1–2. Round-3 methodology-factor specifics (how to
group `dataset_id` collectors, fixed vs random) to be decided when that work starts._
```
