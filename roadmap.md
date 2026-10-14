# Roadmap: species means and community layers

This is the planning document for `austraits-means`. It records what the repository delivers, the decisions taken and why, the milestones and the issues that reach them, and the invariants a future edit could break. Keep it current or retire it; do not start a second plan beside it. Issue numbers are on the [AusTraits board](https://github.com/orgs/traitecoevo/projects/9) (view "means").

Last revised: 14 October 2026.

## Products

Three data products, defined in `terminology.md`:

1. **Observed means**: per species × trait, computed directly from the records. Two methods exist. The *portal method* (`method = portal_arithmetic`) reproduces exactly the species averages the AusTraits data portal computes for AusTraits 7.0.0, so that product can be archived and cited. The *geometric method* (`method = geometric`) is its successor from the next AusTraits build onward.
2. **Modelled means**: per species × trait, the model-estimated geometric mean at the centre of the species' climate distribution (the within-species-centred model in `R/model/`), built here.
3. **Community layers**: the community mean and community variance of the species means, for woody and non-woody plants separately, per spatial unit across Australia. First as observed statistics; later modelled across space.

Each product carries the **same version tag as the AusTraits release it derives from** (for example "species means 7.0.0"), so the two are easily related. A re-issue for the same release appends a fourth component (7.0.0.1). Products derived from a development build carry that build's version string until the corresponding release exists.

## Milestones

| Milestone | Due | What it delivers |
|---|---|---|
| first draft | 3 Nov 2026 | A first version of all three products, for review: portal-method observed means reproduced and documented; geometric observed means and modelled means on the latest austraits.build development build; community layers on provisional choices; an observed-vs-modelled comparison; the decision issues opened. |
| release candidate | 23 Nov 2026 | The revised products after review: decisions closed and applied; refreshed occurrence data and climate medoids; products packaged and published; the portal reading the published product; hooks filed in the sibling repositories. |
| paper submission | undated | Model extensions (unconverged BRMS refits, collection method as a covariate, a joint multivariate model), a product-description methods and figure set, modelled community layers. |

The milestones are review points, not product versions. Both deliver everything; the second delivers it revised.

### first draft (3 Nov 2026)

Epics and issues, with dependencies:

- Epic [#6](https://github.com/traitecoevo/austraits-means/issues/6) Observed means, portal method: [#7](https://github.com/traitecoevo/austraits-means/issues/7) port the portal method (pinned to the commit that built the deployed data) → [#8](https://github.com/traitecoevo/austraits-means/issues/8) reproduction test → [#9](https://github.com/traitecoevo/austraits-means/issues/9) workflow document → [#10](https://github.com/traitecoevo/austraits-means/issues/10) draft Zenodo record. [austraits.portal#139](https://github.com/traitecoevo/austraits.portal/issues/139) (portal repo) confirms provenance and restores the core-trait list; it does not block the test. [#7](https://github.com/traitecoevo/austraits-means/issues/7) also feeds [#15](https://github.com/traitecoevo/austraits-means/issues/15) and [#24](https://github.com/traitecoevo/austraits-means/issues/24).
- Epic [#11](https://github.com/traitecoevo/austraits-means/issues/11) Data basis: [#12](https://github.com/traitecoevo/austraits-means/issues/12) read trait and height records from an austraits.build development build and retire the AVH CSV loader → [#13](https://github.com/traitecoevo/austraits-means/issues/13) rebuild the name crosswalk with the APC release pinned → [#14](https://github.com/traitecoevo/austraits-means/issues/14) refit modelled means and compare with the June 2026 outputs.
- [#15](https://github.com/traitecoevo/austraits-means/issues/15) geometric observed means on the development build → [#16](https://github.com/traitecoevo/austraits-means/issues/16) reliability indicators; [#17](https://github.com/traitecoevo/austraits-means/issues/17) rename exports to the agreed terminology.
- Decisions for the community layers, taken by the project lead: [#18](https://github.com/traitecoevo/austraits-means/issues/18) weighting (decides the names), [#19](https://github.com/traitecoevo/austraits-means/issues/19) woody / non-woody assignment, [#20](https://github.com/traitecoevo/austraits-means/issues/20) spatial unit, minimum species, occurrence source and species-means input. Wanted before the first-draft milestone. Precedent to weigh in all three: Andrew et al. 2025, *Global Ecology and Biogeography* 34: e70136 ([doi:10.1111/geb.70136](https://doi.org/10.1111/geb.70136)), which mapped community trait means and variability for Australia from gap-filled AusTraits species means and species distribution models, on a 10 km equal-area grid, and validated them against AusPlots; compared choice by choice in `docs/community-layers-design.md`.
- Epic [#21](https://github.com/traitecoevo/austraits-means/issues/21) Community layers: [#22](https://github.com/traitecoevo/austraits-means/issues/22) provisional woody lookup → [#23](https://github.com/traitecoevo/austraits-means/issues/23) first draft of the layers on a provisional 0.25° grid, with a design document.
- [#24](https://github.com/traitecoevo/austraits-means/issues/24) observed-vs-modelled comparison → [#25](https://github.com/traitecoevo/austraits-means/issues/25) figure set and coverage summaries → [#26](https://github.com/traitecoevo/austraits-means/issues/26) refresh reports, README and AGENTS.md.

Dependencies: the longest chain is [#12](https://github.com/traitecoevo/austraits-means/issues/12) → [#13](https://github.com/traitecoevo/austraits-means/issues/13) → [#14](https://github.com/traitecoevo/austraits-means/issues/14) → [#23](https://github.com/traitecoevo/austraits-means/issues/23) → [#25](https://github.com/traitecoevo/austraits-means/issues/25); [#23](https://github.com/traitecoevo/austraits-means/issues/23) also needs [#7](https://github.com/traitecoevo/austraits-means/issues/7) → [#15](https://github.com/traitecoevo/austraits-means/issues/15) → [#22](https://github.com/traitecoevo/austraits-means/issues/22). Alongside: [#7](https://github.com/traitecoevo/austraits-means/issues/7) → [#8](https://github.com/traitecoevo/austraits-means/issues/8) → [#9](https://github.com/traitecoevo/austraits-means/issues/9) → [#10](https://github.com/traitecoevo/austraits-means/issues/10); [#17](https://github.com/traitecoevo/austraits-means/issues/17) and [#24](https://github.com/traitecoevo/austraits-means/issues/24) → [#26](https://github.com/traitecoevo/austraits-means/issues/26). With one person doing most of the work, total effort rather than chain length is the constraint, and there is no slack. If it does not fit, slip in this order: [#16](https://github.com/traitecoevo/austraits-means/issues/16), [#10](https://github.com/traitecoevo/austraits-means/issues/10), [#25](https://github.com/traitecoevo/austraits-means/issues/25) beyond the comparison.

### release candidate (23 Nov 2026)

- [#27](https://github.com/traitecoevo/austraits-means/issues/27) woody lookup per the agreed rule → [#28](https://github.com/traitecoevo/austraits-means/issues/28) community layers per the agreed weighting, grid and inputs.
- [#32](https://github.com/traitecoevo/austraits-means/issues/32) fresh GBIF occurrence download, cleaning and climate extraction → [#33](https://github.com/traitecoevo/austraits-means/issues/33) align names to the pinned APC release, rebuild the crosswalk and the climate medoids (closes [#1](https://github.com/traitecoevo/austraits-means/issues/1)) → [#34](https://github.com/traitecoevo/austraits-means/issues/34) refit modelled means and community layers.
- [#29](https://github.com/traitecoevo/austraits-means/issues/29) address review comments.
- [#30](https://github.com/traitecoevo/austraits-means/issues/30) package and publish the products (Zenodo deposit "AusTraits species means").
- [#31](https://github.com/traitecoevo/austraits-means/issues/31) community consultation on the spatial layers, scheduled.
- Cross-repository work, tracked in those repositories without a milestone here: [austraits.portal#139](https://github.com/traitecoevo/austraits.portal/issues/139) and [austraits.portal#140](https://github.com/traitecoevo/austraits.portal/issues/140) (portal: provenance, then reading the published product), [austraits.build#867](https://github.com/traitecoevo/austraits.build/issues/867) (austraits.build: products in the release workflow), [austraits-api#16](https://github.com/traitecoevo/austraits-api/issues/16) (austraits-api: serve stored means), [austraits#184](https://github.com/traitecoevo/austraits/issues/184) (austraits R package: loader).

### paper submission (undated)

[#35](https://github.com/traitecoevo/austraits-means/issues/35) BRMS refits for unconverged traits; [#3](https://github.com/traitecoevo/austraits-means/issues/3) collection method as a covariate, with [#36](https://github.com/traitecoevo/austraits-means/issues/36) importing the factor tables and testing on taxa measured under more than one method; [#37](https://github.com/traitecoevo/austraits-means/issues/37) methods description and figure set for the products; [#38](https://github.com/traitecoevo/austraits-means/issues/38) modelled community means across space. [#40](https://github.com/traitecoevo/austraits-means/issues/40) joint multivariate model with gap-filling. Deferred without a milestone: [#39](https://github.com/traitecoevo/austraits-means/issues/39) leaf_delta13C transform.

## Decisions and their reasons

Dated decisions, newest first. Earlier decisions from the refactor plan are kept below them.

### 14 October 2026

- **The term is "mean", not "average".** Both species estimates are geometric means, so "mean" is accurate; "average" is ambiguous and clashes with `average_temperature`. Flexibility of method is carried by the `method` column, not by the term. An earlier preference for "average" (to leave room for medians and trimmed means) was considered and set aside.
- **The portal-method observed mean reproduces the deployed portal output for AusTraits 7.0.0.** The portal computes its species averages offline (`estimate_species_trait_means()` in `austraits.portal`, at the commit that built the deployed data); that is what users have been seeing, so it is the product to archive. Exactness is at key level (taxon × trait × unit) with numeric tolerance; row order is excluded because the portal has no lockfile.
- **That product is an arithmetic mean and is labelled as such.** The geometric observed mean is the successor method from the next AusTraits build onward, with the arithmetic mean kept in its own column. The portal method's quirks (a geometric-mean column that equals the mean; no life-stage or basis-of-record filter on raw/mean/median records; lone flora bounds; no taxonomic roll-up) are documented, not fixed, in that product.
- **Derived products get their own Zenodo deposit**, "AusTraits species means", with its own concept DOI and one version per AusTraits release (or per development build, tagged with its version string), each linked to the release or build it derives from. Alternative considered: re-releasing AusTraits 7.0.0 with the product added. Set aside because a derived product is not a database release and the two would drift.
- **Versioning follows the AusTraits release tag** (see Products).
- **Trait and height records come from an austraits.build development build** (which now includes the AVH height dataset), replacing the separate AVH CSV. The GBIF refresh and medoid rebuild follow in the release candidate; the first draft uses the existing occurrence extract and medoids, with the crosswalk rebuilt from the new names.
- **Unconverged BRMS fits stay postponed** to the paper milestone; TMB estimates serve the first two milestones.
- **The name crosswalk is built inside the pipeline, against APC release 2026-07-30** (`apc_version` in `_targets.R`), from every taxon name in the development build and every GBIF name ([#13](https://github.com/traitecoevo/austraits-means/issues/13)). Changing the release rebuilds everything downstream. The release candidate re-pins it when the GBIF extract is refreshed ([#33](https://github.com/traitecoevo/austraits-means/issues/33)).
- **The field-only per-record geometric mean (`estimate_simple_means()`) is retired** ([#17](https://github.com/traitecoevo/austraits-means/issues/17)). The geometric observed mean ([#15](https://github.com/traitecoevo/austraits-means/issues/15)) supersedes it as the observed mean, and the combined product (modelled where available, otherwise observed) is built from that. Exports carry `estimate_type` and `method` and name their values `value`, `log_value`, `value_lower`, `value_upper` (`terminology.md`).
- **Infraspecific taxa are rolled up to species** for the modelled means and the geometric observed means (`roll_up_to_species()`; [#13](https://github.com/traitecoevo/austraits-means/issues/13), [#15](https://github.com/traitecoevo/austraits-means/issues/15)). The AVH height dataset in austraits.build keeps subspecies, varieties and forms, which the earlier AVH CSV had rolled up to species, and GBIF names (and so medoids) are almost all at species rank. Rolling up gives every record a species-level climate centre and makes both species estimates join on species. The build's own name is kept in `taxon_name_aligned`. Alternative considered and set aside: keeping infraspecific taxa as their own groups, centred on their parent species' medoid. The portal-method product is not rolled up; that is one of its documented quirks.
- **Weighting and woody / non-woody assignment for the community layers are decisions, not defaults.** The first draft uses provisional choices, flagged as such; the release candidate applies the decisions.
- **The spatial layers get a community consultation**, scheduled after the first-draft review.

### Earlier (refactor plan, 2026)

- Centring on the medoid: `x_c = obs − medoid`, `x_mn = medoid` per climate variable; consistent across traits.
- Model structure: `log_value ~ v1_mn + v2_mn + v1_c + v2_c + (1 | taxon) + (0 + v1_c | taxon) + (0 + v2_c | taxon)`, with v1 = log annual precipitation (WorldClim BIO12) and v2 = mean annual temperature (BIO1). Random slopes are uncorrelated and there is no climate interaction; both choices were driven by convergence. Species need at least 5 records (`n_min = 5`).
- One fitting interface for both backends: `fit_trait_climate_model(backend = "tmb" | "brms")` returns a standard fit object, so downstream code is backend-blind. BRMS uses a non-centred parameterisation. Parallelism is `furrr` multisession, never multicore.
- Per-trait models now (TMB default, BRMS for the full posterior); joint multivariate and gap-filling deferred.
- Species with no medoid are dropped from the climate-response models and routed to the observed means.
- AusTraits records for the modelled means: `basis_of_record == "field"`, `value_type` in raw, mean, median; per-dataset exceptions are allowed where a dataset carries only bounds (the AVH height dataset).
- BRMS validated against TMB (centres r = 1.000); production path set up (cmdstanr, threading, on-disk cache, diagnostics); larger traits need `iter ≥ 4000`.
- Medoid quality gated by a reasonableness check (`medoid_ok`): the medoid sits inside its own occurrence cloud, agrees with simple summaries (median and mean of occurrences), and falls within the climate span of the species' trait records. Species that fail are excluded from the models and get only an observed mean.
- Height methodology: one `dataset_id` per AVH collector was carried for a later methodology factor; with the move to the austraits.build dataset this becomes a single `dataset_id`, and the collector is recoverable from the record identifier if needed.
- Orchestration by `targets`, dependencies pinned by `renv`, all logic in `R/` functions, notebooks as reporting layers, figures saved as PNG and PDF.
- Trait list: plant height plus LMA, leaf area, leaf N and P per mass and per area, seed mass, wood density. `leaf_delta13C` deferred.

## Invariants

Things a future edit could break without noticing:

- The portal-method observed mean must stay reproducible from the AusTraits 7.0.0 release and the pinned portal commit. Do not "fix" its quirks; they are the product. Fixes go in the geometric method.
- Value-type rules are per dataset: datasets that carry only minimum/maximum records (the AVH height dataset) enter the modelled means by an explicit rule, not by the default raw/mean/median filter.
- Centring is on the GBIF medoid, never on the trait-sample mean.
- Both species estimates are geometric means from the geometric method onward; an arithmetic mean is always a separately labelled column.
- Column names follow `terminology.md`: `estimate_type` and `method`, never a bare `estimate` or `geom_mean`.
- Community statistics are unweighted unless the weighting decision says otherwise; unweighted values are never called CWM or CWV.
- Use "mean", never "average", in code, exports, docs and metadata.

## Later work

On the paper milestone:

- Modelled community layers across space ([#38](https://github.com/traitecoevo/austraits-means/issues/38)), after the observed layers are settled.
- Collection-method covariate ([#3](https://github.com/traitecoevo/austraits-means/issues/3), [#36](https://github.com/traitecoevo/austraits-means/issues/36)): waits for the factor tables and a test on taxa measured under more than one method.
- Joint multivariate model with gap-filling ([#40](https://github.com/traitecoevo/austraits-means/issues/40)).

Deferred, no milestone:

- `leaf_delta13C` ([#39](https://github.com/traitecoevo/austraits-means/issues/39)): needs the discrimination transform.
