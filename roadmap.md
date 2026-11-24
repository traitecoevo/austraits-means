# Roadmap: modelled species means and community means from climate

This is the planning document for `austraits-means`. It records what the repository delivers, the decisions taken and why, the milestones and the issues that reach them, and the invariants a future edit could break. Keep it current or retire it; do not start a second plan beside it. Issue numbers are on the [AusTraits board](https://github.com/orgs/traitecoevo/projects/9) (view "means").

Last revised: 24 November 2026.

## Products

Re-scoped 24 November 2026 (decisions below). This repository publishes one product, in a Zenodo record of its own, for a short list of traits; terms are defined in `terminology.md`:

1. **Modelled means**: per species × trait, the model-estimated geometric mean at the centre of the species' climate distribution (the within-species-centred model in `R/model/`).
2. **Modelled community statistics**, in the same record: for a location, predicted from its climate alone, the woody and non-woody **community mean** (mean trait across species) and **community range** (spread of the species' trait means there). They are a function of climate, which can be evaluated at any location; they are not maps of the records. The model that gives an honest community range is not built yet (Later work).

Not published from here:

- **Observed means** (the simple species means) are published by austraits.build with each AusTraits release and its Zenodo record ([austraits.build#899](https://github.com/traitecoevo/austraits.build/pull/899)). This repository reads them.
- The **observed community layers** (`observed_local`, `observed_presence`) and the current GBIF-presence `modelled` layer stay in the pipeline as unpublished diagnostics: the within-cell spread of `observed_local` is the benchmark for the modelled community range.
- The **portal-method observed means for AusTraits 7.0.0** are retired, unpublished. The code and its reproduction check are kept in the pipeline but not maintained.

The trait list is open until the location-model prototype has been checked trait by trait (Later work). How the new record is named and versioned is open.

## Milestones

| Milestone | Due | What it delivers |
|---|---|---|
| first draft | 3 Nov 2026 | A first version of all three products, for review: portal-method observed means reproduced and documented; geometric observed means and modelled means on the latest austraits.build development build; community layers on provisional choices; an observed-vs-modelled comparison; the decision issues opened. |
| release candidate | 23 Nov 2026 | The revised products after review: decisions closed and applied; refreshed occurrence data and climate medoids; products packaged and published; the portal reading the published product; hooks filed in the sibling repositories. |
| extensions | undated | Model extensions (unconverged BRMS refits, collection method as a covariate, a joint multivariate model), a product-description methods and figure set, modelled community layers. |

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

Status, 14 October 2026: every first-draft issue is done except [#10](https://github.com/traitecoevo/austraits-means/issues/10), whose draft deposit is created but still uploading, and the decisions [#18](https://github.com/traitecoevo/austraits-means/issues/18)–[#20](https://github.com/traitecoevo/austraits-means/issues/20), which stay with the project lead. The community layers were built on the provisional choices, in three variants for review (`combined`, `modelled_centre`, `model_surface`), since replaced by the three products of [#28](https://github.com/traitecoevo/austraits-means/issues/28) (`docs/community-layers-design.md`). Nothing slipped.

### release candidate (23 Nov 2026)

- [#27](https://github.com/traitecoevo/austraits-means/issues/27) woody lookup per the agreed rule → [#28](https://github.com/traitecoevo/austraits-means/issues/28) community layers per the agreed weighting, grid and inputs.
- [#32](https://github.com/traitecoevo/austraits-means/issues/32) fresh GBIF occurrence download, cleaning and climate extraction → [#33](https://github.com/traitecoevo/austraits-means/issues/33) align names to the pinned APC release, rebuild the crosswalk and the climate medoids (closes [#1](https://github.com/traitecoevo/austraits-means/issues/1)) → [#34](https://github.com/traitecoevo/austraits-means/issues/34) refit modelled means and community layers.
- [#29](https://github.com/traitecoevo/austraits-means/issues/29) address review comments.
- [#30](https://github.com/traitecoevo/austraits-means/issues/30) package and publish the products (Zenodo deposit "AusTraits species means").
- [#31](https://github.com/traitecoevo/austraits-means/issues/31) community consultation on the spatial layers, scheduled.
- Cross-repository work, tracked in those repositories without a milestone here: [austraits.portal#139](https://github.com/traitecoevo/austraits.portal/issues/139) and [austraits.portal#140](https://github.com/traitecoevo/austraits.portal/issues/140) (portal: provenance, then reading the published product), [austraits.build#867](https://github.com/traitecoevo/austraits.build/issues/867) (austraits.build: products in the release workflow), [austraits-api#16](https://github.com/traitecoevo/austraits-api/issues/16) (austraits-api: serve stored means), [austraits#184](https://github.com/traitecoevo/austraits/issues/184) (austraits R package: loader).

Status, 20 October 2026:

- Done: the decisions [#18](https://github.com/traitecoevo/austraits-means/issues/18)–[#20](https://github.com/traitecoevo/austraits-means/issues/20) (recorded below) and [#27](https://github.com/traitecoevo/austraits-means/issues/27); [#32](https://github.com/traitecoevo/austraits-means/issues/32) (GBIF download `0027765-260928105237408`, CC0 and CC BY records only) and [#33](https://github.com/traitecoevo/austraits-means/issues/33) (names aligned to APC 2026-07-30 from the providers' names, crosswalk and medoids rebuilt; closed [#1](https://github.com/traitecoevo/austraits-means/issues/1)). The data basis is the austraits.build development build at `d535de9` with APCalign 2.0.1.
- [#28](https://github.com/traitecoevo/austraits-means/issues/28): the three community products are built, each for woody, non-woody and all species, with cell-by-cell comparisons and a native-only test (`docs/community-layers-design.md`). Open for review, chiefly of the modelled layers' woodiness term, which was chosen over separate models per group.
- Done: [#34](https://github.com/traitecoevo/austraits-means/issues/34). The modelled means were refitted on the new medoids in #33 (fixed effects moved by at most 0.46 standard errors). Against layers rebuilt from the 2025 GBIF extract and its medoids on the same 10 km grid, everything else held fixed, the presence and modelled community means correlate cell by cell at 0.94–1.00 (log scale), with median shifts of at most 0.04 log units: the refresh's main effect is the licence, not the layers.
- [#10](https://github.com/traitecoevo/austraits-means/issues/10) (first draft) is still open: the draft deposit exists on Zenodo but two files are still to upload (`observed_means_portal_arithmetic.csv.gz` and `.parquet`). It is never published from here.
- [#29](https://github.com/traitecoevo/austraits-means/issues/29)–[#31](https://github.com/traitecoevo/austraits-means/issues/31) wait on the review, the go-ahead to publish (licensing is settled: CC BY 4.0, as the NC records are dropped) and the consultation.
- Per-trait figures of within- vs between-species responses are now drawn for every trait (`outputs/traits/figures/<trait>_within_between`). Collection method as a covariate ([#3](https://github.com/traitecoevo/austraits-means/issues/3), [#36](https://github.com/traitecoevo/austraits-means/issues/36)) is untested: the model's `methodology` hook (a `(1 | dataset_id)` intercept) has never been run, height has a single `dataset_id` (`AVH_2026`), and the method-factor tables are not yet imported.

Status, 24 November 2026: the due date passed with the products re-scoped (decisions below), so this milestone does not close as written.

- Closed as not planned: [#10](https://github.com/traitecoevo/austraits-means/issues/10) and its epic [#6](https://github.com/traitecoevo/austraits-means/issues/6): the portal-method product is retired. The draft deposit on Zenodo is left for the maintainer to discard.
- Re-scoped:
  - [#28](https://github.com/traitecoevo/austraits-means/issues/28) and its epic [#21](https://github.com/traitecoevo/austraits-means/issues/21): the observed layers are built and become unpublished diagnostics; no further work on them for publication.
  - [#30](https://github.com/traitecoevo/austraits-means/issues/30): publishes the modelled means and modelled community statistics in a new Zenodo record.
  - [#38](https://github.com/traitecoevo/austraits-means/issues/38): becomes the community statistics predicted from climate, by the model with location effects (new issue [#41](https://github.com/traitecoevo/austraits-means/issues/41)).
  - [#31](https://github.com/traitecoevo/austraits-means/issues/31): whether a consultation is still wanted for climate-derived statistics is open.
- Next: the location-model prototype ([#41](https://github.com/traitecoevo/austraits-means/issues/41)), then the trait list, then [#29](https://github.com/traitecoevo/austraits-means/issues/29) and [#30](https://github.com/traitecoevo/austraits-means/issues/30). No new due date until the prototype reports.

### extensions (undated)

[#35](https://github.com/traitecoevo/austraits-means/issues/35) BRMS refits for unconverged traits; [#3](https://github.com/traitecoevo/austraits-means/issues/3) collection method as a covariate, with [#36](https://github.com/traitecoevo/austraits-means/issues/36) importing the factor tables and testing on taxa measured under more than one method; [#37](https://github.com/traitecoevo/austraits-means/issues/37) methods description and figure set for the products; [#38](https://github.com/traitecoevo/austraits-means/issues/38) modelled community means across space. [#40](https://github.com/traitecoevo/austraits-means/issues/40) joint multivariate model with gap-filling. Deferred without a milestone: [#39](https://github.com/traitecoevo/austraits-means/issues/39) leaf_delta13C transform.

## Decisions and their reasons

Dated decisions, newest first. Earlier decisions from the refactor plan are kept below them.

### 24 November 2026

Product re-scope, decided by Daniel Falster:

- **Simple (observed) means are published by austraits.build with each AusTraits release and its Zenodo record, not here.** They are a direct summary of the release, and the API and portal can then serve one copy. Supersedes the 14 October decision that derived products share one "AusTraits species means" deposit.
- **The portal-method observed means for AusTraits 7.0.0 are retired, unpublished** ([#6](https://github.com/traitecoevo/austraits-means/issues/6), [#10](https://github.com/traitecoevo/austraits-means/issues/10) closed). The successor method ships with the release, so a separate archive of the 7.0.0 portal values adds little. Their code and reproduction check stay in the pipeline, unmaintained. Alternatives considered: finishing the #10 deposit as a one-off archive; handing the table to austraits.build to sit with 7.0.0.
- **Modelled means are a separate product, for far fewer traits, in a new Zenodo record.** Which traits is decided after the location-model prototype, by how well each trait supports a community range. Evidence so far is in Later work.
- **No observed community layers are published: there are too few records.** A trait record collected in a cell is rare: 72–459 sites have 5 or more species measured together for the eight leaf, seed and wood traits (Later work). Spatial output is predicted from climate by the trait model instead. The observed layers already built stay as unpublished diagnostics, because their within-cell spread is the benchmark for the modelled community range. Supersedes the 19 October decision on three community products. Alternative considered: removing them.
- **The community statistics are, for a location and from its climate, the community mean and the community range, for woody and non-woody plants separately** (`terminology.md`). There is no "all species" statistic: mixing the two groups needs each location's woody share, which comes from occurrence data.
- **Community variance is between species only, within one location and within one growth form** (woody or non-woody), so it excludes the difference between the groups as well as variation within species.
- **The community range comes from a model with location effects, never from the pooled species variance.** The current model pools all records from the same climate, so the spread among its species includes how locations with the same climate differ; the spread at one location is smaller (measured: 0.59–0.92 of it, Later work). The next modelling step is a model that separates the two ([#41](https://github.com/traitecoevo/austraits-means/issues/41)).

Observed means:

- **Signed traits get an arithmetic observed mean** (traits.use#51). A trait whose definition allows values below zero (water potentials, temperatures, isotope ratios; 43 in the development build) has no geometric mean, so its `value` is the arithmetic mean of the same groups, labelled `method = arithmetic`. Alternative considered: leaving such traits without an observed mean, as the first version did (23 traits had none).
- **The calculation moves to traits.use** as `summarise_trait_means()` (traits.use#51), after Lizzy Wenk's review: usable on any traits.build database or subset, every trait by default, with `by` for further grouping (AusInverTraits splits by sex). austraits.build#900 calls it.

- **The geometric observed means ship with each AusTraits release, built in austraits.build, not here** ([austraits.build#899](https://github.com/traitecoevo/austraits.build/pull/899)). `post_export()` in austraits.build's `R/post_export.R`, which `traits.build::build_post_export()` runs from `create_release()`, writes `austraits-<version>-observed-means.parquet`, `.csv.gz` and `-dictionary.csv` into the release. The method is this repository's geometric method, ported to austraits.build's `R/observed_means.R`; on the development build `d535de9` the port reproduces `observed_means_geometric` exactly (602,644 rows), less the `source_primary_key` column. This repository keeps the portal-method product, the modelled means and the community layers. Alternative considered: publishing the observed means only in the separate "AusTraits species means" Zenodo deposit; set aside because the simple means are a direct summary of the release and belong with it.
- **The method has one home: austraits.build.** This repository's copy (`R/means/observed_means_geometric.R`) is deleted; the `observed_means_geometric` target reads the file austraits.build's `post_export()` writes for the same build (`austraits_dev_observed_means` in `_targets.R`; `data-raw/austraits/README.md`), and no longer writes `outputs/species_means/<build>/observed_means_geometric.*`. A change to the method is made in austraits.build and reaches here through a new file. `roll_up_to_species()` here, used by the modelled means and the community layers, must keep rolling up the same way as austraits.build's `observed_means_records()`, so the two species estimates join.

### 19 October 2026

- **Community statistics are presence-only: every species present in a unit counts equally**, so they are the *community mean* and *community variance*, never CWM or CWV ([#18](https://github.com/traitecoevo/austraits-means/issues/18)). Abundance-weighted statistics from survey plots were considered; they cover far fewer sites and species.
- **Community layers use a 10 × 10 km equal-area grid, as Andrew et al. 2025** (Lambert Azimuthal Equal Area, `+proj=laea +lat_0=-25.2744 +lon_0=133.7751 +x_0=0 +y_0=0 +datum=WGS84 +units=m`), **at least 10 species with a mean per cell and group, and GBIF occurrence points for presence** ([#20](https://github.com/traitecoevo/austraits-means/issues/20)). Range maps from species distribution models were considered as the occurrence source.
- **A species' woodiness and growth form are their most likely value as the AusTraits data portal computes it** (master's categorical summary: each space-separated value counted once per record, the most frequent kept, ties joined); the geometric observed means' `value_mode` is that same value ([#19](https://github.com/traitecoevo/austraits-means/issues/19)).
- **Woody: `woodiness_detailed` woody or woody_like_stem (palms, cycads, grass trees); non-woody: herbaceous and the other intermediate categories (semi_woody, woody_base, woody_like_inflorescence, woody_root); then growth-form keywords; species still unassigned are left out of the woody and non-woody layers, with no fallback to the coarse `woodiness` trait** ([#19](https://github.com/traitecoevo/austraits-means/issues/19)). Unassigned species are few (23 in the development build) and unclear cases.
- *(Superseded 24 November 2026: none of these is published; they stay as diagnostics.)* **Three community products, replacing the first draft's variants** ([#20](https://github.com/traitecoevo/austraits-means/issues/20), [#28](https://github.com/traitecoevo/austraits-means/issues/28)): *observed, local records*, where each species' records in a cell are averaged and then the species, with no model; *observed, GBIF presence*, the species present in a cell each given its geometric observed mean; and *modelled*, the per-trait model's net average for the cell, b0 + b_within·x + (b_between − b_within)·m̄, from the cell's climate x and the mean climate centre m̄ of the species present, using no species means. The first keeps values tied to where traits were measured; the second gives coverage; the third averages over between- and within-species effects. All three are built for woody and non-woody species separately, and for all species. The modelled product is split through a woodiness term in each trait's model (woody and non-woody intercepts and climate slopes, species random effects shared), so the layer for all species mixes the two by their species in the cell; without it, the model, fitted mostly on woody species, could not follow the woody share of a cell. Separate models per group were the alternative; the shared random effects help the groups with few species. A group with fewer than 20 fitted species gets no modelled layer (non-woody leaf N and P per area, wood density).

- **Only openly licensed GBIF records are used (CC0 and CC BY), so every product can be CC BY 4.0** ([#30](https://github.com/traitecoevo/austraits-means/issues/30), [#32](https://github.com/traitecoevo/austraits-means/issues/32)). GBIF licenses a download as its most restrictive record, and both downloads used are CC BY-NC 4.0. In the October 2026 download the CC BY-NC records are 12.3% (3.36 million, mostly iNaturalist research-grade observations, NSW weeds and the ACT Wildlife Atlas); dropping them leaves 603 of 27,730 species names with no records, of which only 49 have a trait mean. Alternatives considered: licensing the occurrence-derived products CC BY-NC, or treating aggregates as outside the NC terms.
- **The occurrence refresh uses GBIF download `0027765-260928105237408`** ([10.15468/dl.766erw](https://doi.org/10.15468/dl.766erw), 15 October 2026), cleaned with the current cornsters_australian_annuals method, whose capital-city buffer is 2 km rather than the 20 km of the 2025 extract (fewer genuine urban records lost). That download uses GBIF's new taxonomic backbone (Catalogue of Life), so its names are realigned to APC rather than reusing the 2025 alignment.

- **GBIF names are aligned to APC from the provider's name (`verbatimScientificName`), with GBIF's interpreted species name as fallback** ([#32](https://github.com/traitecoevo/austraits-means/issues/32), [#33](https://github.com/traitecoevo/austraits-means/issues/33)). GBIF's backbone now follows the Catalogue of Life, which puts some common Australian species in genera APC does not use (*Corymbia maculata* as *Blakella maculata*, *Leptospermum* species as *Gaudium*), so those species vanished from an extract keyed on GBIF's names; it also lumps some species APC keeps apart (*Callitris glaucophylla* into *C. columellaris*). In the October 2026 download the two routes agree for 95.9% of records; the provider name alone aligns for 2.9% and GBIF's name alone for 0.06%, and the species reached rise from 21,708 to 24,370. Infraspecific names are rolled up to species, as on the trait side. Alternative considered: aligning GBIF's interpreted names only, as before.

### 14 October 2026

- **The term is "mean", not "average".** Both species estimates are geometric means, so "mean" is accurate; "average" is ambiguous and clashes with `average_temperature`. Flexibility of method is carried by the `method` column, not by the term. An earlier preference for "average" (to leave room for medians and trimmed means) was considered and set aside.
- *(Superseded 24 November 2026: the portal-method product is retired, unpublished.)* **The portal-method observed mean reproduces the deployed portal output for AusTraits 7.0.0.** The portal computes its species averages offline (`estimate_species_trait_means()` in `austraits.portal`, at the commit that built the deployed data); that is what users have been seeing, so it is the product to archive. Exactness is at key level (taxon × trait × unit) with numeric tolerance; row order is excluded because the portal has no lockfile.
- **That product is an arithmetic mean and is labelled as such.** The geometric observed mean is the successor method from the next AusTraits build onward, with the arithmetic mean kept in its own column. The portal method's quirks (a geometric-mean column that equals the mean; no life-stage or basis-of-record filter on raw/mean/median records; lone flora bounds; no taxonomic roll-up) are documented, not fixed, in that product.
- *(Superseded 24 November 2026: observed means go out with the AusTraits release; modelled means get a new record.)* **Derived products get their own Zenodo deposit**, "AusTraits species means", with its own concept DOI and one version per AusTraits release (or per development build, tagged with its version string), each linked to the release or build it derives from. Alternative considered: re-releasing AusTraits 7.0.0 with the product added. Set aside because a derived product is not a database release and the two would drift.
- **Versioning follows the AusTraits release tag** (see Products).
- **Trait and height records come from an austraits.build development build** (which now includes the AVH height dataset), replacing the separate AVH CSV. The GBIF refresh and medoid rebuild follow in the release candidate; the first draft uses the existing occurrence extract and medoids, with the crosswalk rebuilt from the new names.
- **Unconverged BRMS fits stay postponed** to the extensions milestone; TMB estimates serve the first two milestones.
- **The name crosswalk is built inside the pipeline, against APC release 2026-07-30** (`apc_version` in `_targets.R`), from every taxon name in the development build and every GBIF name ([#13](https://github.com/traitecoevo/austraits-means/issues/13)). Changing the release rebuilds everything downstream. The release candidate re-pins it when the GBIF extract is refreshed ([#33](https://github.com/traitecoevo/austraits-means/issues/33)).
- **The field-only per-record geometric mean (`estimate_simple_means()`) is retired** ([#17](https://github.com/traitecoevo/austraits-means/issues/17)). The geometric observed mean ([#15](https://github.com/traitecoevo/austraits-means/issues/15)) supersedes it as the observed mean, and the combined table (modelled where available, otherwise observed) is built from that; since 24 November 2026 that table is an internal input, not a product. Exports carry `estimate_type` and `method` and name their values `value`, `log_value`, `value_lower`, `value_upper` (`terminology.md`).
- **Infraspecific taxa are rolled up to species** for the modelled means and the geometric observed means (`roll_up_to_species()`; [#13](https://github.com/traitecoevo/austraits-means/issues/13), [#15](https://github.com/traitecoevo/austraits-means/issues/15)). The AVH height dataset in austraits.build keeps subspecies, varieties and forms, which the earlier AVH CSV had rolled up to species, and GBIF names (and so medoids) are almost all at species rank. Rolling up gives every record a species-level climate centre and makes both species estimates join on species. The build's own name is kept in `taxon_name_aligned`. Alternative considered and set aside: keeping infraspecific taxa as their own groups, centred on their parent species' medoid. The portal-method product is not rolled up; that is one of its documented quirks.
- **Weighting and woody / non-woody assignment for the community layers are decisions, not defaults.** The first draft uses provisional choices, flagged as such; the release candidate applies the decisions.
- **The spatial layers get a community consultation**, scheduled after the first-draft review. (Open since 24 November 2026: whether it still applies to statistics predicted from climate.)

### Earlier (refactor plan, 2026)

- Centring on the medoid: `x_c = obs − medoid`, `x_mn = medoid` per climate variable; consistent across traits.
- Model structure: `log_value ~ v1_mn + v2_mn + v1_c + v2_c + (1 | taxon) + (0 + v1_c | taxon) + (0 + v2_c | taxon)`, with v1 = log annual precipitation (WorldClim BIO12) and v2 = mean annual temperature (BIO1). Random slopes are uncorrelated and there is no climate interaction; both choices were driven by convergence. Species need at least 5 records (`n_min = 5`).
- One fitting interface for both backends: `fit_trait_climate_model(backend = "tmb" | "brms")` returns a standard fit object, so downstream code is backend-blind. BRMS uses a non-centred parameterisation. Parallelism is `furrr` multisession, never multicore.
- Per-trait models now (TMB default, BRMS for the full posterior); joint multivariate and gap-filling deferred.
- Species with no medoid are dropped from the climate-response models and routed to the observed means.
- AusTraits records for the modelled means: `basis_of_record == "field"`, `value_type` in raw, mean, median; per-dataset exceptions are allowed where a dataset carries only bounds (the AVH height dataset).
- BRMS validated against TMB (centres r = 1.000); production path set up (cmdstanr, threading, on-disk cache, diagnostics); larger traits need `iter ≥ 4000`.
- Medoid quality gated by a reasonableness check (`medoid_ok`): the medoid sits inside its own occurrence cloud, agrees with simple summaries (median and mean of occurrences), and falls within the climate span of the species' trait records. Species that fail are excluded from the models and get only an observed mean. (Not so in the pipeline, found 24 November 2026: `check_medoid()` is only called by `reports/dataset_construction.qmd`, and no species is excluded by it. Applying the gate is open.)
- Height methodology: one `dataset_id` per AVH collector was carried for a later methodology factor; with the move to the austraits.build dataset this becomes a single `dataset_id`. The build carries no record identifiers for that dataset, so the collector is recoverable only by mapping the identifiers in its `data.csv` into the build ([#3](https://github.com/traitecoevo/austraits-means/issues/3)).
- Orchestration by `targets`, dependencies pinned by `renv`, all logic in `R/` functions, notebooks as reporting layers, figures saved as PNG and PDF.
- Trait list: plant height plus LMA, leaf area, leaf N and P per mass and per area, seed mass, wood density. `leaf_delta13C` deferred.

## Invariants

Things a future edit could break without noticing:

- A community range or community variance is never computed from the pooled species variance of a model without location effects: that includes how locations with the same climate differ, and overstates the spread at one location.
- Community statistics are given per growth form (woody, non-woody); the variance is between species within one form and one location.
- The portal-method code is retired and unmaintained. If it is ever revived, it must reproduce the deployed portal output for 7.0.0, quirks included.
- Value-type rules are per dataset: datasets that carry only minimum/maximum records (the AVH height dataset) enter the modelled means by an explicit rule, not by the default raw/mean/median filter.
- Centring is on the GBIF medoid, never on the trait-sample mean.
- Both species estimates are geometric means from the geometric method onward, except the observed mean of a signed trait, which is arithmetic and labelled `method = arithmetic` (24 November 2026); otherwise an arithmetic mean is always a separately labelled column.
- Column names follow `terminology.md`: `estimate_type` and `method`, never a bare `estimate` or `geom_mean`.
- Community statistics are unweighted unless the weighting decision says otherwise; unweighted values are never called CWM or CWV.
- Use "mean", never "average", in code, exports, docs and metadata.

## Later work

Next, before the trait list and publication:

- **Location model for the community range** ([#41](https://github.com/traitecoevo/austraits-means/issues/41), feeding [#38](https://github.com/traitecoevo/austraits-means/issues/38)). Extend the woodiness model (`fit_woodiness_model()`) with a location random effect and a species × location term, variances per growth form:

  ```
  log_value ~ 0 + g + g:(P_mn + T_mn + P_c + T_c)
            + diag(0 + g | taxon_name) + (0 + P_c | taxon_name) + (0 + T_c | taxon_name)
            + diag(0 + g | location) + (1 | taxon_name:location)
  ```

  `g` is the growth form and `location` a coarse equal-area cell, 50 km by default (coarse so that a location pools several studies; see the result below). The location variance (same climate, different place) is left out of the community range. The species variance, the species × location variance and the spread of co-occurring species' climate centres around the location's climate go into it. The prototype (`R/model/local_range.R`, run by `R/scripts/run_location_prototype.R`, reported in `docs/local-range-model.qmd`) is kept out of the product targets until it passes all four checks:
  1. A simulation in which species sort among locations: the model recovers the known community SD, and the current model overstates it.
  2. The model's community SD agrees with the observed within-location SD, per growth form.
  3. Likelihood and AIC with and without the location terms, and the convergence and fit time for height (~320k records, ~29k cells).
  4. The shift in the modelled species means, old against new.

  Risk to test, not assume: one location term covers both species sorting among locations and species responding to a location; if sorting dominates, the split between species and location variance is biased.

  Evidence (24 November 2026, `austraits_climate` target, natural-log SDs). "Within a location" is the SD among species' means at sites with 5 or more species. "Pooled" is the SD of species means of residuals from a climate-only regression, a crude benchmark not yet split by growth form:

  | trait | locations with ≥5 species | within a location | pooled, given climate | ratio |
  |---|---|---|---|---|
  | leaf_mass_per_area | 459 sites | 0.45 | 0.76 | 0.59 |
  | leaf_P_per_dry_mass | 178 sites | 0.47 | 0.70 | 0.67 |
  | leaf_N_per_dry_mass | 259 sites | 0.37 | 0.49 | 0.76 |
  | leaf_area | 324 sites | 1.45 | 1.86 | 0.78 |
  | leaf_P_per_area | 72 sites | 0.41 | 0.51 | 0.80 |
  | wood_density | 120 sites | 0.21 | 0.25 | 0.84 |
  | leaf_N_per_area | 101 sites | 0.40 | 0.45 | 0.89 |
  | plant_height | 10,924 cells (10 km) | 1.19 | 1.30 | 0.92 |
  | seed_dry_mass | 53 sites | 2.19 | 2.40 | 0.91 |

  Site, exact-coordinate and 10 km-cell groupings give near-identical SDs for the eight leaf, seed and wood traits. Height is 99% herbarium specimens with one record per location, so it has replication only at the cell grain.

  Prototype result, 24 November 2026, all nine traits at 50 km (`docs/local-range-model.qmd`):

  - **Checks 1, 3 and 4 pass.**
    - Simulation: the model recovers known variances. With independent placement it gets the local spread within about 8%; the published model overstates it by 1.4–1.7×.
    - Every trait converges, and AIC favours the location model for every trait.
    - The modelled species means move by 0.06–0.35 log units (SD), correlation 0.85–1.00 with the old ones.
  - **Check 2 mostly passes.** Model-to-observed community SD ratios:
    - location model: 0.97–1.31, mostly 1.05–1.20, against 1.00–1.97 for the published model;
    - the exception is wood density at 1.47.
  - **A 10–20% overstatement remains, and its cause is open.** Either co-occurring species really are more alike than random (the true range is smaller), or studies co-sample similar species (the observed spread is too small). These data cannot separate the two.
  - **Coarse cells cost nothing.** Cell size made little difference to LMA between 10, 25 and 50 km.
  - **Daniel's concern about sites with few species.** The range is estimated from variance parameters across all the data, not from each site's own spread. A site with one or two species adds little information but does not shrink the range. Selective sampling of similar species does bias it low, and coarse cells reduce that.

  Next:
  1. Predict to a location from its climate alone. This needs the mean and spread of co-occurring species' climate centres as a smooth function of climate, from GBIF co-occurrence (extend `cell_species_centres()`).
  2. Choose the trait list. Benchmark coverage of both growth forms is good for height, LMA and leaf area. The other six traits have woody cells only, or fewer than 10 non-woody cells.
  3. Decide whether a genus effect should take up the similarity of co-sampled congeners.

On the extensions milestone:
- Collection-method covariate ([#3](https://github.com/traitecoevo/austraits-means/issues/3), [#36](https://github.com/traitecoevo/austraits-means/issues/36)): waits for the factor tables and a test on taxa measured under more than one method.
- Joint multivariate model with gap-filling ([#40](https://github.com/traitecoevo/austraits-means/issues/40)).

Deferred, no milestone:

- `leaf_delta13C` ([#39](https://github.com/traitecoevo/austraits-means/issues/39)): needs the discrimination transform.
