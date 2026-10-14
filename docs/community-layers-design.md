# Community layers: design of the first draft

The community layers summarise the species means of the plant species present in each spatial unit across Australia, separately for woody and non-woody plants. Per trait, unit and group they give the **community mean**, exp(mean of the log species means), and the **community variance**, the variance of the log species means, in natural-log units (`terminology.md`). The community variance is between species only: it leaves out variation within species.

This first draft is built on **provisional choices**. Each is open in a decision issue, and the release candidate rebuilds the layers on the decisions (#27, #28). Modelled community layers across space are later work (#38). Every output carries a metadata file stating these choices.

## Provisional choices

| Choice | First draft | Decided in | Alternatives |
|---|---|---|---|
| Weighting of species | presence only: every species present counts equally; hence "community mean" and "community variance", never CWM or CWV | [#18](https://github.com/traitecoevo/austraits-means/issues/18) | abundance- or cover-weighted from survey plots (CWM, CWV); both |
| Woody / non-woody | the provisional lookup of #22: modal `woodiness_detailed`, intermediate categories to non-woody (flagged), then growth-form keywords; unassigned taxa left out | [#19](https://github.com/traitecoevo/austraits-means/issues/19) | other placements of intermediate categories; `woodiness` as a fallback; a third class |
| Spatial unit | regular 0.25° grid in EPSG:4326 over 112–154° E, 44–9° S (168 × 140 cells) | [#20](https://github.com/traitecoevo/austraits-means/issues/20) | 0.1° grid; bioregions; survey plots |
| Minimum species | statistics only where at least 10 species of the group present in the cell have a species mean | [#20](https://github.com/traitecoevo/austraits-means/issues/20) | other thresholds; a coverage threshold |
| Occurrence source | the cleaned GBIF occurrences (download `0001494-250914085247600`); a species is present in a cell if it has any occurrence there | [#20](https://github.com/traitecoevo/austraits-means/issues/20) | AVH; both; range estimates (alpha hulls) instead of points |
| Species-means input | three variants are built for review: `combined` (the modelled mean, TMB, where there is one, otherwise the geometric observed mean), `modelled_centre` (modelled means only) and `model_surface` (the hierarchical model's own prediction, below); recorded in metadata, not in the name | [#20](https://github.com/traitecoevo/austraits-means/issues/20) | observed means only |

Two further choices follow decisions already taken (`roadmap.md`): species means are per species, with subspecies, varieties and forms rolled up, and GBIF names are matched to them through the name crosswalk aligned to APC release 2026-07-30.

## Method

`R/community/community_layers.R`, run by the pipeline:

1. `combine_species_means()`: for each of the nine modelled traits, the modelled mean where there is one, otherwise the geometric observed mean (`species_means_combined` target). `select_species_means()` gives each variant its input: all of it (`combined`) or the modelled means only (`modelled_centre`).
2. `gbif_presence_cells()`: each GBIF occurrence is assigned to its 0.25° cell; one row per GBIF name × cell (`gbif_presence`; 2.06 million rows over 11,143 cells).
3. `gbif_to_species()`: each GBIF name is mapped to one species-means name, directly or through the crosswalk with its names rolled up to species (22,130 of 22,136 GBIF names map; 22,115 to a taxon with a woody or non-woody class).
4. `compute_community_layers()`: species present in a cell are split by group; for each trait, the species with a mean give `n_species`, the community mean and the community variance, and `coverage` = `n_species / n_species_present`. Below 10 species the statistics are empty but the counts are kept. It runs once per variant (targets `community_layers_combined`, `community_layers_modelled_centre`). Besides `woody` and `non_woody`, a third group, `all`, pools every species present, including those the lookup leaves unassigned.

## Outputs

In `outputs/community/<build version>/<variant>/`, where the variant is `combined` or `modelled_centre`:

- `community_layers.parquet`: one row per trait × cell × group with species present: `trait_name`, `cell_x`, `cell_y`, `longitude`, `latitude` (cell centre), `woodiness_group`, `n_species_present`, `n_species`, `coverage`, `community_mean`, `community_variance`.
- `community_layers_metadata.json`: definitions, the provisional choices, the species-means input, the occurrence source, the spatial unit, and the data basis (build version and commit, APC release, GBIF extract, medoid build).
- `geotiff/<trait>_<group>_<statistic>.tif`: one single-band GeoTIFF per trait × group × statistic (`community_mean`, `community_variance`, `n_species`, `coverage`) on the 0.25° grid, EPSG:4326.
- `figures/<trait>_community.png` and `.pdf`: quick-look maps of community mean and variance, woody and non-woody. Grey cells have species present but fewer than 10 with a mean.
- `../woody_lookup_provisional.csv`: the group of each taxon (#22), shared by both variants.

## The model surface

A third, model-based variant (`R/community/model_surface.R`, folder `model_surface/`). For each trait, the community mean of a cell is the population-level prediction of the trait's hierarchical model at the cell's mean climate, exp(b0 + b_P log P + b_T T): the expected trait of a species whose climate centre is that climate, with species effects and within-species terms set to zero, so averaged over species. It covers every land cell, needs no occurrence data, and has a 94% interval from the fixed-effect uncertainty. Under the model the community variance is the between-species intercept variance, the same in every cell; it is recorded in the metadata (`sd_species`), not mapped. The per-trait models are not split into woody and non-woody, so the surface is for all species. Cells whose climate lies outside the range of the species centres in the fit are flagged `extrapolated` (grey on the maps).

Outputs: `model_surface.parquet` (trait × cell: climate, `community_mean`, `lower`, `upper`, `extrapolated`, `sd_species`), `geotiff/<trait>_all_community_mean_modelled.tif`, metadata, and maps placing the surface beside the observed community mean of all species (combined input) on one colour scale.

Against the observed community mean of all species (cells with both), the surface agrees well for wood density and leaf N per area (log correlation 0.87 and 0.90, median offset near zero), moderately for leaf area, seed mass and plant height (0.83, 0.62, 0.52) but with the surface higher by a median 0.77, 1.25 and 0.18 log units, and poorly for leaf N and P per mass (correlation −0.28 and −0.19). The surface rests on the species with a modelled mean, which are few and well studied for most traits (212 for seed mass), and their climate centres cover only part of the continent: 55% of cells are extrapolated for seed mass and 38% for wood density, against 1.4% for plant height.

## Coverage in this draft

Cells with statistics (of about 10,700 cells per group with any species present), and in brackets the median coverage of present species by a species mean, for each variant:

| trait | woody, combined | woody, modelled | non-woody, combined | non-woody, modelled |
|---|---|---|---|---|
| plant_height | 8,317 (1.00) | 8,158 (0.93) | 8,382 (1.00) | 8,092 (0.86) |
| seed_dry_mass | 7,848 (0.78) | 1,249 (0.02) | 7,910 (0.73) | 601 (0.00) |
| leaf_area | 6,711 (0.50) | 4,043 (0.16) | 6,992 (0.48) | 2,218 (0.04) |
| leaf_mass_per_area | 6,951 (0.53) | 4,905 (0.22) | 5,742 (0.30) | 3,364 (0.07) |
| leaf_N_per_dry_mass | 6,526 (0.39) | 3,846 (0.13) | 4,668 (0.12) | 217 (0.02) |
| leaf_P_per_dry_mass | 4,552 (0.19) | 1,707 (0.06) | 3,192 (0.07) | 190 (0.02) |
| wood_density | 5,863 (0.35) | 2,209 (0.07) | 244 (0.00) | 0 |
| leaf_N_per_area | 5,555 (0.25) | 1,589 (0.06) | 1,927 (0.03) | 0 |
| leaf_P_per_area | 3,555 (0.12) | 828 (0.04) | 115 (0.00) | 0 |

The modelled-only variant is close to the combined one for plant height, where most species present have enough records and a medoid, and much sparser for the other traits, where a modelled mean needs more than 5 field records per species.

## Precedent: Andrew et al. 2025

Andrew, Martín-Forés, Guerin, Coleman, Falster, Wenk, Wright and Gallagher (2025), "Mapping plant functional traits using gap-filled datasets to inform management and modelling", *Global Ecology and Biogeography* 34: e70136 ([doi:10.1111/geb.70136](https://doi.org/10.1111/geb.70136)), mapped the community mean and variability of maximum height, seed mass, leaf area and LMA across Australia from AusTraits species means. It is the closest precedent for these layers, and each of its choices bears on an open decision:

| Choice | Andrew et al. 2025 | This draft | Bears on |
|---|---|---|---|
| Presence | binary range maps from species distribution models (Poisson point process models, range bagging or area of occupancy) fitted to vouchered ALA specimens collected since 1950, limited to bioregions where the species occurs | raw GBIF occurrence points | #20: range estimates instead of points would remove the collecting-effort gaps in the interior |
| Taxa | native angiosperms and gymnosperms only (17,251 species) | all species with a GBIF record, native or naturalised | #20 |
| Species trait values | AusTraits species means, gap-filled with BHPMF and Rphylopars to 95–100% coverage; gap-filled and raw maps agree (R² 0.8–0.99) except mean LMA (0.37) | combined species means, no gap-filling | coverage of the leaf traits; gap-filling is #40 |
| Weighting | presence only, every species equal | presence only | #18 |
| Spatial unit | 10 × 10 km equal-area grid (Lambert Azimuthal Equal Area) | 0.25° grid, EPSG:4326 | #20 |
| Statistics | mean of log trait values; variability as the log of their SD | community mean and community variance of log species means | naming in `terminology.md` |
| Groups | all species together | woody and non-woody, plus all | #19 |
| Validation | community means and variability from 747 AusPlots plot inventories, compared on their trait-climate slopes | not yet | review of these layers; plot data also carry the abundances #18 would need |

Its maps, gap-filled traits and species range maps are published with the paper, so the range maps are a candidate occurrence source for #20 and its maps a direct check on these layers.

## Known limitations of the draft

- Where coverage is low, the statistic describes the well-studied fraction of the community, which is unlikely to be a random sample of it.
- A cell mixes modelled and observed species means; their difference is examined in #24.
- Presence comes from raw occurrence points, so it follows collecting effort: well-collected cells have more species present, and the interior has many cells below the threshold.
- Wood density for non-woody plants is nearly empty and of doubtful meaning; the layer is kept for completeness.
- The GBIF extract dates from September 2025; it is refreshed in the release candidate (#32).
