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
| Species-means input | two variants are built: `combined` (the modelled mean, TMB, where there is one, otherwise the geometric observed mean) and `modelled` (modelled means only); recorded in metadata, not in the name | [#20](https://github.com/traitecoevo/austraits-means/issues/20) | one of the two; observed means only |

Two further choices follow decisions already taken (`roadmap.md`): species means are per species, with subspecies, varieties and forms rolled up, and GBIF names are matched to them through the name crosswalk aligned to APC release 2026-07-30.

## Method

`R/community/community_layers.R`, run by the pipeline:

1. `combine_species_means()`: for each of the nine modelled traits, the modelled mean where there is one, otherwise the geometric observed mean (`species_means_combined` target). `select_species_means()` gives each variant its input: all of it (`combined`) or the modelled means only (`modelled`).
2. `gbif_presence_cells()`: each GBIF occurrence is assigned to its 0.25° cell; one row per GBIF name × cell (`gbif_presence`; 2.06 million rows over 11,143 cells).
3. `gbif_to_species()`: each GBIF name is mapped to one species-means name, directly or through the crosswalk with its names rolled up to species (22,130 of 22,136 GBIF names map; 22,115 to a taxon with a woody or non-woody class).
4. `compute_community_layers()`: species present in a cell are split by group; for each trait, the species with a mean give `n_species`, the community mean and the community variance, and `coverage` = `n_species / n_species_present`. Below 10 species the statistics are empty but the counts are kept. It runs once per variant (targets `community_layers_combined`, `community_layers_modelled`).

## Outputs

In `outputs/community/<build version>/<variant>/`, where the variant is `combined` or `modelled`:

- `community_layers.parquet`: one row per trait × cell × group with species present: `trait_name`, `cell_x`, `cell_y`, `longitude`, `latitude` (cell centre), `woodiness_group`, `n_species_present`, `n_species`, `coverage`, `community_mean`, `community_variance`.
- `community_layers_metadata.json`: definitions, the provisional choices, the species-means input, the occurrence source, the spatial unit, and the data basis (build version and commit, APC release, GBIF extract, medoid build).
- `geotiff/<trait>_<group>_<statistic>.tif`: one single-band GeoTIFF per trait × group × statistic (`community_mean`, `community_variance`, `n_species`, `coverage`) on the 0.25° grid, EPSG:4326.
- `figures/<trait>_community.png` and `.pdf`: quick-look maps of community mean and variance, woody and non-woody. Grey cells have species present but fewer than 10 with a mean.
- `../woody_lookup_provisional.csv`: the group of each taxon (#22), shared by both variants.

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

## Known limitations of the draft

- Where coverage is low, the statistic describes the well-studied fraction of the community, which is unlikely to be a random sample of it.
- A cell mixes modelled and observed species means; their difference is examined in #24.
- Presence comes from raw occurrence points, so it follows collecting effort: well-collected cells have more species present, and the interior has many cells below the threshold.
- Wood density for non-woody plants is nearly empty and of doubtful meaning; the layer is kept for completeness.
- The GBIF extract dates from September 2025; it is refreshed in the release candidate (#32).
