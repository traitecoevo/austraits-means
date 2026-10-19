# Community layers: design

The community layers summarise trait values of the plant species in each cell of a 10 × 10 km equal-area grid of Australia, for all species and for woody and non-woody species separately. Per trait, cell and group they give the **community mean**, exp(mean of the species' log values), a geometric mean in trait units, and for the observed products the **community variance**, the variance of the species' log values, in natural-log units (`terminology.md`). The community variance is between species only: it leaves out variation within species.

Three products are built, which differ in where the species and their values come from:

| Product | Species in a cell | Value of each species | Model |
|---|---|---|---|
| `observed_local` | species with a trait record collected in the cell | geometric mean of its records in that cell | none |
| `observed_presence` | species with a GBIF occurrence in the cell | its observed mean (geometric), wherever it was measured | none |
| `modelled` | species with a GBIF occurrence in the cell | none: the trait-climate model's net average over those species | the per-trait model |

Every output carries a metadata file stating its construction and data basis.

## Decisions

All decided 19 October 2026 (`roadmap.md`).

| Choice | Decision | Issue |
|---|---|---|
| Weighting of species | presence only: every species in a cell counts equally; hence "community mean" and "community variance", never CWM or CWV | [#18](https://github.com/traitecoevo/austraits-means/issues/18) |
| Woody / non-woody | the species' most likely `woodiness_detailed` value: `woody` and `woody_like_stem` (palms, cycads, grass trees) are woody; `herbaceous`, `semi_woody`, `woody_base`, `woody_like_inflorescence` and `woody_root` are non-woody (the intermediate ones flagged); otherwise growth-form keywords; otherwise unassigned. Unassigned species are left out of the woody and non-woody layers and count only in `all`; the coarse `woodiness` trait is not used as a fallback | [#19](https://github.com/traitecoevo/austraits-means/issues/19) |
| Spatial unit | 10 × 10 km cells of a Lambert Azimuthal Equal Area grid centred on Australia, as in Andrew et al. 2025 (below) | [#20](https://github.com/traitecoevo/austraits-means/issues/20) |
| Minimum species | statistics only where at least 10 species of the group in the cell have a value; counts are always given | [#20](https://github.com/traitecoevo/austraits-means/issues/20) |
| Occurrence source | GBIF occurrence points: a species is present in a cell if it has any cleaned occurrence there | [#20](https://github.com/traitecoevo/austraits-means/issues/20) |
| Species values | `observed_local`: the species' records in the cell; `observed_presence`: the geometric observed means; `modelled`: no species values | [#20](https://github.com/traitecoevo/austraits-means/issues/20), [#28](https://github.com/traitecoevo/austraits-means/issues/28) |

Species are per species, with subspecies, varieties and forms rolled up, and GBIF names are matched to AusTraits names through the name crosswalk aligned to APC release 2026-07-30 (`roadmap.md`).

## The grid

`community_layer_settings()` in `R/community/community_layers.R`: projection `+proj=laea +lat_0=-25.2744 +lon_0=133.7751 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs`, 10,000 m cells, extent x −2,300,000 to 2,200,000 m and y −2,200,000 to 1,800,000 m (450 × 400 cells). Cell (`cell_x`, `cell_y`) spans [`cell_x`, `cell_x` + 1) × 10 km east and likewise north; outputs give the cell centre as `x`, `y` (metres) and `longitude`, `latitude`.

## Method

### Observed, local records (`observed_local`)

`community_layers_local()`. The field records of the modelled traits (`records_all`: value types raw, mean and median, plus the AVH height bounds, infraspecific taxa rolled up) are placed in cells. In each cell, each species' records are averaged on the log scale, then the species are averaged. No model and no occurrence data are used, so only cells where the trait was measured have values, and the value describes the species measured there, at that place.

### Observed, GBIF presence × species means (`observed_presence`)

`gbif_presence_cells()` places each cleaned GBIF occurrence in its cell (one row per GBIF name × cell). `gbif_to_species()` maps each GBIF name to one AusTraits species, directly or through the crosswalk. `community_layers_presence()` then gives each species present its geometric observed mean and averages them. `n_species_present` counts all GBIF species in the cell, including those with no AusTraits name; `n_species` those with an observed mean of the trait; `coverage` their ratio.

### Modelled (`modelled`)

`R/community/community_modelled.R`. Each trait's model (`log_value ~ P_mn + T_mn + P_c + T_c` with species random intercepts and within-species slopes; P is log annual precipitation, T mean annual temperature, `_mn` the species' climate centre (medoid) and `_c` the record's climate minus that centre) gives a record at climate x from a species centred at m the expected log value

b0 + b_between · m + b_within · (x − m) + species effects.

The species effects average to zero over species, so the average over the species in a cell is

b0 + b_within · x + (b_between − b_within) · m̄,

where x is the cell's climate (WorldClim averaged onto the grid, `cell_climate()`) and m̄ the mean climate centre of the GBIF species recorded in the cell (`cell_species_centres()`), for each climate variable. Back-transformed, this is the model's net average trait value for the cell: it averages over between- and within-species effects, uses no species means, and needs no trait records from the species in the cell, only their climate centres. The 94% interval (`lower`, `upper`) comes from the fixed-effect uncertainty (delta method on the log scale). A cell is flagged `extrapolated` where m̄ or x − m̄ lies outside the range of the species centres or within-species deviations the model was fitted on. There is no community variance: the model's between-species variance is the same in every cell.

The per-trait models are fitted on all species together, so this product is for `all` species only. A woody / non-woody split needs either separate models per group or woodiness in the model; that choice is open.

## Outputs

In `outputs/community/<build version>/`:

- `woody_lookup.csv`, `woody_lookup_summary.csv`: the group of each taxon and how it was assigned, shared by all products.
- `<product>/community_layers.parquet`: one row per trait × cell × group: `trait_name`, `cell_x`, `cell_y`, `x`, `y`, `longitude`, `latitude`, `woodiness_group` (`woody`, `non_woody`, `all`), the counts, and the statistics. `observed_local` adds `n_records`; `observed_presence` `n_species_present` and `coverage`; `modelled` the cell climate, the mean species centre, `lower`, `upper` and `extrapolated`.
- `<product>/community_layers_metadata.json`: definitions, construction, decisions, and the data basis (build version and commit, APC release, GBIF download, medoid build).
- `<product>/geotiff/<trait>_<group>_<statistic>.tif`: one single-band GeoTIFF per trait × group × statistic on the grid.
- `<product>/figures/`: quick-look maps per trait. For the observed products, community mean and variance, woody and non-woody (`<trait>_observed_local`, `<trait>_observed_presence`); for `modelled`, the modelled community mean beside the `observed_presence` one for all species, on one colour scale (`<trait>_modelled`). Grey cells have species but fewer than 10 with a value, or (modelled) are extrapolated.
- `comparison/`: the products against each other, cell by cell, for all species (`R/community/compare_layers.R`): `community_layer_agreement.csv` (per pair and trait: cells, correlation of the log community means, median and SD of their difference, and each product's spatial SD), a scatter figure per pair (`figures/pair_presence_modelled`, `pair_local_presence`, `pair_local_modelled`), and maps of modelled minus presence per trait (`figures/difference_modelled_minus_presence`).

## Coverage

Development build 7.0.0.9000 and GBIF download `0027765-260928105237408`. GBIF species occur in 55,176 cells, 33,879 of them with at least 10 species. Cells with statistics per product and group; in brackets, for `observed_presence` the median coverage of the species present by an observed mean, and for `modelled` the share of cells flagged `extrapolated`:

| trait | local: all | local: woody | local: non-woody | presence: all | presence: woody | presence: non-woody | modelled: all |
|---|---|---|---|---|---|---|---|
| plant_height | 5,781 | 3,681 | 1,940 | 33,730 (0.99) | 25,993 (1.00) | 27,341 (1.00) | 33,574 (0%) |
| seed_dry_mass | 30 | 21 | 17 | 31,008 (0.78) | 23,377 (0.82) | 24,704 (0.77) | 33,574 (25%) |
| leaf_area | 160 | 126 | 33 | 26,617 (0.61) | 19,355 (0.69) | 20,831 (0.61) | 33,574 (0.1%) |
| leaf_mass_per_area | 208 | 158 | 45 | 25,100 (0.54) | 20,283 (0.69) | 16,913 (0.52) | 33,574 (0%) |
| leaf_N_per_dry_mass | 100 | 80 | 12 | 21,871 (0.27) | 17,519 (0.47) | 11,385 (0.16) | 33,574 (0.2%) |
| leaf_P_per_dry_mass | 67 | 56 | 16 | 16,255 (0.16) | 10,989 (0.28) | 7,782 (0.10) | 33,574 (0.8%) |
| wood_density | 48 | 48 | 1 | 15,817 (0.19) | 15,685 (0.52) | 402 (0.05) | 33,574 (7.2%) |
| leaf_N_per_area | 50 | 41 | 2 | 17,008 (0.16) | 14,184 (0.34) | 4,543 (0.06) | 33,574 (0.6%) |
| leaf_P_per_area | 34 | 29 | 1 | 11,398 (0.10) | 9,754 (0.24) | 198 (0.02) | 33,574 (0.9%) |

The modelled product has a value wherever at least 10 GBIF species with a climate centre occur and the cell has WorldClim climate, whatever the trait.

How the products compare, all species, cells with both (modelled cells not extrapolated): correlation of log community means, and median difference in log units (from `comparison/community_layer_agreement.csv`):

| trait | modelled vs presence | modelled − presence | presence vs local (cells) |
|---|---|---|---|
| plant_height | 0.54 | 0.24 | 0.67 (5,707) |
| seed_dry_mass | 0.48 | 0.99 | 0.84 (30) |
| leaf_area | 0.78 | 0.72 | 0.83 (157) |
| leaf_mass_per_area | 0.25 | 0.24 | 0.47 (203) |
| leaf_N_per_dry_mass | −0.21 | −0.21 | 0.70 (94) |
| leaf_P_per_dry_mass | −0.27 | −0.38 | 0.77 (61) |
| wood_density | 0.86 | 0.01 | 0.72 (45) |
| leaf_N_per_area | 0.88 | 0.10 | 0.92 (48) |
| leaf_P_per_area | 0.66 | 0.04 | 0.76 (32) |

Two things drive the differences between the modelled and presence products. The modelled level is set by the species the model was fitted on, which for most traits are a small, well-studied subset (773 species for leaf area, 214 for seed mass); for leaf area and plant height they average 0.66 and 0.26 log units above all species with an observed mean. And the modelled pattern comes only from the two climate variables, so shifts in composition that climate does not predict are missing: for leaf N and P per mass, the presence product is highest across the south-east, while the model's climate slopes are weak and its surface nearly flat (spatial SD of the log community mean 0.05 for leaf N per mass, against 0.19 for presence). The modelled community mean varies less across cells than the presence one for every trait except leaf area.

## Precedent: Andrew et al. 2025

Andrew, Martín-Forés, Guerin, Coleman, Falster, Wenk, Wright and Gallagher (2025), "Mapping plant functional traits using gap-filled datasets to inform management and modelling", *Global Ecology and Biogeography* 34: e70136 ([doi:10.1111/geb.70136](https://doi.org/10.1111/geb.70136)), mapped the community mean and variability of maximum height, seed mass, leaf area and LMA across Australia from AusTraits species means. It is the closest precedent for these layers:

| Choice | Andrew et al. 2025 | These layers |
|---|---|---|
| Presence | binary range maps from species distribution models (Poisson point process models, range bagging or area of occupancy) fitted to vouchered ALA specimens collected since 1950, limited to bioregions where the species occurs | GBIF occurrence points |
| Taxa | native angiosperms and gymnosperms only (17,251 species) | all species with a GBIF record, native or naturalised |
| Species trait values | AusTraits species means, gap-filled with BHPMF and Rphylopars to 95–100% coverage; gap-filled and raw maps agree (R² 0.8–0.99) except mean LMA (0.37) | observed means, no gap-filling (`observed_presence`); records in the cell (`observed_local`); none (`modelled`) |
| Weighting | presence only, every species equal | presence only |
| Spatial unit | 10 × 10 km equal-area grid (Lambert Azimuthal Equal Area) | the same |
| Statistics | mean of log trait values; variability as the log of their SD | community mean and community variance of log values |
| Groups | all species together | woody and non-woody, plus all |
| Validation | community means and variability from 747 AusPlots plot inventories, compared on their trait-climate slopes | not yet |

Its maps, gap-filled traits and species range maps are published with the paper, so its maps are a direct check on these layers and its range maps an alternative occurrence source.

## Known limitations

- Presence comes from occurrence points, so it follows collecting effort: well-collected cells have more species, and much of the interior has fewer than 10.
- In `observed_presence`, where coverage is low the statistic describes the well-studied fraction of the community, which is unlikely to be a random sample of it; and each species carries one value from wherever it was measured, not the value it would have in that cell.
- `observed_local` covers few cells for most traits, and its species in a cell are those someone chose to measure there.
- `modelled` depends on the model's two climate variables only, and on the species with enough records and a medoid to fit it; cells far from that climate space are flagged.
- Wood density for non-woody plants is nearly empty and of doubtful meaning; the layer is kept for completeness.
