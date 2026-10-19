# Community layers: design

The community layers summarise trait values of the plant species in each cell of a 10 × 10 km equal-area grid of Australia, for woody and non-woody species separately, and for all species together. Per trait, cell and group they give the **community mean**, exp(mean of the species' log values), a geometric mean in trait units, and for the observed products the **community variance**, the variance of the species' log values, in natural-log units (`terminology.md`). The community variance is between species only: it leaves out variation within species.

Three products are built, which differ in where the species and their values come from:

| Product | Species in a cell | Value of each species | Model |
|---|---|---|---|
| `observed_local` | species with a trait record collected in the cell | geometric mean of its records in that cell | none |
| `observed_presence` | species with a GBIF occurrence in the cell | its observed mean (geometric), wherever it was measured | none |
| `modelled` | species with a GBIF occurrence in the cell | none: the trait-climate model's net average over those species | each trait's model, with woody and non-woody coefficients |

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

`R/community/community_modelled.R`. Each trait gets its own model for the layers (`fit_woodiness_model()`, targets `woody_fit_<trait>`): the trait-climate model of the modelled species means, with a separate intercept and climate slopes for woody and non-woody species and the species random effects shared,

`log_value ~ 0 + group + group:(P_mn + T_mn + P_c + T_c) + (1 | taxon) + (0 + P_c | taxon) + (0 + T_c | taxon)`,

where P is log annual precipitation, T mean annual temperature, `_mn` the species' climate centre (medoid) and `_c` the record's climate minus that centre. It is fitted on the species with more than 5 records, a medoid and a woody or non-woody class. For a group g, a record at climate x from a species centred at m has expected log value

b0_g + b_between,g · m + b_within,g · (x − m) + species effects.

The species effects average to zero over species, so the average over the group's species in a cell is

b0_g + b_within,g · x + (b_between,g − b_within,g) · m̄_g,

where x is the cell's climate (WorldClim averaged onto the grid, `cell_climate()`) and m̄_g the mean climate centre of the group's GBIF species recorded in the cell (`cell_species_centres()`), for each climate variable. Back-transformed, this is the model's net average trait value for the group in that cell: it averages over between- and within-species effects, uses no species means, and needs no trait records from the species in the cell, only their climate centres. The value for all species is the woody and non-woody values mixed on the log scale in proportion to their species in the cell, which is the same average taken over the species of both groups. Species without a woody or non-woody class, and GBIF species with no AusTraits name, do not enter the modelled layers. The 94% interval (`lower`, `upper`) comes from the fixed-effect uncertainty (delta method on the log scale). A cell is flagged `extrapolated` where m̄_g or x − m̄_g lies outside the range of the group's species centres or within-species deviations the model was fitted on. There is no community variance: the model's between-species variance is the same in every cell.

A group is fitted only if it has at least 20 species with enough records. Non-woody falls short for leaf N per area (8 species), leaf P per area (7) and wood density (2); for these traits the model is fitted on woody species alone, and the non-woody and all-species modelled layers are empty. All nine fits converge (positive-definite Hessian). Seed mass (27 non-woody species), leaf P per mass (37) and leaf N per mass (42) have few non-woody species, so their non-woody coefficients are uncertain and many cells are extrapolated (below). Coefficients per trait and group are in `modelled/woodiness_model_fixed_effects.csv`.

Why a woodiness term rather than one model for all species: the first version of this product used the pooled model of the species means, which is fitted mostly on woody species (81–99% of fitted species for most traits), and could not follow differences in the woody share of the species in a cell (below). Separate woody and non-woody models were the alternative; one model with group coefficients shares the random-effect variances, which helps the groups with few species.

## Outputs

In `outputs/community/<build version>/`:

- `woody_lookup.csv`, `woody_lookup_summary.csv`: the group of each taxon and how it was assigned, shared by all products.
- `<product>/community_layers.parquet`: one row per trait × cell × group: `trait_name`, `cell_x`, `cell_y`, `x`, `y`, `longitude`, `latitude`, `woodiness_group` (`woody`, `non_woody`, `all`), the counts, and the statistics. `observed_local` adds `n_records`; `observed_presence` `n_species_present` and `coverage`; `modelled` the cell climate, the group's mean species centre, `lower`, `upper` and `extrapolated`, and `modelled/woodiness_model_fixed_effects.csv` the model coefficients.
- `<product>/community_layers_metadata.json`: definitions, construction, decisions, and the data basis (build version and commit, APC release, GBIF download, medoid build).
- `<product>/geotiff/<trait>_<group>_<statistic>.tif`: one single-band GeoTIFF per trait × group × statistic on the grid.
- `<product>/figures/`: quick-look maps per trait. For the observed products, community mean and variance, woody and non-woody (`<trait>_observed_local`, `<trait>_observed_presence`); for `modelled`, the modelled community mean beside the `observed_presence` one, woody and non-woody, on one colour scale (`<trait>_modelled`). Grey cells have species but fewer than 10 with a value, or (modelled) are extrapolated or in a group the model could not fit.
- `comparison/`: the products against each other, cell by cell, for woody, non-woody and all species (`R/community/compare_layers.R`): `community_layer_agreement.csv` (per pair, group and trait: cells, correlation of the log community means, median and SD of their difference, and each product's spread across cells), a scatter figure per pair and group (`figures/pair_<pair>_<group>`, pairs `presence_modelled`, `local_presence`, `local_modelled`), and maps of modelled minus presence per group (`figures/difference_presence_modelled_<group>`).
- Native-only test (below): `observed_presence_native/` and `modelled_native/` in the same layout, `comparison_native/` comparing them as above, and `comparison_native_vs_all/` comparing each with its version that includes introduced species (scatter figures and difference maps per group).

## Coverage

Development build 7.0.0.9000 and GBIF download `0027765-260928105237408`. GBIF species occur in 55,176 cells, 33,879 of them with at least 10 species.

Observed products, cells with statistics per group; in brackets, the median coverage of the species present by an observed mean:

| trait | local: woody | local: non-woody | local: all | presence: woody | presence: non-woody | presence: all |
|---|---|---|---|---|---|---|
| plant_height | 3,681 | 1,940 | 5,781 | 25,993 (1.00) | 27,341 (1.00) | 33,730 (0.99) |
| seed_dry_mass | 21 | 17 | 30 | 23,377 (0.82) | 24,704 (0.77) | 31,008 (0.78) |
| leaf_area | 126 | 33 | 160 | 19,355 (0.69) | 20,831 (0.61) | 26,617 (0.61) |
| leaf_mass_per_area | 158 | 45 | 208 | 20,283 (0.69) | 16,913 (0.52) | 25,100 (0.54) |
| leaf_N_per_dry_mass | 80 | 12 | 100 | 17,519 (0.47) | 11,385 (0.16) | 21,871 (0.27) |
| leaf_P_per_dry_mass | 56 | 16 | 67 | 10,989 (0.28) | 7,782 (0.10) | 16,255 (0.16) |
| wood_density | 48 | 1 | 48 | 15,685 (0.52) | 402 (0.05) | 15,817 (0.19) |
| leaf_N_per_area | 41 | 2 | 50 | 14,184 (0.34) | 4,543 (0.06) | 17,008 (0.16) |
| leaf_P_per_area | 29 | 1 | 34 | 9,754 (0.24) | 198 (0.02) | 11,398 (0.10) |

Modelled product: a group has a value wherever at least 10 of its GBIF species with a climate centre occur and the cell has WorldClim climate, whatever the trait: 25,891 cells for woody, 27,284 for non-woody, 33,565 for all species. In brackets, the share of those cells flagged `extrapolated`; "not fitted" marks a group with too few species:

| trait | woody | non-woody | all |
|---|---|---|---|
| plant_height | 0% | 0% | 0% |
| seed_dry_mass | 23% | 47% | 52% |
| leaf_area | 0.1% | 20% | 22% |
| leaf_mass_per_area | 0% | 17% | 19% |
| leaf_N_per_dry_mass | 0.3% | 8.2% | 7.8% |
| leaf_P_per_dry_mass | 0.7% | 17% | 18% |
| wood_density | 8.4% | not fitted | not fitted |
| leaf_N_per_area | 0.7% | not fitted | not fitted |
| leaf_P_per_area | 0.8% | not fitted | not fitted |

How the modelled and presence products compare, cells with both (modelled cells not extrapolated): correlation of log community means, with the median difference (modelled minus presence, log units) in brackets, from `comparison/community_layer_agreement.csv`. The last column is the first version of the modelled product, from the pooled model without a woodiness term:

| trait | woody | non-woody | all | all, pooled model |
|---|---|---|---|---|
| plant_height | 0.72 (−0.18) | 0.60 (−0.03) | 0.87 (−0.08) | 0.54 (0.24) |
| seed_dry_mass | 0.60 (−0.18) | 0.32 (1.24) | 0.53 (0.66) | 0.48 (0.99) |
| leaf_area | 0.86 (0.87) | 0.75 (−0.41) | 0.78 (0.08) | 0.78 (0.72) |
| leaf_mass_per_area | 0.77 (0.15) | −0.12 (−0.10) | 0.74 (−0.08) | 0.25 (0.24) |
| leaf_N_per_dry_mass | −0.30 (−0.15) | 0.42 (−0.17) | 0.11 (−0.18) | −0.21 (−0.21) |
| leaf_P_per_dry_mass | 0.02 (−0.15) | 0.04 (−0.34) | 0.27 (−0.25) | −0.27 (−0.38) |
| wood_density | 0.86 (−0.01) | | | 0.86 (0.01) |
| leaf_N_per_area | 0.92 (0.03) | | | 0.88 (0.10) |
| leaf_P_per_area | 0.68 (0.03) | | | 0.66 (0.04) |

With the woodiness term, the all-species modelled layer follows the presence layer much more closely for plant height and LMA, and loses most of its offset for leaf area and plant height: the pooled model's level was set by its mostly woody species. Within groups, agreement is good for woody plant height, leaf area, LMA, wood density and leaf N per area, and weak for leaf N and P per mass in either group, whose presence pattern climate does not predict. The model's level still differs from presence where its fitted species are a biased subset of a group: woody leaf area (+0.87) and non-woody seed mass (+1.24, from 27 species).

Observed products against each other, local records vs presence (cells with both): correlation 0.69–0.96 for woody species and 0.37–0.83 for non-woody where at least 9 cells compare.

### Native species only

Before the woodiness term was added, the gap between the pooled modelled product and the presence product was largest in the south-east. To test whether introduced species drove it, the presence and modelled products are also built from native species only (`observed_presence_native/`, `modelled_native/`): GBIF species native in at least one state or territory per the APC (`APCalign::native_anywhere_in_australia()`, APC release 2026-07-30), which drops 2,971 introduced species and 13 names not in the APC. Introduced species are 15% of species × cell presences, a median of 6.5% of the species in a cell (90th percentile 29%). The models are still fitted on all species with a woody or non-woody class.

With the pooled model, dropping them changed little: each native-only layer correlated with its version including introduced species at 0.88–1.00 (log scale), with median shifts of at most 0.06 log units, and the modelled-vs-presence gap stayed. Introduced species are a small share of the species in most cells, so they move a presence-only mean over all species little, even though their own distribution differs: Wenk et al. (2026, "Contrasting climate signals between native and introduced annual and perennial floras", bioRxiv, [doi:10.64898/2026.07.31.741935](https://doi.org/10.64898/2026.07.31.741935)) find introduced annuals concentrated in wet, populated areas with low precipitation seasonality, and native annual richness tracking climate far more tightly. The current comparisons, per group with the woodiness model, are in `comparison_native/` and `comparison_native_vs_all/`.

### Woody share

What the pooled model's gap did track is the share of woody species in the cell. Across cells, the presence community mean of all species follows the woody share of the species present (correlation 0.70 for plant height, 0.64 LMA, 0.60 seed mass, −0.44 leaf N per mass, −0.55 leaf P per mass); the pooled modelled one did not (−0.13 to 0.34), so the gap correlated with woody share at −0.80, −0.69, −0.56, 0.40 and 0.58 for those traits. That model was fitted mostly on woody species (81–99% of fitted species for the leaf traits, seed mass and wood density; 62% for plant height), while the cells' species are a median 44% woody. This is why the modelled product now has the woodiness term (above). Life form shifting with climate and with human influence is the pattern Wenk et al. (2026) report for annual versus perennial species.

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
- `modelled` depends on the model's two climate variables and the woody / non-woody split only, and on the species with enough records and a medoid to fit it; cells far from that climate space are flagged. Its non-woody layers rest on few species for seed mass (27), leaf P per mass (37) and leaf N per mass (42), with wide intervals and many cells extrapolated, and are absent for leaf N and P per area and wood density.
- Wood density for non-woody plants is nearly empty and of doubtful meaning; the layer is kept for completeness.
