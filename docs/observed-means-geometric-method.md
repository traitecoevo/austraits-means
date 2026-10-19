# Observed means, geometric method

## What this product is

The geometric observed means are per taxon × trait summaries of the trait records in an AusTraits build, with `estimate_type = observed_mean` and `method = geometric`. They succeed the portal method (`docs/observed-means-portal-method.md`) from the next AusTraits build onward. The product is built today from the austraits.build development build 7.0.0.9000 (commit `d535de9e`; see `data-raw/austraits/README.md`), and carries that build's version string until the corresponding release exists.

For a numeric trait the value is a **geometric mean**, as `terminology.md` requires of both species estimates. The arithmetic mean of the same groups is kept in its own column, `arithmetic_mean`.

The current build gives 601,959 rows: 153,994 numeric (one per taxon × trait × unit) and 447,965 categorical (one per taxon × trait), for 30,553 taxa and 93 traits (40 numeric, 53 categorical).

## Traits

The portal's 96 core traits (`data-raw/portal/trait_groups_for_portal.csv`), plus three traits the core list lacks: `leaf_N_per_area` and `leaf_P_per_area`, so that all nine modelled traits are covered, and `woodiness`, an input to the woody / non-woody lookup. Traits without records in the build drop out. The list is `geometric_mean_traits()` in `R/means/observed_means_geometric.R`. It is a starting point, not a decision about the eventual scope.

## Method

The recipe is the portal method's, on the log scale. `estimate_observed_means_geometric()` implements it.

**Record selection** (numeric traits).

- Records with `value_type` raw, mean or median are used, whatever their `basis_of_record` and `life_stage`, as in the portal method.
- Datasets listed in `dataset_value_type_rules()` use their own value types instead. Today that is the AVH height dataset `AVH_2026`, which carries one label bound per specimen (a minimum or a maximum) and no raw values; its bounds are used as observations, like raw records. The modelled means use the same rule.
- Minimum and maximum records from other datasets are used only when `basis_of_record` is `preserved_specimen` or `literature` (floras), as in the portal method.

**Groups.** Raw, mean and median records (and the AVH bounds) are grouped by taxon × trait × dataset × location; the group value is the geometric mean of its records. Flora bounds are grouped by observation; the group value is the geometric midpoint, sqrt(minimum × maximum), or the lone bound when only one is given.

**Species value.** `value` is the geometric mean of the group values, each group counting once whatever its number of records; `log_value` is its natural log. `arithmetic_mean` is the arithmetic mean of the groups' arithmetic means and midpoints, which for the same records equals the portal method's `value_mean`. `value_min` and `value_max` are the smallest and largest record values used.

**Zero and negative values.** A geometric mean needs positive values, so records with a zero or negative value are dropped before anything is computed, as for the modelled means; `arithmetic_mean` uses the same records. In the current build this drops 6,869 records: all 5,918 of `leaf_delta13C`, whose values are negative, so that trait has no geometric observed mean (it needs a transform first, as for the modelled means, #39); 895 of `petiole_length`, where zero means a sessile leaf; and a few of `leaf_stomatal_density_adaxial`, `seed_viability` and `leaf_inclination_angle`.

**Taxa.** Subspecies, varieties and forms are rolled up to their species (`roll_up_to_species()`), as for the modelled means, so the two species estimates join on species. Taxon details (rank, establishment means, APC identifier) come from the species' row in the build's taxa table. Genus-level names keep their own rows. In the current build 186,692 of 1,939,853 records were rolled up.

**Categorical traits.** One row per taxon × trait. Each record's value is split at spaces and every value counted; `value_count` lists `value (n)` for every value and `value_mode` gives the most frequent one, with tied values joined by a space.

**Reliability indicators** (#16; `R/means/reliability.R`). `n_records` is the number of records used (each bound counts as a record), `n_datasets` the number of datasets, `n_locations` the number of dataset × location groups, `n_flora` the number of flora minimum/maximum records, `n_field` the number of field records, and `sd_log` the standard deviation of the log record values. `reliability` is a provisional class: high with at least 3 datasets and 10 records, medium with at least 2 datasets or 5 records, low otherwise. The same columns are added to the portal-method product (after its own columns) and to the modelled means.

## Differences from the portal method

| | portal method | geometric method |
|---|---|---|
| species value | arithmetic mean | geometric mean, arithmetic mean alongside |
| `value_geom_mean` quirk | equals the mean | not carried; `value` is geometric |
| AVH height dataset | not in 7.0.0 | bounds used as observations |
| infraspecific taxa | kept as their own rows | rolled up to species |
| categorical traits | one row per distinct value's sources, no modal value | one row per taxon × trait, with `value_mode` |
| zero or negative values | used; `value_geom_mean` is NaN where the mean is negative | records dropped |
| record selection, grouping, group weighting | | unchanged |

The geometric mean sits below the arithmetic mean by an amount that grows with within-taxon spread. Across taxa the median ratio of `arithmetic_mean` to `value` is 1.10 for plant height and 1.12 for leaf length, and 1.00 for leaf mass per area and seed mass, where a third and a half of taxa respectively have a single record.

## Known limitations

- No life-stage or basis-of-record filter on raw, mean and median records, as in the portal method: seedlings, glasshouse and lab records are pooled with field records of adults.
- Groups count equally: a site with one record counts as much as one with a hundred, and a flora midpoint as much as a site mean.
- A lone flora bound stands in for the midpoint, biasing the value low (a minimum) or high (a maximum).
- Minimum and maximum records with a basis other than preserved specimen or literature are not used, except under a per-dataset rule.
- Zero values that carry meaning are dropped with the rest (a sessile leaf has a petiole length of zero), so the value describes the records with a positive value.
- The value does not adjust for where in the species' range the records came from; that is what the modelled mean is for.

## Columns

The column dictionary is written with the product (`observed_means_geometric_dictionary.csv`) and defined in `observed_means_geometric_dictionary()`.

## Outputs and rebuilding

`targets::tar_make()` builds `dev_records` → `observed_means_geometric` → `observed_means_geometric_files`, which writes `outputs/species_means/<build version>/observed_means_geometric` as `.parquet`, `.csv.gz`, `_dictionary.csv` and `_metadata.json`. The metadata records the build version and commit and the trait list.
