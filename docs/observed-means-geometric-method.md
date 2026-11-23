# Observed means, geometric method

## What this product is

The geometric observed means are per taxon × trait summaries of the trait records in an AusTraits build, with `estimate_type = observed_mean` and `method = geometric`. They succeed the portal method (`docs/observed-means-portal-method.md`) from the next AusTraits build onward. **austraits.build computes and publishes them**, from the release after 7.0.0 onward, as `austraits-<version>-observed-means.parquet` (`.csv.gz`, `-dictionary.csv`) beside each release (`R/observed_means.R` there, written by its `post_export()`; [austraits.build#899](https://github.com/traitecoevo/austraits.build/pull/899)). The method moved there from this repository, which now only reads the file; this document describes the method, and austraits.build's code is the authority if the two disagree. The pipeline here reads the file published for the development build 7.0.0.9000 (commit `d535de9e`; see `data-raw/austraits/README.md`).

For a numeric trait the value is a **geometric mean**, as `terminology.md` requires of both species estimates. The arithmetic mean of the same groups is kept in its own column, `arithmetic_mean`.

## Traits

The portal's core traits (`data-raw/portal/trait_groups_for_portal.csv`), the traits the AusTraits API summarises, and three more: `leaf_N_per_area` and `leaf_P_per_area`, so that all nine modelled traits are covered, and `woodiness`, an input to the woody / non-woody lookup. Traits without records in the build drop out. The list is `observed_means_traits()` in austraits.build's `R/observed_means.R`.

## Method

The recipe is the portal method's, on the log scale. `observed_means()` in austraits.build's `R/observed_means.R` implements it.

**Record selection** (numeric traits).

- Records with `value_type` raw, mean or median are used, whatever their `basis_of_record` and `life_stage`, as in the portal method.
- Datasets listed in `dataset_value_type_rules()` use their own value types instead. Today that is the AVH height dataset `AVH_2026`, which carries one label bound per specimen (a minimum or a maximum) and no raw values; its bounds are used as observations, like raw records. The modelled means use the same rule.
- Minimum and maximum records from other datasets are used only when `basis_of_record` is `preserved_specimen` or `literature` (floras), as in the portal method.

**Groups.** Raw, mean and median records (and the AVH bounds) are grouped by taxon × trait × dataset × location; the group value is the geometric mean of its records. Flora bounds are grouped by observation; the group value is the geometric midpoint, sqrt(minimum × maximum), or the lone bound when only one is given.

**Species value.** `value` is the geometric mean of the group values, each group counting once whatever its number of records; `log_value` is its natural log. `arithmetic_mean` is the arithmetic mean of the groups' arithmetic means and midpoints, which for the same records equals the portal method's `value_mean`. `value_min` and `value_max` are the smallest and largest record values used.

**Zero and negative values.** A geometric mean needs positive values, so records with a zero or negative value are dropped before anything is computed, as for the modelled means; `arithmetic_mean` uses the same records. In the current build this drops 6,869 records: all 5,918 of `leaf_delta13C`, whose values are negative, so that trait has no geometric observed mean (it needs a transform first, as for the modelled means, #39); 895 of `petiole_length`, where zero means a sessile leaf; and a few of `leaf_stomatal_density_adaxial`, `seed_viability` and `leaf_inclination_angle`.

**Taxa.** Subspecies, varieties and forms are rolled up to their species (the name before its rank marker), as for the modelled means (`roll_up_to_species()` here), so the two species estimates join on species. Hybrid formulae (*A x B var. c*) are not rolled up and keep their own rows. Taxon details (rank, establishment means, APC identifier) come from the species' row in the build's taxa table. Genus-level names keep their own rows. In the current build 186,692 of 1,939,853 records were rolled up.

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

The column dictionary is published with the product (`austraits-<version>-observed-means-dictionary.csv`) and defined in `observed_means_dictionary()` in austraits.build. It has the columns above except `source_primary_key`: a dataset's references are in the release's `sources` table under its `dataset_id`.

## Outputs and rebuilding

austraits.build writes the files when it exports a build (`create_release()`, or `traits.build::build_post_export()` on a development build; `data-raw/austraits/README.md`). Here, the `observed_means_geometric` target reads them with `read_observed_means()` (`R/means/observed_means_release.R`), which checks the columns the pipeline uses and that every trait it needs is present; this repository no longer writes its own copy.
