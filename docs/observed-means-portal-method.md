# Observed means, portal method: AusTraits species means 7.0.0

## What this product is

AusTraits species means 7.0.0 is a table summarising the trait records of each taxon and trait in AusTraits 7.0.0. It reproduces the species means that the AusTraits data portal computes offline for that release, so that what portal users have been seeing can be archived and cited. Every row carries `estimate_type = observed_mean` and `method = portal_arithmetic`.

**This is an arithmetic mean.** For a numeric trait it is the unweighted arithmetic mean of per-location means and per-observation flora midpoints, not of the individual records. It is not a geometric mean, and the column `value_geom_mean` is not a geometric mean either: it equals `value_mean`, or is `NaN` where that is negative (see [Known limitations](#known-limitations)). These quirks are the product. They are documented here and not fixed. Fixes belong to the geometric observed mean, a later product with `method = geometric`.

The product carries the same version tag as the AusTraits release it derives from, so "species means 7.0.0" derives from AusTraits 7.0.0. A re-issue for the same release appends a fourth component, such as 7.0.0.1.

The product in numbers:

- 697,079 rows: 146,211 numeric rows (one per taxon × trait × unit) and 550,868 categorical rows.
- 33,287 taxa, of which 27,115 have at least one numeric trait.
- 92 traits: 39 numeric and 53 categorical.
- 1,479,903 core-trait records go in.
- Of the numeric rows, 102,390 draw on flora midpoints only, 38,375 on location means only, and 5,446 on both.

## Inputs

The product is built from three inputs.

- **The AusTraits 7.0.0 release**, Zenodo DOI [10.5281/zenodo.15718081](https://doi.org/10.5281/zenodo.15718081), file `austraits-7.0.0.rds`. The pipeline checks its MD5 (`c040dada5776177a1fabaced7ac0afdb`) before use and stops if it differs.
- **The flattened database.** The release is flattened with `austraits::flatten_database()` (austraits R package 3.1.1), which joins traits, locations, contexts, methods, taxa and sources into one table of records.
- **The core-trait list.** The portal computes means only for its 96 core traits, the rows of its `trait_groups_for_portal.csv` where `core_trait` is not `NA`. Four of them have no records in 7.0.0 (`leaf_hairs_juvenile_leaves`, `leaf_axil_angle`, `stem_density` and `climbing_plant_height`), so the product covers 92 traits. The file was deleted from the portal later. This repository holds the copy from immediately before its deletion (portal commit `dc8a965^`), which is identical to the file at commit `13479b9`. It is kept in `data-raw/portal/` and must not be edited.

The records are prepared in this order, reproducing the portal's data build of 9 February 2026:

1. Read the release and check its MD5.
2. Set `source_primary_key` to `dataset_id` for the Bryant_2021_2 and Bryant_2021_3 datasets. In 7.0.0 the key for these datasets is wrong, and the portal's data build applied this fix before flattening. The code matches any `dataset_id` that contains `Bryant_2021`.
3. Flatten with `austraits::flatten_database()`.
4. Keep the records of the core traits (`trait_name %in% core_traits`).

The portal build also added a `row_id` to the records and re-encoded `measurement_remarks`. Neither reaches the means, so both steps are left out.

The method reads these columns of the flattened records: `taxon_name`, `trait_name`, `unit`, `value`, `value_type`, `basis_of_record`, `dataset_id`, `location_id`, `observation_id`, `original_name`, `source_primary_key`, `source_primary_citation`, and the taxon columns from `taxon_rank` to `scientific_name_id`. Latitude, longitude and location name are carried through the location step but dropped at the species step, so they are not in the output.

## Method

The method is a port of `estimate_species_trait_means()` and its helpers in `R/species_means.R` of `traitecoevo/austraits.portal` at commit `13479b9`, as run by `prepare_data_for_portal()`. The logic is copied as it was. Only namespaces were made explicit and functions renamed.

A trait is numeric if it has a `unit` and categorical if `unit` is `NA`. The two kinds are summarised differently and their rows are stacked in one table.

### Numeric traits

Numeric traits have one row per taxon × trait × unit. The row is built in five steps.

**Step 1: record selection.** Two groups of records are used. Records of any other kind are ignored.

- Records with `value_type` raw, mean or median are used. They are used whatever their `basis_of_record` (field, lab, preserved_specimen, literature, captive_cultivated, experiment and so on) and whatever their `life_stage` (adult, seedling, sapling and so on). There is no filter on either.
- Records with `value_type` minimum or maximum are used only when `basis_of_record` is exactly `preserved_specimen` or `literature`. Minimum and maximum records with any other basis are not used at all. For example, 11,820 minimum or maximum records with basis `field` are dropped.

Records whose `value` cannot be read as a number are dropped. A `value_type` other than raw, mean, median, minimum or maximum (for example mode) is not used.

**Step 2: location means.** The raw, mean and median records are grouped by taxon × trait × dataset × location (`taxon_name`, `trait_name`, `dataset_id`, `location_id`, `unit`). Each group becomes one location mean:

- its value is the arithmetic mean of the records in the group, with raw, mean and median records counted alike;
- its minimum, maximum and median are those of the records in the group;
- its record count is the number of records in the group.

Records with no `location_id` are grouped together as one location per taxon × trait × dataset.

**Step 3: flora midpoints.** The selected minimum and maximum records are grouped by taxon × trait × dataset × observation × `original_name`. Each group becomes one flora midpoint:

- its value is the arithmetic mean of the records in the group, normally one minimum and one maximum, so the midpoint of the range;
- its minimum and maximum are the smallest and largest record in the group;
- its median is the midpoint itself;
- it counts as one flora observation, whether it has one bound or two.

A group with only a minimum, or only a maximum, is accepted. The lone bound is its own midpoint. If a group has more than one minimum or maximum record, all of them enter the mean.

**Step 4: species-level summary.** The location means and the flora midpoints are stacked and grouped by taxon × trait × unit. Within each group:

- `value_mean` is the arithmetic mean of the location means and flora midpoints. Every location mean and every flora midpoint has the same weight, whatever the number of records behind it.
- `value_min` is the smallest of the group minimums, and `value_max` the largest of the group maximums. They are the smallest and largest record values used, including flora bounds.
- `value_median` is the median of the location medians and flora midpoints. It is a median of group medians, not the median of all records.
- `value_geom_mean` equals `value_mean`. The code computes it as `10^mean(log10(value_mean))` in the same `summarise()` that has just overwritten `value_mean` with the species mean, so the logarithm and its inverse cancel. Where `value_mean` is negative the logarithm is undefined and `value_geom_mean` is `NaN`. The location-level geometric means calculated in step 2 are discarded.

**Step 5: replicate counts and sources.**

- `location_replicates` is the number of location means used, and `flora_replicates` the number of flora midpoints used.
- `all_replicates` is the number of records used. Each raw, mean or median record counts once. Each flora observation counts once, however many bound records it has. A mean or median record counts as one record whatever the sample size behind it.
- `dataset_id`, `observation_id` and `source_primary_key` list the distinct contributing datasets, observations and primary source keys, joined by `"; "`. Observation identifiers are unique only within a dataset.
- `source_primary_citation` is built last. Each key in `source_primary_key` is looked up in the flattened records and its primary citation is taken. The citations are joined by `"; "` in the order of the keys. Only primary sources are listed.

### Categorical traits

Categorical traits are counted; no mean is taken. Every record of the trait is used, whatever its `value_type`, `basis_of_record` or `life_stage`.

1. The `value` of each record is split at spaces. A record listing several space-separated values counts once for each value.
2. For each taxon × trait, the records are counted per distinct value.
3. `value_count` lists every value with its count, as `value (n)`, joined by `"; "`, in sorted order of the value (byte order, so upper case sorts before lower case). For example, `herbaceous (2); woody (3)`.
4. Each distinct value of a taxon × trait first gets its own row, and every row carries the same `value_count`. Its `all_replicates` is the count for its own value. Its `dataset_id`, `observation_id` and `source_primary_key` list the sources of its own value.

There is no separate column for the modal value. It is the value with the largest count in `value_count`. Ties are not resolved.

The value itself is dropped from each row. Rows that are identical in every column are then collapsed to one, so two values with the same count and the same sources share a row. `value_count` always lists all the values, so use it, not the per-row columns, to read what was recorded. This is why the product has more categorical rows than taxon × trait pairs, and why the per-row source columns cannot always be traced back to a single value.

### Worked example

The numbers below are made up. They show one taxon and one numeric trait, `leaf_mass_per_area` in g/m2.

| Dataset | Location or observation | `value_type` | `basis_of_record` | Values | Used? |
|---|---|---|---|---|---|
| A | location L1 | raw | field | 100, 120, 140 | yes |
| A | location L2 | raw | field | 80 | yes |
| B | location L1 | mean, median | lab | 150, 130 | yes |
| C | observation O1 | minimum, maximum | preserved_specimen | 60, 200 | yes |
| C | observation O2 | maximum | literature | 90 | yes |
| D | observation O3 | minimum | field | 50 | no: wrong basis |

Step 2 gives three location means and step 3 gives two flora midpoints:

| Group | Kind | Value | Minimum | Maximum | Median | Records |
|---|---|---|---|---|---|---|
| A, L1 | location | 120 | 100 | 140 | 120 | 3 |
| A, L2 | location | 80 | 80 | 80 | 80 | 1 |
| B, L1 | location | 140 | 130 | 150 | 140 | 2 |
| C, O1 | flora | 130 | 60 | 200 | 130 | 1 |
| C, O2 | flora | 90 | 90 | 90 | 90 | 1 |

Step 4 gives the species row:

- `value_mean` = (120 + 80 + 140 + 130 + 90) / 5 = 112.
- `value_min` = 60 and `value_max` = 200.
- `value_median` = median(120, 80, 140, 130, 90) = 120.
- `value_geom_mean` = 112, the same as `value_mean`.
- `all_replicates` = 3 + 1 + 2 + 1 + 1 = 8, `location_replicates` = 3 and `flora_replicates` = 2.

Dataset D's field minimum is not used. The lone maximum of 90 from dataset C is used as its own midpoint; had its minimum been recorded, the midpoint would have been lower, so a lone maximum biases the mean up.

For a categorical trait, suppose a taxon has four records of `woodiness`: `woody`, `woody`, `herbaceous` and `woody herbaceous`. The last record counts once for each of its two values. The counts are woody 3 and herbaceous 2, so `value_count` is `herbaceous (2); woody (3)`. There are two rows, one with `all_replicates` = 2 and one with `all_replicates` = 3.

## Known limitations

These follow from the method and apply to every row. Read them before using the values.

- **Arithmetic mean, not geometric.** Trait values are often right-skewed. For such traits an arithmetic mean sits above the geometric mean and the median.
- **Unweighted by record.** A location with one record counts as much as a location with 100 records. A dataset with many locations counts more than a dataset with few. A flora midpoint counts as much as a location mean.
- **`value_geom_mean` equals `value_mean`.** It carries no information beyond `value_mean`. In the code it is `10^mean(log10(value_mean))`. For a positive `value_mean` this returns `value_mean` up to floating-point rounding. For a negative `value_mean`, `log10()` gives `NaN`, so the code gives `NaN` and not `value_mean`. This affects 1,671 rows, all of `leaf_delta13C`, whose values are negative.
- **No life-stage or basis-of-record filter on raw, mean and median records.** Seedlings, saplings, lab, glasshouse, cultivated and experimental records are mixed with field records of adults. Use the `observation_id` and `dataset_id` columns with the AusTraits release if you need to filter.
- **Lone flora bounds are accepted.** A minimum without a maximum, or the reverse, is used as the midpoint. A lone minimum biases the value low and a lone maximum biases it high.
- **Minimum and maximum records are used only from preserved specimens and the literature.** Those with any other basis are dropped, including 11,820 field records.
- **Group medians and replicate counts.** `value_median` is a median of group medians. `all_replicates` counts records, so a mean from 50 individuals counts as one record.
- **No roll-up of infraspecific or genus-level names.** Each `taxon_name` is summarised on its own. Records of a subspecies or variety do not contribute to its species row, and records identified only to genus have their own rows.
- **Core traits only.** The product covers the portal's 92 core traits that have records in 7.0.0, not every trait in AusTraits.
- **Categorical traits.** There is no modal-value column, rows do not name their value, and every record of the trait is used whatever its `value_type`, basis or life stage.

## Column dictionary

Missing values are `NA` in the parquet file and empty in the CSV. "Numeric" and "categorical" in the second column give the rows where the column is filled. The `_dictionary.csv` file gives the same definitions and adds the R class of each column.

| Column | Filled for | Definition |
|---|---|---|
| `taxon_name` | all | Taxon name as aligned in the AusTraits release (APC-accepted where possible). Infraspecific and genus-level names are kept as they are, not rolled up to species. |
| `trait_name` | all | AusTraits trait name (APD). |
| `unit` | numeric | Unit of the numeric trait values. Empty for categorical traits. |
| `estimate_type` | all | Always `observed_mean`: computed directly from the records. |
| `method` | all | Always `portal_arithmetic`: the AusTraits data portal's method for release 7.0.0, an arithmetic mean. |
| `value_count` | categorical | Each value with the number of records carrying it, as `value (n)`, joined by `"; "`. A record listing several space-separated values counts once for each. |
| `value_mean` | numeric | Arithmetic mean of the taxon's location means and flora midpoints. A location mean is the arithmetic mean of the raw, mean and median records of one dataset at one location. A flora midpoint is the mean of the minimum and maximum records of one observation from a preserved specimen or the literature. |
| `value_min` | numeric | Smallest record value across the location records and flora bounds used. |
| `value_max` | numeric | Largest record value across the location records and flora bounds used. |
| `value_median` | numeric | Median of the location medians and flora midpoints. A median of group medians, not the median of all records. |
| `value_geom_mean` | numeric | Equal to `value_mean` where that is positive, and `NaN` where it is negative (1,671 `leaf_delta13C` rows). Despite its name this is not a geometric mean: the portal method computes it from the species mean that has already been calculated. Kept for fidelity to the deployed portal output. |
| `all_replicates` | all | Number of records used. Each raw, mean or median record counts once, and each flora observation counts once whether it has one bound or two. For categorical traits, the number of records carrying the value of this row. |
| `location_replicates` | numeric | Number of dataset × location groups contributing a location mean. |
| `flora_replicates` | numeric | Number of flora observations contributing a midpoint. |
| `dataset_id` | all | AusTraits datasets contributing, joined by `"; "`. |
| `observation_id` | all | AusTraits observation identifiers contributing, joined by `"; "`. Identifiers are unique only within a dataset. |
| `source_primary_key` | all | Primary source keys of the contributing datasets, joined by `"; "`. |
| `value_type` | categorical | `value_count` for categorical rows. Empty for numeric rows. |
| `type` | all | `numerical` or `categorical`. |
| `taxon_rank` | all | Taxonomic rank of `taxon_name`, from the AusTraits taxa table. |
| `taxonomic_status` | all | Taxonomic status of `taxon_name` in the taxonomic dataset. |
| `taxonomic_dataset` | all | Taxonomic reference the name was aligned to (APC or APNI). |
| `genus` | all | Genus, from the AusTraits taxa table. |
| `family` | all | Family, from the AusTraits taxa table. |
| `binomial` | all | Binomial of `taxon_name`, from the AusTraits taxa table. |
| `trinomial` | all | Trinomial of `taxon_name`, where applicable. |
| `taxon_distribution` | all | Distribution by state or territory as given by APC. |
| `establishment_means` | all | Native, naturalised or both, as given by APC. |
| `scientific_name` | all | Scientific name with authority. |
| `taxon_id` | all | APC taxon identifier. |
| `taxon_id_genus` | all | APC identifier of the genus. |
| `taxon_id_family` | all | APC identifier of the family. |
| `scientific_name_id` | all | APNI scientific name identifier. |
| `row_id` | all | Row number in this table. It depends on row order and is not a stable identifier. |
| `source_primary_citation` | all | Citations for `source_primary_key`, joined by `"; "`. |

## Reproduction

The product is tested against the species means that the portal deployed for AusTraits 7.0.0, the file `austraits-species-averages.parquet`, built on 9 February 2026. The method is pinned to portal commit `13479b9`, the commit that built the deployed data.

The comparison is made by `compare_portal_reproduction()`:

- Numeric rows are compared by key (`taxon_name`, `trait_name`, `unit`). Numbers must agree within a relative tolerance of 1e-8 and text must agree exactly.
- Categorical rows have no unique key, because the portal keeps one row per value's sources. They are compared as multisets of whole rows.
- `row_id` and row order are left out of the comparison, because the portal had no lockfile and so the package versions of its build are not known. `estimate_type` and `method` are left out because the deployed file does not have them.

The result is a match in every row and column:

| Part | Check | Result |
|---|---|---|
| all | rows in the product and in the deployed file | 697,079 and 697,079 |
| all | columns on one side only | 0 |
| numeric | rows in the product and in the deployed file | 146,211 and 146,211 |
| numeric | duplicate keys on either side | 0 |
| numeric | keys on one side only | 0 |
| numeric | keys with any value mismatch | 0 |
| categorical | rows in the product and in the deployed file | 550,868 and 550,868 |
| categorical | rows on one side only | 0 |
| categorical | taxon × trait pairs with any mismatch | 0 |

The match is in fact stronger than the test requires. With dplyr 1.2.1 the pipeline's output is identical to the deployed file, row order and `row_id` included.

The check runs as the `portal_reproduction` target and is also asserted by `tests/testthat/test-observed-means-portal.R`. The results are written to `outputs/checks/` as `portal_reproduction_summary.csv`, `portal_reproduction_columns.csv` and `portal_reproduction_examples.csv`.

## Relation to the other species means

The project publishes two species-level estimates, defined in `terminology.md`, which fixes the terms used here.

- **Observed mean, portal method** (this product): `estimate_type = observed_mean`, `method = portal_arithmetic`. An arithmetic mean, labelled as such.
- **Observed mean, geometric method**: `estimate_type = observed_mean`, `method = geometric`. The successor, used from the next AusTraits build onward. It is a geometric mean, it keeps the arithmetic mean in its own clearly labelled column, and it is where the limitations above are addressed.
- **Modelled mean**: `estimate_type = modelled_mean`, `method = tmb` or `brms`. The model-estimated geometric mean at the centre of the species' climate distribution. Only species with a usable climate centre and enough records get one.

`terminology.md` states that both species estimates are geometric means. The portal-method product is the one exception, and it is labelled as arithmetic wherever it appears. Do not pool or compare rows of different `method` values without checking what each one is. Filter on `method`.

## How to rebuild it

Run the pipeline from the repository root:

```r
targets::tar_make()
```

The targets run in this order:

1. `portal_records`: read the AusTraits 7.0.0 release, check its MD5, apply the Bryant_2021 fix, flatten it and keep the core-trait records.
2. `observed_means_portal`: compute the means with `estimate_observed_means_portal()`.
3. `observed_means_portal_files`: label the rows and write the output files.

Two input files are git-ignored and must be in place first: the release at `data/austraits/austraits-7.0.0.rds` (how to fetch it: `data-raw/austraits/README.md`), and, for the reproduction check, the deployed portal file at `data/portal/austraits-7.0.0-full/austraits-species-averages.parquet`. The `portal_reproduction` target compares the result with that file.

The output goes to `outputs/species_means/7.0.0/`:

- `observed_means_portal_arithmetic.parquet`: the table.
- `observed_means_portal_arithmetic.csv.gz`: the same table as gzipped CSV, with missing values empty.
- `observed_means_portal_arithmetic_dictionary.csv`: the column dictionary.
- `observed_means_portal_arithmetic_metadata.json`: the product, version, method labels, source release DOI, method source (`traitecoevo/austraits.portal@13479b9, R/species_means.R`), trait list, this document, the repository and its commit, the number of rows and the date written.

The code that implements the method is in `R/` of the repository `https://github.com/traitecoevo/austraits-means`, in `observed_means_portal.R`.
