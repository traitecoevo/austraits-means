# Observed means

## What they are and who makes them

The observed means are per taxon × trait summaries of the trait records in an AusTraits release, computed directly from the records, with no model (`estimate_type = observed_mean`). **austraits.build computes and publishes them** with each release from 8.0.0, as `austraits-<version>-observed-means.parquet` (also `.rds`, `.csv.gz` and `-dictionary.csv`) in the release's Zenodo record. This repository reads that file and does not compute or publish its own.

The method has one home, and this document does not restate it in full:

- the calculation is `traits.use::summarise_trait_means()`; its help page describes every step;
- what is specific to AusTraits is in austraits.build's `R/observed_means.R`: the rule that the AVH height dataset `AVH_2026` enters by its minimum and maximum records, and the files the release carries;
- the columns are defined in the dictionary published beside the file (`traits.use::trait_means_dictionary()`).

A change of method is made there and reaches this repository through a new release. The method began here as the "geometric method" (#15) and moved to austraits.build on 24 November 2026, then to traits.use (`roadmap.md`).

## The method in brief

- **Numeric traits** (`method = geometric`): records with `value_type` raw, mean or median are grouped by dataset × location, and minimum and maximum records from preserved specimens and the literature (floras) by observation. A group's value is the geometric mean of its records, or the geometric midpoint of a flora range. The species value is the geometric mean of the groups, each counting once. Records with a zero or negative value are not used. The arithmetic mean of the same groups is kept alongside in `arithmetic_mean`.
- **Signed traits** (`method = arithmetic`): a trait whose definition allows values below zero has no geometric mean, so its `value` is the arithmetic mean of the same groups, using every record.
- **Categorical traits** (`method = mode`): `value_mode` is the most frequent value, with ties joined by a space, and `value_count` lists every value with its count. Traits that record months are summarised month by month (`method = majority`).
- **Taxa**: subspecies, varieties and forms are rolled up to their species. Hybrid formulae and genus-level names keep their own rows.
- **Evidence**: `n_records`, `n_datasets`, `n_locations`, `n_flora`, `n_field` and `sd_log` count what is behind each value, and `reliability` is a provisional class from them.

In AusTraits 8.0.0 the file has 767,333 rows for 34,115 taxa and 628 traits: 432 geometric, 43 arithmetic, 149 by mode and 4 by majority.

## What this repository relies on

- `read_observed_means()` (`R/means/observed_means_release.R`), behind the `observed_means_geometric` target, stops if the file lacks a column the pipeline reads, lacks one of the nine modelled traits, `woodiness_detailed` or `plant_growth_form`, or gives an arithmetic mean for a modelled trait. The file must come from the same release as the trait records (`austraits_rds` and `austraits_observed_means` in `_targets.R`; `data/austraits/README.md`).
- The modelled traits are fitted on the log scale, so their observed means must be geometric for the two species estimates to be comparable.
- `roll_up_to_species()` here must keep rolling taxa up the same way as `summarise_trait_means()`, so the modelled and observed means join on species.
- The woody / non-woody lookup reads `value_mode` of `woodiness_detailed` and `plant_growth_form` (`R/community/woody_lookup.R`).
- The observed mean uses a broader set of records than the modelled mean, which uses field records only (raw, mean or median, plus the AVH bounds) of species with a climate centre (`terminology.md`).

## Known limitations

- There is no life-stage or basis-of-record filter on raw, mean and median records: seedlings, glasshouse and lab records are pooled with field records of adults.
- Groups count equally: a site with one record counts as much as one with a hundred, and a flora midpoint as much as a site mean.
- A lone flora bound stands in for the midpoint, biasing the value low (a minimum) or high (a maximum).
- Zero values that carry meaning are dropped from a geometric mean (a sessile leaf has a petiole length of zero), so the value describes the records with a positive value.
- The value does not adjust for where in the species' range the records came from; that is what the modelled mean is for.
