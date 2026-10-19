# Baseline: modelled means, June 2026

A snapshot of the June 2026 TMB outputs, taken from the targets store before the modelled means were refitted on the austraits.build development build (#14). It is kept in git because `outputs/` is git-ignored and the targets store is overwritten by a refit.

Data basis of these outputs:

- traits other than height: AusTraits 7.0.0 release, field records with value type raw, mean or median;
- plant height: the AVH height CSV (extract of 11 March 2026), with one `dataset_id` per AVH institution, names effectively at species rank;
- occurrences and medoids: the cleaned GBIF extract (download `0001494-250914085247600`) and the medoids built from it in February 2026;
- name crosswalk: built in June 2026 from the AVH CSV names, with APCalign's then-current APC release (not pinned).

Files:

- `trait_model_summary.csv`: the pipeline's fixed-effect summary, one row per trait, as written to `outputs/traits/`.
- `brms_diagnostics.csv`: BRMS convergence for the eight traits other than height, as written to `outputs/traits/`.
- `fixed_effects.csv`: the fixed effects of each TMB fit with standard error and 95% interval (from the `tmb_fits` target).
- `trait_counts.csv`: records and species per trait in `unified` and in the fits (species with more than 5 records), and the number of species with an observed mean.
- `species_centres.csv.gz`: per trait × species, the modelled mean on the log scale with its 94% interval, the number of records, and the centre climate.

`compare_with_baseline()` (`R/process/compare_baseline.R`) reads these files; the comparison is reported in `reports/analysis.qmd`. Do not edit them.
