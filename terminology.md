# Terminology: species and community means

This project publishes two species-level trait estimates, plus community-level summaries of them mapped across space. This file fixes what they are called and what each one means. The paper, the exported data, the Zenodo and portal metadata, and the code should all use these terms. If a definition changes, change it here first.

Species terms decided 14 October 2026. Community terms proposed 14 October 2026; two choices still open (see "Open choices" below).

## The terms

| Role | Term | Data product (`estimate_type`) |
|---|---|---|
| Umbrella for both | **species mean** (species trait mean) | one table, with an `estimate_type` column |
| Simple | **observed mean** | `observed_mean` |
| Model-based | **modelled mean** | `modelled_mean` |

**Observed mean:** the geometric mean of a species' trait records.

**Modelled mean:** the model-estimated geometric mean trait value at the centre (medoid) of the species' climate distribution. This is the value from the within-species-centred (Mundlak) model, predicted at the medoid, where centred climate = 0. Where the species' climate centre needs to be named, "climate-centred" is the qualifier: "the modelled mean is the climate-centred estimate".

## Rules

- **Use "mean" for both, never "average".** Readers treat mean and average as synonyms, so the two words can't tell the two estimates apart. "Average" is also vague: it can mean the arithmetic mean, the median or the mode. It also clashes with the climate variable `average_temperature`. In lay text such as portal blurbs, "typical value" is acceptable, as long as the metadata gives the exact definition.
- **Avoid these alternative names:** "naive mean" sounds dismissive in a published product. "Raw mean" clashes with `value_type = raw`. "Sample mean" is ambiguous, because in trait work a "sample" is a physical specimen.
- **Both are geometric means.** Traits are treated as log-normal: values are averaged or modelled on the log scale, then back-transformed with `exp()`. For a log-normal, the geometric mean equals the median. It falls below the arithmetic mean, exp(μ + σ²/2), by an amount that grows with within-species spread. Keeping both estimates geometric means their difference comes from method alone, not from the type of summary. State "geometric" in definitions and metadata. It need not appear in the headline term.
- **If an arithmetic mean is ever published,** it goes in a separate column with its own clear label, never in place of the observed mean.

## What each estimate does and does not account for

**The observed mean** (`estimate_simple_means()` in `R/means/simple_means.R`) pools:

- height records from AVH specimens
- AusTraits records with `basis_of_record == "field"` and `value_type` in raw, mean or median

These records are not filtered by medoid, which gives maximum coverage. The mean is unweighted: a study with many records counts for more than one with few. A mean or median reported by a study counts the same as a single raw record. It does not adjust for where in the species' range the records came from.

**The modelled mean** corrects for records that cluster away from the centre of the species' climate range. Through partial pooling, it draws species with few records toward the across-species climate trend. A species is only included if it has a medoid that passes the reasonableness check and enough records to fit. All others get only an observed mean.

That last difference is why both are published. The observed mean gives coverage; the modelled mean corrects for where the records were collected.

## Point estimates and intervals

- **TMB:** the point estimate is `exp(fit)`. The interval is `exp(fit ± z·se)`, computed on the log scale.
- **BRMS:** the point estimate is `exp(posterior mean of log value)`. The interval is the exp of the posterior quantiles on the log scale. This is not the posterior mean of the trait on its original scale, so the methods section should say so.

## Community means

These summarise the species means of all species present at a location. They are computed for **woody** and **non-woody** plants separately, and then modelled across space.

| Term | Definition |
|---|---|
| **community mean** | the mean, on the log scale, of the species means of the species present at a location; back-transformed, so it is a geometric mean |
| **community variance** | the variance, on the log scale, of those species means |
| **woody** / **non-woody** | qualifiers on either term, for example "woody community mean" |
| **observed** / **modelled** | as for species: observed = computed directly from the species present in a spatial unit; modelled = the prediction of the spatial model |

Full forms follow one order: estimate type, then growth form, then level, then statistic. For example: "modelled woody community mean", "observed non-woody community variance".

Rules:

- **Always say "community".** It is what separates these from the species terms. "Modelled mean" on its own always refers to a species. "Modelled community mean" refers to the spatial surface.
- **Community variance is between-species only.** It is the spread of species means and leaves out variation within species. Say this wherever the term is defined. Report it in natural-log units. If the SD is also published, call it **community SD**.
- **Name the input species means in the metadata, not in the term.** Record whether the community statistics use the modelled means only, or the combined product (modelled where available, otherwise observed). Do not build that into the name.
- **Use "non-woody", not "herbaceous".** In AusTraits, "herbaceous" is one category of the `woodiness` trait, and non-woody is meant to cover everything that is not woody. The exact assignment of categories is an open choice below.
- **The spatial unit is part of the definition.** "Location" means a stated unit, such as a grid cell of a given resolution or a survey plot. Give it wherever the community terms are defined.

## Open choices for community terms

1. **Weighting. This decides the name.** If species are counted equally (presence only, as occurrence data such as AVH or GBIF allow), the terms above stand. If species are weighted by abundance or cover (for example from plot surveys), the field-standard terms apply instead: **community-weighted mean (CWM)** and **community-weighted variance (CWV)**. Do not call unweighted statistics CWM, because readers will assume abundance weighting.
2. **Woody / non-woody assignment.** Decide where intermediate AusTraits woodiness categories go (for example semi-woody, or woody only at the base), and how to handle species with conflicting or missing woodiness. Record the mapping here once decided.

## Code still to rename

The current outputs predate these terms:

- `combine_estimates()` writes a bare `estimate` column plus a `method` column with values `climate_centre` and `simple_mean`. The round-1 exports call it `estimate`, and the round-2 column is `geom_mean`.
- These should become `estimate_type` with values `modelled_mean` and `observed_mean`, matching the table above.

Not yet done.
