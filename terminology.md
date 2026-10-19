# Terminology: species and community means

This project publishes two species-level trait estimates, plus community-level summaries of them mapped across space. This file fixes what they are called and what each one means. The paper, the exported data, the Zenodo and portal metadata, and the code should all use these terms. If a definition changes, change it here first.

Species terms decided 14 October 2026. Community terms proposed 14 October 2026; their open choices were decided 19 October 2026 (see "Community choices" below). The choice of "mean" over "average" was confirmed on 14 October 2026: "average" had been suggested to leave room for medians and trimmed means, but both published estimates are geometric means, so "mean" is accurate, and flexibility of method is carried by the `method` column (below), not by the term.

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

**The observed mean** (geometric method: `estimate_observed_means_geometric()` in `R/means/observed_means_geometric.R`; full description in `docs/observed-means-geometric-method.md`) uses, from an AusTraits build:

- records with `value_type` raw, mean or median, whatever their basis of record or life stage;
- the AVH height dataset in austraits.build (`AVH_2026`), whose specimen label bounds (one minimum or maximum per specimen) are used as observations, like raw records;
- minimum and maximum records from preserved specimens and the literature (floras), as one geometric midpoint per observation.

Records are first summarised per dataset × location (the geometric mean of the records there) and per flora observation (the geometric midpoint of its bounds); the species value is the geometric mean of those groups, each counting once. So a study that measured many individuals at one site counts once for that site, and a mean or median reported by a study counts the same as a raw record within its group. Records with a zero or negative value are not used (a geometric mean needs positive values). Subspecies, varieties and forms are rolled up to species. The records are not filtered by medoid, which gives maximum coverage, and the mean does not adjust for where in the species' range the records came from. The arithmetic mean of the same groups is kept in its own column, `arithmetic_mean`.

The record selection differs from the modelled mean's, which uses field records only (raw, mean or median, plus the AVH bounds): the observed mean keeps the portal method's broader selection so that the two observed methods stay comparable.

**The modelled mean** corrects for records that cluster away from the centre of the species' climate range. Through partial pooling, it draws species with few records toward the across-species climate trend. A species is only included if it has a medoid that passes the reasonableness check and enough records to fit. All others get only an observed mean.

That last difference is why both are published. The observed mean gives coverage; the modelled mean corrects for where the records were collected.

## Point estimates and intervals

- **TMB:** the point estimate is `exp(fit)`. The interval is `exp(fit ± z·se)`, computed on the log scale.
- **BRMS:** the point estimate is `exp(posterior mean of log value)`. The interval is the exp of the posterior quantiles on the log scale. This is not the posterior mean of the trait on its original scale, so the methods section should say so.

## Product versions and methods

The species means product carries the **same version tag as the AusTraits release it derives from**, so the two are easily related: "species means 7.0.0" derives from AusTraits 7.0.0. A re-issue for the same release appends a fourth component (7.0.0.1). Products built from a development build carry that build's version string until the corresponding release exists.

Each row records how it was computed in a `method` column:

| `estimate_type` | `method` | Meaning |
|---|---|---|
| `observed_mean` | `portal_arithmetic` | The method the AusTraits data portal used for 7.0.0: an unweighted **arithmetic** mean of per-location means (raw, mean and median records) and per-observation flora midpoints (minimum and maximum records). Reproduced exactly from the portal code at the commit that built the deployed data, quirks included, so that product can be archived as what users have been seeing. |
| `observed_mean` | `geometric` | The observed mean as defined above (geometric), used from the next AusTraits build onward. The arithmetic mean is kept in its own clearly labelled column. |
| `modelled_mean` | `tmb` or `brms` | The modelled mean as defined above, by fitting backend. |

The rule that both estimates are geometric means therefore holds from the `geometric` method onward. The `portal_arithmetic` product is the one exception, and it is labelled as arithmetic wherever it appears.

## Community means

These summarise the trait values of the species present at a location. They are computed for all species, and for **woody** and **non-woody** plants separately.

| Term | Definition |
|---|---|
| **community mean** | the mean, on the log scale, of the trait values of the species present at a location, one value per species; back-transformed, so it is a geometric mean |
| **community variance** | the variance, on the log scale, of those species values |
| **woody** / **non-woody** | qualifiers on either term, for example "woody community mean" |
| **observed** / **modelled** | observed = computed directly from the species in a spatial unit, with no model; modelled = the trait-climate model's net average for the unit over the species present, averaging over between- and within-species effects |

An observed community mean takes each species' value either from its trait records collected in the unit (the *local records* product) or from its observed mean (the *GBIF presence* product). The modelled community mean uses no species values: the model's prediction for a species centred at climate m, observed at the unit's climate x, averaged over the species present, is b0 + b_within·x + (b_between − b_within)·m̄. It has no community variance. The three products are described in `docs/community-layers-design.md`.

Full forms follow one order: estimate type, then growth form, then level, then statistic. For example: "modelled woody community mean", "observed non-woody community variance".

Rules:

- **Always say "community".** It is what separates these from the species terms. "Modelled mean" on its own always refers to a species. "Modelled community mean" refers to the spatial surface.
- **Community variance is between-species only.** It is the spread of species means and leaves out variation within species. Say this wherever the term is defined. Report it in natural-log units. If the SD is also published, call it **community SD**.
- **Name the source of the species values in the metadata, not in the term.** Record whether an observed community statistic uses the species' records in the unit or their observed means. Do not build that into the name; the product name (`observed_local`, `observed_presence`) carries it.
- **Use "non-woody", not "herbaceous".** In AusTraits, "herbaceous" is one category of the `woodiness` trait, and non-woody is meant to cover everything that is not woody. The exact assignment of categories is an open choice below.
- **The spatial unit is part of the definition.** "Location" means a stated unit; for these layers, a 10 × 10 km cell of an equal-area (Lambert Azimuthal Equal Area) grid of Australia (decided 19 October 2026, #20). Give it wherever the community terms are defined.

## Community choices

1. **Weighting: decided 19 October 2026 (#18).** Presence only: every species present in a unit counts equally, as occurrence data such as AVH or GBIF allow, so the terms above stand: community mean and community variance. Abundance- or cover-weighted statistics from plot surveys would be **community-weighted mean (CWM)** and **community-weighted variance (CWV)**; do not use those names for these layers, because readers assume abundance weighting.
2. **Woody / non-woody assignment: decided 19 October 2026 (#19).** From the species' most likely `woodiness_detailed` value: `woody` and `woody_like_stem` (palms, cycads, grass trees) are woody; `herbaceous`, `semi_woody`, `woody_base`, `woody_like_inflorescence` and `woody_root` are non-woody. A species with no such value, or a tie between the two classes, is assigned from its growth forms (any herb, fern, tussock, geophyte, lycophyte, basal_large, graminoid or hummock form makes it non-woody; any other recorded form woody). Species still unassigned are left out of the woody and non-woody statistics and count only for all species; the coarse `woodiness` trait is not used.

## Column names in the exports

Every species means export carries `estimate_type` and `method` (values in "Product versions and methods"), and names its values the same way:

| Column | Meaning |
|---|---|
| `value` | the species mean (observed or modelled, as `estimate_type` says) in trait units |
| `log_value` | its natural log |
| `value_lower`, `value_upper` | the 94% interval of a modelled mean (and `log_value_lower`, `log_value_upper` on the log scale) |
| `arithmetic_mean` | the arithmetic mean, where an observed product gives it alongside the geometric `value` |
| `n_records` | records behind the value |

A bare `estimate`, `geom_mean`, or the older `method` values `climate_centre` and `simple_mean` are not used (renamed in #17). The portal-method product keeps the portal's own column names (`value_mean`, `value_geom_mean`, ...), because it reproduces the deployed portal output; its rows carry `estimate_type` and `method` like the others.
