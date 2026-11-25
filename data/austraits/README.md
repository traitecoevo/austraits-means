# AusTraits inputs

The pipeline reads three AusTraits files from `data/austraits/`. They are git-ignored, so they are not in the repository; this README is.

## `austraits-7.0.0.rds`: the release

The AusTraits 7.0.0 release from Zenodo ([10.5281/zenodo.15718081](https://doi.org/10.5281/zenodo.15718081)), input to the portal-method observed means (species means 7.0.0). `read_austraits_release()` checks its MD5 against the Zenodo record (`c040dada5776177a1fabaced7ac0afdb`). To fetch it:

```r
austraits::load_austraits(version = "7.0.0", path = "data/austraits")
```

## `austraits-<version>-<commit>.rds`: a development build

An export of austraits.build `develop`, input to everything built on the development build (geometric observed means, modelled means, community layers). It carries the AVH height dataset `AVH_2026`, which replaced the separate AVH height CSV. The file in use is set by `austraits_dev_rds` in `_targets.R`; its version and commit are read from the build's `build_info` into the product metadata.

Current file: `austraits-7.0.0.9000-d535de9.rds`, built from austraits.build `develop` at `d535de9e` (19 October 2026; 470 datasets, 2,325,498 trait records) with traits.build 2.2.0, as `develop` requires, from an unmodified checkout.

Earlier files: `austraits-7.0.0.9000-62aa0e7.rds` (`62aa0e79`) has the same data, built with a locally regenerated `build.R` before that fix reached `develop` ([traitecoevo/austraits.build#868](https://github.com/traitecoevo/austraits.build/issues/868)); `austraits-7.0.0.9000-cd26bb2.rds` (`cd26bb21`, 442 datasets) also needed traits.build from its `add-date-parsing` branch.

To rebuild, check out austraits.build `develop` in a separate worktree (not a working checkout you curate in), install the traits.build version its `DESCRIPTION` requires, run `source("build.R")`, and copy `export/data/curr/austraits.rds` to `data/austraits/austraits-<version>-<short commit>.rds`. Then point `austraits_dev_rds` at it, and make the observed means for the same build (next section).

## `austraits-<version>-<commit>-observed-means.parquet`: the observed means

The geometric observed means of that same build, as austraits.build publishes them (`R/observed_means.R` there, written by its `post_export()`; [austraits.build#899](https://github.com/traitecoevo/austraits.build/pull/899)). This repository reads them rather than computing its own (`read_observed_means()`), so the method has one home. The file in use is set by `austraits_dev_observed_means` in `_targets.R`, and must come from the same build as `austraits_dev_rds`.

For a release, take `austraits-<version>-observed-means.parquet` from the release files. For a development build, in the same austraits.build worktree, after `source("build.R")`:

```r
traits.build::build_post_export(austraits, "export/data/curr")
```

which writes the flattened exports and `austraits-<version>-observed-means.parquet` into `export/data/curr/`; copy the latter to `data/austraits/austraits-<version>-<short commit>-observed-means.parquet`. (To write only the observed means: `source("R/observed_means.R"); write_observed_means(observed_means(austraits), dir, austraits$build_info$version)`.)

Current file: `austraits-7.0.0.9000-d535de9-observed-means.parquet`, written on 24 November 2026 by austraits.build#899 (`2b3607c0`) from `austraits-7.0.0.9000-d535de9.rds`.
