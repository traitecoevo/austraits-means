# AusTraits inputs

The pipeline reads two AusTraits files from `data/austraits/` (git-ignored, so they are not in the repository).

## `austraits-7.0.0.rds`: the release

The AusTraits 7.0.0 release from Zenodo ([10.5281/zenodo.15718081](https://doi.org/10.5281/zenodo.15718081)), input to the portal-method observed means (species means 7.0.0). `read_austraits_release()` checks its MD5 against the Zenodo record (`c040dada5776177a1fabaced7ac0afdb`). To fetch it:

```r
austraits::load_austraits(version = "7.0.0", path = "data/austraits")
```

## `austraits-<version>-<commit>.rds`: a development build

An export of austraits.build `develop`, input to everything built on the development build (geometric observed means, modelled means, community layers). It carries the AVH height dataset `AVH_2026`, which replaced the separate AVH height CSV. The file in use is set by `austraits_dev_rds` in `_targets.R`; its version and commit are read from the build's `build_info` into the product metadata.

Current file: `austraits-7.0.0.9000-62aa0e7.rds`, built from austraits.build `develop` at `62aa0e79` (16 October 2026; 470 datasets, 2,325,498 trait records) with traits.build 2.2.0, as `develop` requires. One local adjustment was needed, reported in [traitecoevo/austraits.build#868](https://github.com/traitecoevo/austraits.build/issues/868): `build.R` was regenerated, because the committed file lists 410 of the 470 datasets, `AVH_2026` among those missing:

```r
traits.build::build_setup_pipeline(method = "base", database_name = "austraits")
```

The previous file, `austraits-7.0.0.9000-cd26bb2.rds` (`cd26bb21`, 442 datasets), was built the same way but with traits.build from its `add-date-parsing` branch, which some datasets then needed; that is no longer necessary.

To rebuild, check out austraits.build `develop` in a separate worktree (not a working checkout you curate in), regenerate `build.R` as above until #868 is resolved, run `source("build.R")`, and copy `export/data/curr/austraits.rds` to `data/austraits/austraits-<version>-<short commit>.rds`. Then point `austraits_dev_rds` at it.
