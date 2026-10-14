# AusTraits inputs

The pipeline reads two AusTraits files from `data/austraits/` (git-ignored, so they are not in the repository).

## `austraits-7.0.0.rds`: the release

The AusTraits 7.0.0 release from Zenodo ([10.5281/zenodo.15718081](https://doi.org/10.5281/zenodo.15718081)), input to the portal-method observed means (species means 7.0.0). `read_austraits_release()` checks its MD5 against the Zenodo record (`c040dada5776177a1fabaced7ac0afdb`). To fetch it:

```r
austraits::load_austraits(version = "7.0.0", path = "data/austraits")
```

## `austraits-<version>-<commit>.rds`: a development build

An export of austraits.build `develop`, input to everything built on the development build (geometric observed means, modelled means, community layers). It carries the AVH height dataset `AVH_2026`, which replaced the separate AVH height CSV. The file in use is set by `austraits_dev_rds` in `_targets.R`; its version and commit are read from the build's `build_info` into the product metadata.

Current file: `austraits-7.0.0.9000-cd26bb2.rds`, built from austraits.build `develop` at `cd26bb21` (442 datasets, 2,294,019 trait records). Building it needed two local adjustments, both reported in [traitecoevo/austraits.build#868](https://github.com/traitecoevo/austraits.build/issues/868):

- `build.R` was regenerated, because the committed file omits 32 datasets, `AVH_2026` among them:

  ```r
  traits.build::build_setup_pipeline(method = "base", database_name = "austraits")
  ```

- traits.build was installed from the branch `add-date-parsing` (`3f21320`), because eight datasets call `util_parse_collection_date()`, which traits.build `develop` does not yet have.

To rebuild, check out austraits.build `develop` in a separate worktree (not a working checkout you curate in), apply the same two steps until #868 is resolved, run `source("build.R")`, and copy `export/data/curr/austraits.rds` to `data/austraits/austraits-<version>-<short commit>.rds`. Then point `austraits_dev_rds` at it.
