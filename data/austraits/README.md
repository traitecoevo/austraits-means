# AusTraits inputs

The pipeline reads two AusTraits files from `data/austraits/`, both from the same release. They are git-ignored, so they are not in the repository; this README is.

The release in use is **AusTraits 8.0.0** ([10.5281/zenodo.24543275](https://doi.org/10.5281/zenodo.24543275), 26 November 2026; 503 datasets, 2,450,777 trait records, built from austraits.build `489360f0`). It is set in `_targets.R` by five values: the two paths (`austraits_rds`, `austraits_observed_means`), the release's version DOI (`austraits_doi`) and the MD5 of each file as its Zenodo record gives it (`austraits_rds_md5`, `austraits_observed_means_md5`). The pipeline stops if a file's MD5 differs (`check_file_md5()`), and writes the version, the austraits.build commit and the DOI into the product metadata.

| File | What it is | MD5 |
|---|---|---|
| `austraits-8.0.0.rds` | the database: trait records for every product, including the AVH height dataset `AVH_2026` | `7fb6c4bc0500247b3b0cf17b62414020` |
| `austraits-8.0.0-observed-means.parquet` | the observed means published with the release (`docs/observed-means-geometric-method.md`); read here, not computed | `597f9c53488172e603987da9e9fdecc1` |

To fetch them:

```sh
cd data/austraits
for f in austraits-8.0.0.rds austraits-8.0.0-observed-means.parquet; do
  curl -L -o "$f" "https://zenodo.org/api/records/24543275/files/$f/content"
done
md5 austraits-8.0.0.rds austraits-8.0.0-observed-means.parquet
```

(`traits.use::load_database()` and `austraits::load_austraits(version = "8.0.0", means = TRUE)` fetch the same files, but need traits.use 1.0.1 or austraits 4.0.1, which are newer than this project's `renv.lock` pins.)

## Moving to a new release

Fetch the new release's two files, then change the five values in `_targets.R`. `tar_make()` rebuilds everything from the records on, and writes the products under the new version (`outputs/species_means/<version>/`, `outputs/community/<version>/`, `outputs/forecast/<version>/`). Then update the release named above and in the main `README.md`.

## Using a development build instead

To build on austraits.build `develop` between releases: check out `develop` in a separate worktree (not a working checkout you curate in), install the traits.build version its `DESCRIPTION` requires, and run `source("build.R")`, then

```r
traits.build::build_post_export(austraits, "export/data/curr")
```

which writes the observed means beside the database. Copy `export/data/curr/austraits.rds` to `data/austraits/austraits-<version>-<short commit>.rds` and the observed means to `data/austraits/austraits-<version>-<short commit>-observed-means.parquet`, point the two paths in `_targets.R` at them, and set the DOI and the two MD5s to `NULL`. The products then carry the build's version string (such as `8.0.0.9000`) and its commit.

## Earlier files

Before 26 November 2026 the pipeline read development builds of what became 8.0.0: `austraits-7.0.0.9000-d535de9.rds` (austraits.build `d535de9e`, 19 October 2026; 470 datasets, 2,325,498 trait records) with its observed means, and before that `-62aa0e7` and `-cd26bb2`. The products made from them are under `outputs/*/7.0.0.9000/`. `austraits-7.0.0.rds`, the 7.0.0 release, was the input to the portal-method product, removed on 25 November 2026. None of these files is read now.
